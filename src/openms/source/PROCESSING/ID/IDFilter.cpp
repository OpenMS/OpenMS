// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <regex>

using namespace std;

namespace OpenMS
{

  struct IDFilter::HasMinPeptideLength {
    typedef PeptideHit argument_type; // for use as a predicate

    Size length_;

    explicit HasMinPeptideLength(Size length) : length_(length)
    {
    }

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getSequence().size() >= length_;
    }
  };


  struct IDFilter::HasMinCharge {
    typedef PeptideHit argument_type; // for use as a predicate

    Int charge_;

    explicit HasMinCharge(Int charge) : charge_(charge)
    {
    }

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getCharge() >= charge_;
    }
  };


  struct IDFilter::HasLowMZError {
    typedef PeptideHit argument_type; // for use as a predicate

    double precursor_mz_, tolerance_;

    HasLowMZError(double precursor_mz, double tolerance, bool unit_ppm) : precursor_mz_(precursor_mz), tolerance_(tolerance)
    {
      if (unit_ppm)
        this->tolerance_ *= precursor_mz / 1.0e6;
    }

    bool operator()(const PeptideHit& hit) const
    {
      Int z = hit.getCharge();
      if (z == 0)
        z = 1;
      double peptide_mz = hit.getSequence().getMZ(z);
      return fabs(precursor_mz_ - peptide_mz) <= tolerance_;
    }
  };


  struct IDFilter::HasMatchingModification {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& mods_;

    explicit HasMatchingModification(const set<String>& mods) : mods_(mods)
    {
    }

    bool operator()(const PeptideHit& hit) const
    {
      const AASequence& seq = hit.getSequence();
      if (mods_.empty())
      {
        return seq.isModified();
      }
      for (Size i = 0; i < seq.size(); ++i)
      {
        if (seq[i].isModified())
        {
          String mod_name = seq[i].getModification()->getFullId();
          if (mods_.count(mod_name) > 0)
            return true;
        }
      }

      // terminal modifications:
      if (seq.hasNTerminalModification())
      {
        String mod_name = seq.getNTerminalModification()->getFullId();
        if (mods_.count(mod_name) > 0)
          return true;
      }
      if (seq.hasCTerminalModification())
      {
        String mod_name = seq.getCTerminalModification()->getFullId();
        if (mods_.count(mod_name) > 0)
          return true;
      }

      return false;
    }
  };


  struct IDFilter::HasMatchingSequence {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& sequences_;
    bool ignore_mods_;

    explicit HasMatchingSequence(const set<String>& sequences, bool ignore_mods = false) : sequences_(sequences), ignore_mods_(ignore_mods)
    {
    }

    bool operator()(const PeptideHit& hit) const
    {
      const String& query = (ignore_mods_ ? hit.getSequence().toUnmodifiedString() : hit.getSequence().toString());
      return (sequences_.count(query) > 0);
    }
  };


  struct IDFilter::HasNoEvidence {
    typedef PeptideHit argument_type; // for use as a predicate

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getPeptideEvidences().empty();
    }
  };

  struct IDFilter::HasRTInRange {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double rt_min_, rt_max_;

    HasRTInRange(double rt_min, double rt_max) : rt_min_(rt_min), rt_max_(rt_max)
    {
    }

    bool operator()(const PeptideIdentification& id) const
    {
      double rt = id.getRT();
      return (rt >= rt_min_) && (rt <= rt_max_);
    }
  };


  struct IDFilter::HasMZInRange {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double mz_min_, mz_max_;

    HasMZInRange(double mz_min, double mz_max) : mz_min_(mz_min), mz_max_(mz_max)
    {
    }

    bool operator()(const PeptideIdentification& id) const
    {
      double mz = id.getMZ();
      return (mz >= mz_min_) && (mz <= mz_max_);
    }
  };


  void IDFilter::extractPeptideSequences(const vector<PeptideIdentification>& peptides, set<String>& sequences, bool ignore_mods)
  {
    for (const PeptideIdentification& pep : peptides)
    {
      for (const PeptideHit& hit : pep.getHits())
      {
        if (ignore_mods)
        {
          sequences.insert(hit.getSequence().toUnmodifiedString());
        }
        else
        {
          sequences.insert(hit.getSequence().toString());
        }
      }
    }
  }

  map<String, vector<ProteinHit>> IDFilter::extractUnassignedProteins(ConsensusMap& cmap)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String>> run_to_accessions;

    for (const auto& f : cmap)
    {
      for (const auto& pepid : f.getPeptideIdentifications())
      {
        const String& run_id = pepid.getIdentifier();
        // extract protein accessions of each peptide hit:
        for (const PeptideHit& hit : pepid.getHits())
        {
          const set<String>& current_accessions = hit.extractProteinAccessionsSet();
          run_to_accessions[run_id].insert(current_accessions.begin(), current_accessions.end());
        }
      }
    }

    vector<ProteinIdentification>& prots = cmap.getProteinIdentifications();

    map<String, vector<ProteinHit>> result {};
    for (ProteinIdentification& prot : prots)
    {
      const String& run_id = prot.getIdentifier();
      auto target = result.emplace(run_id, vector<ProteinHit> {});
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      moveMatchingItems(prot.getHits(), std::not_fn(acc_filter), target.first->second);
    }
    return result;
  }

  void IDFilter::removeUnreferencedProteins(ConsensusMap& cmap, bool include_unassigned)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String>> run_to_accessions;

    auto add_references_to_map = [&run_to_accessions](const PeptideIdentification& pepid) {
      const String& run_id = pepid.getIdentifier();
      // extract protein accessions of each peptide hit:
      for (const PeptideHit& hit : pepid.getHits())
      {
        const set<String>& current_accessions = hit.extractProteinAccessionsSet();
        run_to_accessions[run_id].insert(current_accessions.begin(), current_accessions.end());
      }
    };
    cmap.applyFunctionOnPeptideIDs(add_references_to_map, include_unassigned);

    vector<ProteinIdentification>& prots = cmap.getProteinIdentifications();

    for (ProteinIdentification& prot : prots)
    {
      const String& run_id = prot.getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      keepMatchingItems(prot.getHits(), acc_filter);
    }
  }

  void IDFilter::removeUnreferencedProteins(ProteinIdentification& proteins, const vector<PeptideIdentification>& peptides)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String>> run_to_accessions;
    for (const PeptideIdentification& pep : peptides)
    {
      const String& run_id = pep.getIdentifier();
      // extract protein accessions of each peptide hit:
      for (const PeptideHit& hit : pep.getHits())
      {
        const set<String>& current_accessions = hit.extractProteinAccessionsSet();
        run_to_accessions[run_id].insert(current_accessions.begin(), current_accessions.end());
      }
    }

    const String& run_id = proteins.getIdentifier();
    const unordered_set<String>& accessions = run_to_accessions[run_id];
    struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
    keepMatchingItems(proteins.getHits(), acc_filter);
  }

  void IDFilter::removeUnreferencedProteins(vector<ProteinIdentification>& proteins, const vector<PeptideIdentification>& peptides)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String>> run_to_accessions;
    for (const PeptideIdentification& pep : peptides)
    {
      const String& run_id = pep.getIdentifier();
      // extract protein accessions of each peptide hit:
      for (const PeptideHit& hit : pep.getHits())
      {
        const set<String>& current_accessions = hit.extractProteinAccessionsSet();
        run_to_accessions[run_id].insert(current_accessions.begin(), current_accessions.end());
      }
    }

    for (ProteinIdentification& prot : proteins)
    {
      const String& run_id = prot.getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      keepMatchingItems(prot.getHits(), acc_filter);
    }
  }

  // TODO write version where you look up in a specific run (e.g. first inference run)
  void IDFilter::updateProteinReferences(ConsensusMap& cmap, bool remove_peptides_without_reference)
  {
    vector<ProteinIdentification>& proteins = cmap.getProteinIdentifications();
    // collect valid protein accessions for each ID run:
    map<String, unordered_set<String>> run_to_accessions;
    for (const ProteinIdentification& prot : proteins)
    {
      const String& run_id = prot.getIdentifier();
      for (const ProteinHit& hit : prot.getHits())
      {
        run_to_accessions[run_id].insert(hit.getAccession());
      }
    }

    auto check_prots_avail = [&run_to_accessions, &remove_peptides_without_reference](PeptideIdentification& pep) -> void {
      const String& run_id = pep.getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      // check protein accessions of each peptide hit
      for (PeptideHit& hit : pep.getHits())
      {
        // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
        // can't use "keepMatchingItems":
        vector<PeptideEvidence> evidences;
        remove_copy_if(hit.getPeptideEvidences().begin(), hit.getPeptideEvidences().end(), back_inserter(evidences), std::not_fn(acc_filter));
        hit.setPeptideEvidences(evidences);
      }

      if (remove_peptides_without_reference)
      {
        removeMatchingItems(pep.getHits(), HasNoEvidence());
      }
    };

    cmap.applyFunctionOnPeptideIDs(check_prots_avail);
  }

  void IDFilter::updateProteinReferences(ConsensusMap& cmap, const ProteinIdentification& ref_run, bool remove_peptides_without_reference)
  {
    // collect valid protein accessions for each ID run:
    unordered_set<String> accessions_avail;

    for (const ProteinHit& hit : ref_run.getHits())
    {
      accessions_avail.insert(hit.getAccession());
    }

    // TODO could be refactored and pulled out
    auto check_prots_avail = [&accessions_avail, &remove_peptides_without_reference](PeptideIdentification& pep) -> void {
      const unordered_set<String>& accessions = accessions_avail;
      struct HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      // check protein accessions of each peptide hit
      for (PeptideHit& hit : pep.getHits())
      {
        // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
        // can't use "keepMatchingItems":
        vector<PeptideEvidence> evidences;
        remove_copy_if(hit.getPeptideEvidences().begin(), hit.getPeptideEvidences().end(), back_inserter(evidences), std::not_fn(acc_filter));
        hit.setPeptideEvidences(evidences);
      }

      if (remove_peptides_without_reference)
      {
        removeMatchingItems(pep.getHits(), HasNoEvidence());
      }
    };

    cmap.applyFunctionOnPeptideIDs(check_prots_avail);
  }

  void IDFilter::updateProteinReferences(vector<PeptideIdentification>& peptides, const vector<ProteinIdentification>& proteins, bool remove_peptides_without_reference)
  {
    // collect valid protein accessions for each ID run:
    map<String, unordered_set<String>> run_to_accessions;
    for (const ProteinIdentification& prot : proteins)
    {
      const String& run_id = prot.getIdentifier();
      for (const ProteinHit& hit : prot.getHits())
      {
        run_to_accessions[run_id].insert(hit.getAccession());
      }
    }

    for (PeptideIdentification& pep : peptides)
    {
      const String& run_id = pep.getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      // check protein accessions of each peptide hit
      for (PeptideHit& hit : pep.getHits())
      {
        // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
        // can't use "keepMatchingItems":
        vector<PeptideEvidence> evidences;
        remove_copy_if(hit.getPeptideEvidences().begin(), hit.getPeptideEvidences().end(), back_inserter(evidences), std::not_fn(acc_filter));
        hit.setPeptideEvidences(evidences);
      }

      if (remove_peptides_without_reference)
      {
        removeMatchingItems(pep.getHits(), HasNoEvidence());
      }
    }
  }


  bool IDFilter::updateProteinGroups(vector<ProteinIdentification::ProteinGroup>& groups, const vector<ProteinHit>& hits)
  {
    if (groups.empty())
      return true; // nothing to update

    // we'll do lots of look-ups, so use a suitable data structure:
    unordered_set<String> valid_accessions;
    for (const ProteinHit& hit : hits)
    {
      valid_accessions.insert(hit.getAccession());
    }

    bool valid = true;
    vector<ProteinIdentification::ProteinGroup> filtered_groups;
    for (ProteinIdentification::ProteinGroup& group : groups)
    {
      ProteinIdentification::ProteinGroup filtered;
      for (const String& acc : group.accessions)
      {
        if (valid_accessions.find(acc) != valid_accessions.end())
        {
          filtered.accessions.push_back(acc);
        }
      }
      if (!filtered.accessions.empty())
      {
        if (filtered.accessions.size() < group.accessions.size())
        {
          valid = false; // some proteins removed from group
        }
        filtered.probability = group.probability;
        filtered_groups.push_back(filtered);
      }
    }
    groups.swap(filtered_groups);

    return valid;
  }

  void IDFilter::removeUngroupedProteins(const vector<ProteinIdentification::ProteinGroup>& groups, vector<ProteinHit>& hits)
  {
    if (hits.empty())
    {
      return; // nothing to update
    }
    // we'll do lots of look-ups, so use a suitable data structure:
    unordered_set<String> valid_accessions;
    for (const auto& grp : groups)
    {
      valid_accessions.insert(grp.accessions.begin(), grp.accessions.end());
    }

    hits.erase(std::remove_if(hits.begin(), hits.end(), std::not_fn(HasMatchingAccessionUnordered<ProteinHit>(valid_accessions))), hits.end());
  }

  void IDFilter::keepBestPeptideHits(vector<PeptideIdentification>& peptides, bool strict)
  {
    for (PeptideIdentification& pep : peptides)
    {
      vector<PeptideHit>& hits = pep.getHits();
      if (hits.size() > 1)
      {
        pep.sort();
        double top_score = hits[0].getScore();
        bool higher_better = pep.isHigherScoreBetter();
        struct HasGoodScore<PeptideHit> good_score(top_score, higher_better);
        if (strict) // only one best score allowed
        {
          if (good_score(hits[1])) // two (or more) best-scoring hits
          {
            hits.clear();
          }
          else
          {
            hits.resize(1);
          }
        }
        else
        {
          // we could use keepMatchingHits() here, but it would be less
          // efficient (since the hits are already sorted by score):
          for (vector<PeptideHit>::iterator hit_it = ++hits.begin(); hit_it != hits.end(); ++hit_it)
          {
            if (!good_score(*hit_it))
            {
              hits.erase(hit_it, hits.end());
              break;
            }
          }
        }
      }
    }
  }

  void IDFilter::filterGroupsByScore(std::vector<ProteinIdentification::ProteinGroup>& grps, double threshold_score, bool higher_better)
  {
    const auto& pred = [&threshold_score, &higher_better](ProteinIdentification::ProteinGroup& g) {
      return ((higher_better && (threshold_score >= g.probability)) || (!higher_better && (threshold_score < g.probability)));
    };

    grps.erase(std::remove_if(grps.begin(), grps.end(), pred), grps.end());
  }

  void IDFilter::filterPeptidesByLength(vector<PeptideIdentification>& peptides, Size min_length, Size max_length)
  {
    if (min_length > 0)
    {
      struct HasMinPeptideLength length_filter(min_length);
      for (PeptideIdentification& pep : peptides)
      {
        keepMatchingItems(pep.getHits(), length_filter);
      }
    }

    if (max_length == std::numeric_limits<decltype(max_length)>::max()) return; // no upper end filtering needed

    ++max_length; // the predicate tests for ">=", we need ">"
    if (max_length > min_length)
    {
      struct HasMinPeptideLength length_filter(max_length);
      for (PeptideIdentification& pep : peptides)
      {
        removeMatchingItems(pep.getHits(), length_filter);
      }
    }
  }

  void IDFilter::filterPeptidesByCharge(vector<PeptideIdentification>& peptides, Int min_charge, Int max_charge)
  {
    struct HasMinCharge charge_filter(min_charge);
    
    for (PeptideIdentification& pep : peptides)
    {
      keepMatchingItems(pep.getHits(), charge_filter);
    }

    if (max_charge == std::numeric_limits<decltype(max_charge)>::max()) return; // no upper end filtering needed

    ++max_charge; // the predicate tests for ">=", we need ">"
    if (max_charge > min_charge)
    {     
      charge_filter = HasMinCharge(max_charge);
      for (PeptideIdentification& pep : peptides)
      {
        removeMatchingItems(pep.getHits(), charge_filter);
      }
    }    
  }


  void IDFilter::filterPeptidesByRT(vector<PeptideIdentification>& peptides, double min_rt, double max_rt)
  {
    struct HasRTInRange rt_filter(min_rt, max_rt);
    keepMatchingItems(peptides, rt_filter);
  }


  void IDFilter::filterPeptidesByMZ(vector<PeptideIdentification>& peptides, double min_mz, double max_mz)
  {
    struct HasMZInRange mz_filter(min_mz, max_mz);
    keepMatchingItems(peptides, mz_filter);
  }


  void IDFilter::filterPeptidesByMZError(vector<PeptideIdentification>& peptides, double mass_error, bool unit_ppm)
  {
    for (PeptideIdentification& pep : peptides)
    {
      struct HasLowMZError error_filter(pep.getMZ(), mass_error, unit_ppm);
      keepMatchingItems(pep.getHits(), error_filter);
    }
  }


  void IDFilter::filterPeptidesByRTPredictPValue(vector<PeptideIdentification>& peptides, const String& metavalue_key, double threshold)
  {
    Size n_initial = 0, n_metavalue = 0; // keep track of numbers of hits
    struct HasMetaValue<PeptideHit> present_filter(metavalue_key, DataValue());
    double cutoff = 1 - threshold; // why? - Hendrik
    struct HasMaxMetaValue<PeptideHit> pvalue_filter(metavalue_key, cutoff);
    for (PeptideIdentification& pep : peptides)
    {
      n_initial += pep.getHits().size();
      keepMatchingItems(pep.getHits(), present_filter);
      n_metavalue += pep.getHits().size();

      keepMatchingItems(pep.getHits(), pvalue_filter);
    }

    if (n_metavalue < n_initial)
    {
      OPENMS_LOG_WARN << "Filtering peptides by RTPredict p-value removed " << (n_initial - n_metavalue) << " of " << n_initial << " hits (total) that were missing the required meta value ('"
                      << metavalue_key << "', added by RTPredict)." << endl;
    }
  }


  void IDFilter::removePeptidesWithMatchingModifications(vector<PeptideIdentification>& peptides, const set<String>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (PeptideIdentification& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), mod_filter);
    }
  }

  void IDFilter::removePeptidesWithMatchingRegEx(vector<PeptideIdentification>& peptides, const String& regex)
  {
    const std::regex re(regex);

    // true if regex matches to parts or entire unmodified sequence
    auto regex_matches = [&re](const PeptideHit& ph) -> bool { return std::regex_search(ph.getSequence().toUnmodifiedString(), re); };

    for (auto& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), regex_matches);
    }
  }

  void IDFilter::keepPeptidesWithMatchingModifications(vector<PeptideIdentification>& peptides, const set<String>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (PeptideIdentification& pep : peptides)
    {
      keepMatchingItems(pep.getHits(), mod_filter);
    }
  }


  void IDFilter::removePeptidesWithMatchingSequences(vector<PeptideIdentification>& peptides, const vector<PeptideIdentification>& bad_peptides, bool ignore_mods)
  {
    set<String> bad_seqs;
    extractPeptideSequences(bad_peptides, bad_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(bad_seqs, ignore_mods);
    for (PeptideIdentification& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), seq_filter);
    }
  }


  void IDFilter::keepPeptidesWithMatchingSequences(vector<PeptideIdentification>& peptides, const vector<PeptideIdentification>& good_peptides, bool ignore_mods)
  {
    set<String> good_seqs;
    extractPeptideSequences(good_peptides, good_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(good_seqs, ignore_mods);
    for (PeptideIdentification& pep : peptides)
    {
      keepMatchingItems(pep.getHits(), seq_filter);
    }
  }


  void IDFilter::keepUniquePeptidesPerProtein(vector<PeptideIdentification>& peptides)
  {
    Size n_initial = 0, n_metavalue = 0; // keep track of numbers of hits
    struct HasMetaValue<PeptideHit> present_filter("protein_references", DataValue());
    struct HasMetaValue<PeptideHit> unique_filter("protein_references", DataValue("unique"));
    for (PeptideIdentification& pep : peptides)
    {
      n_initial += pep.getHits().size();
      keepMatchingItems(pep.getHits(), present_filter);
      n_metavalue += pep.getHits().size();

      keepMatchingItems(pep.getHits(), unique_filter);
    }

    if (n_metavalue < n_initial)
    {
      OPENMS_LOG_WARN << "Filtering peptides by unique match to a protein removed " << (n_initial - n_metavalue) << " of " << n_initial << " hits (total) that were missing the required meta value "
                      << "('protein_references', added by PeptideIndexer)." << endl;
    }
  }


  // @TODO: generalize this to protein hits?
  void IDFilter::removeDuplicatePeptideHits(vector<PeptideIdentification>& peptides, bool seq_only)
  {
    for (PeptideIdentification& pep : peptides)
    {
      vector<PeptideHit> filtered_hits;
      if (seq_only)
      {
        set<AASequence> seqs;
        for (PeptideHit& hit : pep.getHits())
        {
          if (seqs.insert(hit.getSequence()).second) // new sequence
          {
            filtered_hits.push_back(hit);
          }
        }
      }
      else
      {
        // there's no "PeptideHit::operator<" defined, so we can't use a set nor
        // "sort" + "unique" from the standard library:
        for (PeptideHit& hit : pep.getHits())
        {
          if (find(filtered_hits.begin(), filtered_hits.end(), hit) == filtered_hits.end())
          {
            filtered_hits.push_back(hit);
          }
        }
      }
      pep.getHits().swap(filtered_hits);
    }
  }

  void IDFilter::keepNBestSpectra(std::vector<PeptideIdentification>& peptides, Size n)
  {
    String score_type;
    for (PeptideIdentification& p : peptides)
    {
      p.sort();
      if (score_type.empty())
      {
        score_type = p.getScoreType();
      }
      else
      {
        if (p.getScoreType() != score_type)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("PSM score types must be identical to allow proper filtering."));
        }
      }
    }

    // there might be fewer spectra identified than n -> adapt
    n = std::min(n, peptides.size());

    auto has_better_peptidehit = [](const PeptideIdentification& l, const PeptideIdentification& r) {
      if (r.getHits().empty())
      {
        return true; // right has no hit? -> left is better
      }
      if (l.getHits().empty())
      {
        return false; // left has no hit but right has a hit? -> right is better
      }
      const bool higher_better = l.isHigherScoreBetter();
      const double l_score = l.getHits()[0].getScore();
      const double r_score = r.getHits()[0].getScore();

      // both have hits? better score of best PSM is better
      if (higher_better)
      {
        return l_score > r_score;
      }
      return l_score < r_score;
    };

    std::partial_sort(peptides.begin(), peptides.begin() + n, peptides.end(), has_better_peptidehit);
    peptides.resize(n);
  }

  void IDFilter::keepBestMatchPerObservation(IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref)
  {
    if (id_data.getObservationMatches().size() <= 1)
      return; // nothing to do

    vector<IdentificationData::ObservationMatchRef> best_matches = id_data.getBestMatchPerObservation(score_ref);

    auto best_match_it = best_matches.begin();

    // predicate to compare the best match(es) to all (ordered) observation matches
    // returns false if the current om is a best match (-> not to be removed)
    // returns true if an inferior om was found (-> will be removed)
    auto has_worse_score = [&best_match_it](IdentificationData::ObservationMatchRef it) -> bool {
      if (it == *best_match_it)
      {
        ++best_match_it;
        return false;
      }
      return true;
    };

    id_data.removeObservationMatchesIf(has_worse_score);
  }

  void IDFilter::filterObservationMatchesByScore(IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref, double cutoff)
  {
    // predicate to compare the score of observation matches to a cutoff
    // returns true if the current score is worse than the cutoff
    // returns false otherwise
    auto is_worse_than_cutoff = [&](IdentificationData::ObservationMatchRef it) -> bool {
      pair<double, bool> score = it->getScore(score_ref);
      return !score.second || score_ref->isBetterScore(cutoff, score.first);
    };

    id_data.removeObservationMatchesIf(is_worse_than_cutoff);
  }

  void IDFilter::removeDecoys(IdentificationData& id_data)
  {
    // predicate to compare the target/decoy status of a parent sequence
    // returns true if decoy
    // returns false if target
    auto is_decoy = [&](IdentificationData::ParentSequenceRef it) -> bool { return it->is_decoy; };

    id_data.removeParentSequencesIf(is_decoy);
  }

} // namespace OpenMS
