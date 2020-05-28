// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS
{
  IDFilter::IDFilter()
  {
  }

  IDFilter::~IDFilter()
  {
  }


  struct IDFilter::HasMinPeptideLength
  {
    typedef PeptideHit argument_type; // for use as a predicate

    Size length_;

    explicit HasMinPeptideLength(Size length):
      length_(length)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getSequence().size() >= length_;
    }
  };


  struct IDFilter::HasMinCharge
  {
    typedef PeptideHit argument_type; // for use as a predicate

    Int charge_;

    explicit HasMinCharge(Int charge):
      charge_(charge)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getCharge() >= charge_;
    }
  };


  struct IDFilter::HasLowMZError
  {
    typedef PeptideHit argument_type; // for use as a predicate

    double precursor_mz_, tolerance_;

    HasLowMZError(double precursor_mz, double tolerance, bool unit_ppm):
      precursor_mz_(precursor_mz), tolerance_(tolerance)
    {
      if (unit_ppm) this->tolerance_ *= precursor_mz / 1.0e6;
    }

    bool operator()(const PeptideHit& hit) const
    {
      Int z = hit.getCharge();
      if (z == 0) z = 1;
      double peptide_mz = (hit.getSequence().getMonoWeight(Residue::Full, z) /
                           double(z));
      return fabs(precursor_mz_ - peptide_mz) <= tolerance_;
    }
  };


  struct IDFilter::HasMatchingModification
  {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& mods_;

    explicit HasMatchingModification(const set<String>& mods):
      mods_(mods)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      const AASequence& seq = hit.getSequence();
      if (mods_.empty()) return seq.isModified();

      for (Size i = 0; i < seq.size(); ++i)
      {
        if (seq[i].isModified())
        {
          String mod_name = seq[i].getModification()->getFullId();
          if (mods_.count(mod_name) > 0) return true;
        }
      }

      // terminal modifications:
      if (seq.hasNTerminalModification())
      {
        String mod_name = seq.getNTerminalModification()->getFullId();
        if (mods_.count(mod_name) > 0) return true;
      }
      if (seq.hasCTerminalModification())
      {
        String mod_name = seq.getCTerminalModification()->getFullId();
        if (mods_.count(mod_name) > 0) return true;
      }

      return false;
    }
  };


  struct IDFilter::HasMatchingSequence
  {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& sequences_;
    bool ignore_mods_;

    explicit HasMatchingSequence(const set<String>& sequences, bool ignore_mods = false):
      sequences_(sequences), ignore_mods_(ignore_mods)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      const String& query = (ignore_mods_ ?
                             hit.getSequence().toUnmodifiedString() :
                             hit.getSequence().toString());
      return (sequences_.count(query) > 0);
    }
  };


  struct IDFilter::HasNoEvidence
  {
    typedef PeptideHit argument_type; // for use as a predicate

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getPeptideEvidences().empty();
    }
  };

  struct IDFilter::HasRTInRange
  {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double rt_min_, rt_max_;

    HasRTInRange(double rt_min, double rt_max):
      rt_min_(rt_min), rt_max_(rt_max)
    {}

    bool operator()(const PeptideIdentification& id) const
    {
      double rt = id.getRT();
      return (rt >= rt_min_) && (rt <= rt_max_);
    }
  };


  struct IDFilter::HasMZInRange
  {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double mz_min_, mz_max_;

    HasMZInRange(double mz_min, double mz_max):
      mz_min_(mz_min), mz_max_(mz_max)
    {}

    bool operator()(const PeptideIdentification& id) const
    {
      double mz = id.getMZ();
      return (mz >= mz_min_) && (mz <= mz_max_);
    }
  };


  void IDFilter::extractPeptideSequences(
    const vector<PeptideIdentification>& peptides, set<String>& sequences,
    bool ignore_mods)
  {
    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      for (vector<PeptideHit>::const_iterator hit_it =
             pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
           ++hit_it)
      {
        if (ignore_mods)
        {
          sequences.insert(hit_it->getSequence().toUnmodifiedString());
        }
        else
        {
          sequences.insert(hit_it->getSequence().toString());
        }
      }
    }
  }

  map<String,vector<ProteinHit>> IDFilter::extractUnassignedProteins(ConsensusMap& cmap)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String> > run_to_accessions;

    for (const auto& f : cmap)
    {
      for (const auto& pepid : f.getPeptideIdentifications())
      {
        const String& run_id = pepid.getIdentifier();
        // extract protein accessions of each peptide hit:
        for (vector<PeptideHit>::const_iterator hit_it =
            pepid.getHits().begin(); hit_it != pepid.getHits().end();
             ++hit_it)
        {

          const set<String>& current_accessions =
              hit_it->extractProteinAccessionsSet();

          run_to_accessions[run_id].insert(current_accessions.begin(),
                                           current_accessions.end());
        }
      }
    }

    vector<ProteinIdentification>& prots = cmap.getProteinIdentifications();

    map<String,vector<ProteinHit>> result{};
    for (vector<ProteinIdentification>::iterator prot_it = prots.begin();
         prot_it != prots.end(); ++prot_it)
    {
      const String& run_id = prot_it->getIdentifier();
      auto target = result.emplace(run_id, vector<ProteinHit>{});
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      moveMatchingItems(prot_it->getHits(), std::not1(acc_filter), target.first->second);
    }
    return result;
  }

  void IDFilter::removeUnreferencedProteins(ConsensusMap& cmap, bool include_unassigned)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String> > run_to_accessions;

    auto add_references_to_map =
        [&run_to_accessions](const PeptideIdentification& pepid)
    {
      const String& run_id = pepid.getIdentifier();
      // extract protein accessions of each peptide hit:
      for (vector<PeptideHit>::const_iterator hit_it =
          pepid.getHits().begin(); hit_it != pepid.getHits().end();
           ++hit_it)
      {

        const set<String>& current_accessions =
            hit_it->extractProteinAccessionsSet();

        run_to_accessions[run_id].insert(current_accessions.begin(),
                                         current_accessions.end());
      }
    };
    cmap.applyFunctionOnPeptideIDs(add_references_to_map,include_unassigned);

    vector<ProteinIdentification>& prots = cmap.getProteinIdentifications();

    for (vector<ProteinIdentification>::iterator prot_it = prots.begin();
         prot_it != prots.end(); ++prot_it)
    {
      const String& run_id = prot_it->getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      keepMatchingItems(prot_it->getHits(), acc_filter);
    }
  }


  void IDFilter::removeUnreferencedProteins(
    vector<ProteinIdentification>& proteins,
    const vector<PeptideIdentification>& peptides)
  {
    // collect accessions that are referenced by peptides for each ID run:
    map<String, unordered_set<String> > run_to_accessions;
    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      const String& run_id = pep_it->getIdentifier();
      // extract protein accessions of each peptide hit:
      for (vector<PeptideHit>::const_iterator hit_it =
             pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
           ++hit_it)
      {
        const set<String>& current_accessions = 
          hit_it->extractProteinAccessionsSet();

        run_to_accessions[run_id].insert(current_accessions.begin(),
                                         current_accessions.end());
      }
    }

    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      const String& run_id = prot_it->getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      keepMatchingItems(prot_it->getHits(), acc_filter);
    }
  }

  void IDFilter::updateProteinReferences(
      ConsensusMap& cmap,
      bool remove_peptides_without_reference)
  {
    vector<ProteinIdentification>& proteins = cmap.getProteinIdentifications();
    // collect valid protein accessions for each ID run:
    map<String, unordered_set<String> > run_to_accessions;
    for (vector<ProteinIdentification>::const_iterator prot_it =
        proteins.begin(); prot_it != proteins.end(); ++prot_it)
    {
      const String& run_id = prot_it->getIdentifier();
      for (vector<ProteinHit>::const_iterator hit_it =
          prot_it->getHits().begin(); hit_it != prot_it->getHits().end();
           ++hit_it)
      {
        run_to_accessions[run_id].insert(hit_it->getAccession());
      }
    }

    auto check_prots_avail = [&run_to_accessions,&remove_peptides_without_reference]
        (PeptideIdentification& pep_it) -> void
    {
      const String& run_id = pep_it.getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      // check protein accessions of each peptide hit
      for (vector<PeptideHit>::iterator hit_it = pep_it.getHits().begin();
           hit_it != pep_it.getHits().end(); ++hit_it)
      {
        // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
        // can't use "keepMatchingItems":
        vector<PeptideEvidence> evidences;
        remove_copy_if(hit_it->getPeptideEvidences().begin(),
                       hit_it->getPeptideEvidences().end(),
                       back_inserter(evidences),
                       not1(acc_filter));
        hit_it->setPeptideEvidences(evidences);
      }

      if (remove_peptides_without_reference)
      {
        removeMatchingItems(pep_it.getHits(), HasNoEvidence());
      }
    };

    cmap.applyFunctionOnPeptideIDs(check_prots_avail);
  }

  void IDFilter::updateProteinReferences(
    vector<PeptideIdentification>& peptides,
    const vector<ProteinIdentification>& proteins,
    bool remove_peptides_without_reference)
  {
    // collect valid protein accessions for each ID run:
    map<String, unordered_set<String> > run_to_accessions;
    for (vector<ProteinIdentification>::const_iterator prot_it =
           proteins.begin(); prot_it != proteins.end(); ++prot_it)
    {
      const String& run_id = prot_it->getIdentifier();
      for (vector<ProteinHit>::const_iterator hit_it =
             prot_it->getHits().begin(); hit_it != prot_it->getHits().end();
           ++hit_it)
      {
        run_to_accessions[run_id].insert(hit_it->getAccession());
      }
    }

    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      const String& run_id = pep_it->getIdentifier();
      const unordered_set<String>& accessions = run_to_accessions[run_id];
      struct HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      // check protein accessions of each peptide hit
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
        // can't use "keepMatchingItems":
        vector<PeptideEvidence> evidences;
        remove_copy_if(hit_it->getPeptideEvidences().begin(),
                       hit_it->getPeptideEvidences().end(),
                       back_inserter(evidences),
                       not1(acc_filter));
        hit_it->setPeptideEvidences(evidences);
      }

      if (remove_peptides_without_reference)
      {
        removeMatchingItems(pep_it->getHits(), HasNoEvidence());
      }
    }
  }


  bool IDFilter::updateProteinGroups(
    vector<ProteinIdentification::ProteinGroup>& groups,
    const vector<ProteinHit>& hits)
  {
    if (groups.empty()) return true; // nothing to update

    // we'll do lots of look-ups, so use a suitable data structure:
    unordered_set<String> valid_accessions;
    for (vector<ProteinHit>::const_iterator hit_it = hits.begin();
         hit_it != hits.end(); ++hit_it)
    {
      valid_accessions.insert(hit_it->getAccession());
    }

    bool valid = true;
    vector<ProteinIdentification::ProteinGroup> filtered_groups;
    for (vector<ProteinIdentification::ProteinGroup>::iterator group_it =
           groups.begin(); group_it != groups.end(); ++group_it)
    {
      ProteinIdentification::ProteinGroup filtered;
      for (const String& acc : group_it->accessions)
      {
	if (valid_accessions.find(acc) != valid_accessions.end())
		filtered.accessions.push_back(acc);
      }
      if (!filtered.accessions.empty())
      {
        if (filtered.accessions.size() < group_it->accessions.size())
        {
          valid = false; // some proteins removed from group
        }
        filtered.probability = group_it->probability;
        filtered_groups.push_back(filtered);
      }
    }
    groups.swap(filtered_groups);

    return valid;
  }


  void IDFilter::keepBestPeptideHits(vector<PeptideIdentification>& peptides,
                                     bool strict)
  {
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      vector<PeptideHit>& hits = pep_it->getHits();
      if (hits.size() > 1)
      {
        pep_it->sort();
        double top_score = hits[0].getScore();
        bool higher_better = pep_it->isHigherScoreBetter();
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
          for (vector<PeptideHit>::iterator hit_it = ++hits.begin();
               hit_it != hits.end(); ++hit_it)
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


  void IDFilter::filterPeptidesByLength(vector<PeptideIdentification>& peptides,
                                        Size min_length, Size max_length)
  {
    if (min_length > 0)
    {
      struct HasMinPeptideLength length_filter(min_length);
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        keepMatchingItems(pep_it->getHits(), length_filter);
      }
    }
    ++max_length; // the predicate tests for ">=", we need ">"
    if (max_length > min_length)
    {
      struct HasMinPeptideLength length_filter(max_length);
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        removeMatchingItems(pep_it->getHits(), length_filter);
      }
    }
  }


  void IDFilter::filterPeptidesByCharge(vector<PeptideIdentification>& peptides,
                                        Int min_charge, Int max_charge)
  {
    struct HasMinCharge charge_filter(min_charge);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      keepMatchingItems(pep_it->getHits(), charge_filter);
    }
    ++max_charge; // the predicate tests for ">=", we need ">"
    if (max_charge > min_charge)
    {
      charge_filter = HasMinCharge(max_charge);
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        removeMatchingItems(pep_it->getHits(), charge_filter);
      }
    }
  }


  void IDFilter::filterPeptidesByRT(vector<PeptideIdentification>& peptides,
                                    double min_rt, double max_rt)
  {
    struct HasRTInRange rt_filter(min_rt, max_rt);
    keepMatchingItems(peptides, rt_filter);
  }


  void IDFilter::filterPeptidesByMZ(vector<PeptideIdentification>& peptides,
                                    double min_mz, double max_mz)
  {
    struct HasMZInRange mz_filter(min_mz, max_mz);
    keepMatchingItems(peptides, mz_filter);
  }


  void IDFilter::filterPeptidesByMZError(
    vector<PeptideIdentification>& peptides, double mass_error, bool unit_ppm)
  {
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      struct HasLowMZError error_filter(pep_it->getMZ(), mass_error, unit_ppm);
      keepMatchingItems(pep_it->getHits(), error_filter);
    }
  }


  void IDFilter::filterPeptidesByRTPredictPValue(
    vector<PeptideIdentification>& peptides, const String& metavalue_key,
    double threshold)
  {
    Size n_initial = 0, n_metavalue = 0; // keep track of numbers of hits
    struct HasMetaValue<PeptideHit> present_filter(metavalue_key, DataValue());
    double cutoff = 1 - threshold; // why? - Hendrik
    struct HasMaxMetaValue<PeptideHit> pvalue_filter(metavalue_key, cutoff);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      n_initial += pep_it->getHits().size();
      keepMatchingItems(pep_it->getHits(), present_filter);
      n_metavalue += pep_it->getHits().size();

      keepMatchingItems(pep_it->getHits(), pvalue_filter);
    }

    if (n_metavalue < n_initial)
    {
      OPENMS_LOG_WARN << "Filtering peptides by RTPredict p-value removed "
               << (n_initial - n_metavalue) << " of " << n_initial
               << " hits (total) that were missing the required meta value ('"
               << metavalue_key << "', added by RTPredict)." << endl;
    }
  }


  void IDFilter::removePeptidesWithMatchingModifications(
    vector<PeptideIdentification>& peptides,
    const set<String>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      removeMatchingItems(pep_it->getHits(), mod_filter);
    }
  }


  void IDFilter::keepPeptidesWithMatchingModifications(
    vector<PeptideIdentification>& peptides,
    const set<String>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      keepMatchingItems(pep_it->getHits(), mod_filter);
    }
  }


  void IDFilter::removePeptidesWithMatchingSequences(
    vector<PeptideIdentification>& peptides,
    const vector<PeptideIdentification>& bad_peptides, bool ignore_mods)
  {
    set<String> bad_seqs;
    extractPeptideSequences(bad_peptides, bad_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(bad_seqs, ignore_mods);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      removeMatchingItems(pep_it->getHits(), seq_filter);
    }
  }


  void IDFilter::keepPeptidesWithMatchingSequences(
    vector<PeptideIdentification>& peptides,
    const vector<PeptideIdentification>& good_peptides, bool ignore_mods)
  {
    set<String> good_seqs;
    extractPeptideSequences(good_peptides, good_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(good_seqs, ignore_mods);
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      keepMatchingItems(pep_it->getHits(), seq_filter);
    }
  }


  void IDFilter::keepUniquePeptidesPerProtein(vector<PeptideIdentification>&
                                              peptides)
  {
    Size n_initial = 0, n_metavalue = 0; // keep track of numbers of hits
    struct HasMetaValue<PeptideHit> present_filter("protein_references",
                                                   DataValue());
    struct HasMetaValue<PeptideHit> unique_filter("protein_references",
                                                  DataValue("unique"));
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      n_initial += pep_it->getHits().size();
      keepMatchingItems(pep_it->getHits(), present_filter);
      n_metavalue += pep_it->getHits().size();

      keepMatchingItems(pep_it->getHits(), unique_filter);
    }

    if (n_metavalue < n_initial)
    {
      OPENMS_LOG_WARN << "Filtering peptides by unique match to a protein removed "
               << (n_initial - n_metavalue) << " of " << n_initial
               << " hits (total) that were missing the required meta value "
               << "('protein_references', added by PeptideIndexer)." << endl;
    }
  }


  // @TODO: generalize this to protein hits?
  void IDFilter::removeDuplicatePeptideHits(vector<PeptideIdentification>&
                                            peptides, bool seq_only)
  {
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      vector<PeptideHit> filtered_hits;
      if (seq_only)
      {
        set<AASequence> seqs;
        for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
             hit_it != pep_it->getHits().end(); ++hit_it)
        {
          if (seqs.insert(hit_it->getSequence()).second) // new sequence
          {
            filtered_hits.push_back(*hit_it);
          }
        }
      }
      else
      {
        // there's no "PeptideHit::operator<" defined, so we can't use a set nor
        // "sort" + "unique" from the standard library:
        for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
             hit_it != pep_it->getHits().end(); ++hit_it)
        {
          if (find(filtered_hits.begin(), filtered_hits.end(), *hit_it) ==
              filtered_hits.end())
          {
            filtered_hits.push_back(*hit_it);
          }
        }
      }
      pep_it->getHits().swap(filtered_hits);
    }
  }


  void IDFilter::keepBestMatchPerQuery(
    IdentificationData& id_data,
    IdentificationData::ScoreTypeRef score_ref)
  {
    if (id_data.getMoleculeQueryMatches().size() <= 1) return; // nothing to do

    vector<IdentificationData::QueryMatchRef> best_matches =
      id_data.getBestMatchPerQuery(score_ref);
    auto best_match_it = best_matches.begin();
    for (auto it = id_data.query_matches_.begin();
         it != id_data.query_matches_.end(); )
    {
      if (it == *best_match_it)
      {
        ++it;
        ++best_match_it;
      }
      else
      {
        it = id_data.query_matches_.erase(it);
      }
    }

    id_data.cleanup();
  }


  void IDFilter::filterQueryMatchesByScore(
    IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref,
    double cutoff)
  {
    id_data.removeFromSetIf_(
      id_data.query_matches_, [&](IdentificationData::QueryMatchRef it) -> bool
      {
        pair<double, bool> score = it->getScore(score_ref);
        return !score.second || score_ref->isBetterScore(cutoff, score.first);
      });

    id_data.cleanup();
  }


  void IDFilter::removeDecoys(IdentificationData& id_data)
  {
    Size n_parents = id_data.getParentMolecules().size();
    id_data.removeFromSetIf_(
      id_data.parent_molecules_,
      [&](IdentificationData::ParentMoleculeRef it) -> bool
      {
        return it->is_decoy;
      });

    if (id_data.getParentMolecules().size() < n_parents) id_data.cleanup();
  }

} // namespace OpenMS
