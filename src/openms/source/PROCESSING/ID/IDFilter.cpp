// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <cmath>
#include <regex>

using namespace std;

namespace OpenMS
{
  namespace // anonymous namespace for internal helpers
  {
    /// Collect protein accessions referenced by peptide hits, grouped by run identifier
    std::map<std::string, std::unordered_set<std::string>> collectReferencedAccessions(const PeptideIdentificationList& peptides)
    {
      std::map<std::string, std::unordered_set<std::string>> run_to_accessions;
      for (const PeptideIdentification& pep : peptides)
      {
        const std::string& run_id = pep.getIdentifier();
        for (const PeptideHit& hit : pep.getHits())
        {
          const set<std::string>& current_accessions = hit.extractProteinAccessionsSet();
          run_to_accessions[run_id].insert(current_accessions.begin(), current_accessions.end());
        }
      }
      return run_to_accessions;
    }

    /// Collect valid protein accessions from protein identifications, grouped by run identifier
    std::map<std::string, std::unordered_set<std::string>> collectProteinAccessions(const std::vector<ProteinIdentification>& proteins)
    {
      std::map<std::string, std::unordered_set<std::string>> run_to_accessions;
      for (const ProteinIdentification& prot : proteins)
      {
        const std::string& run_id = prot.getIdentifier();
        for (const ProteinHit& hit : prot.getHits())
        {
          run_to_accessions[run_id].insert(hit.getAccession());
        }
      }
      return run_to_accessions;
    }

    /// Filter peptide evidences to keep only those referencing proteins in the given accession set
    void filterEvidencesByAccessions(PeptideHit& hit, const std::unordered_set<std::string>& accessions)
    {
      IDFilter::HasMatchingAccessionUnordered<PeptideEvidence> acc_filter(accessions);
      vector<PeptideEvidence> evidences;
      remove_copy_if(hit.getPeptideEvidences().begin(), hit.getPeptideEvidences().end(),
                     back_inserter(evidences), std::not_fn(acc_filter));
      hit.setPeptideEvidences(evidences);
    }

    /// Process a PeptideIdentification to filter evidences and optionally remove hits without references
    void filterPeptideReferences(PeptideIdentification& pep,
                                 const std::unordered_set<std::string>& accessions,
                                 bool remove_peptides_without_reference)
    {
      for (PeptideHit& hit : pep.getHits())
      {
        filterEvidencesByAccessions(hit, accessions);
      }
      if (remove_peptides_without_reference)
      {
        auto has_no_evidence = [](const PeptideHit& hit) { return hit.getPeptideEvidences().empty(); };
        IDFilter::removeMatchingItems(pep.getHits(), has_no_evidence);
      }
    }
  } // anonymous namespace

  namespace // identifications of maps as identification data
  {
    using ID = IdentificationData;

    /**
      The inference result that holds the proteins of @p run, as its legacy protein run: the one that covers the run.
      Without one, the proteins of the run are its database sequences.
    */
    const ID::InferenceResult* proteinResult(const ID& data, const ID::Run& run)
    {
      const ID::InferenceResult* found = nullptr;
      for (const auto& result : data.getInferenceResults())
      {
        for (const auto& input : result.inputs)
        {
          if (input.run_uuid != run.getUuid() || found == &result) continue;
          if (found != nullptr)
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Several inference results cover run '" + run.getIdentifier() + "'; its proteins are ambiguous");
          }
          found = &result;
        }
      }
      return found;
    }

    /// The protein run of a run (its legacy protein run): its inference result, else the run itself
    std::string proteinRunKey(const ID& data, const ID::Run& run)
    {
      const auto* result = proteinResult(data, run);
      return result != nullptr ? "result:" + result->identifier : "run:" + run.getUuid();
    }

    /// The accessions of the proteins of @p run (see proteinResult())
    unordered_set<std::string> proteinAccessions(const ID& data, const ID::Run& run)
    {
      unordered_set<std::string> accessions;
      if (const auto* result = proteinResult(data, run))
      {
        for (const auto& hit : result->proteins.getHits()) accessions.insert(hit.getAccession());
      }
      else if (run.getDatabaseSequences())
      {
        for (const auto& sequence : *run.getDatabaseSequences()) accessions.insert(sequence.accession);
      }
      return accessions;
    }

    /// The matches that the (top-level) features of @p map link
    template<class MapType>
    set<ID::MatchReference> linkedMatches(const MapType& map)
    {
      set<ID::MatchReference> linked;
      for (const auto& feature : map) linked.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
      return linked;
    }

    /**
      Keep the proteins for which @p keep (protein run key, accession) returns true: the hits of inference results,
      and the database sequences of runs without one. The hits that @p removed gets are those of an inference result
      as its legacy protein run has them, and the database sequences as protein hits.
    */
    void keepProteins(ID& data, const std::function<bool(const std::string&, const std::string&)>& keep,
                      map<std::string, vector<ProteinHit>>* removed = nullptr)
    {
      auto results = data.getInferenceResults();
      bool changed = false;
      for (auto& result : results)
      {
        const std::string key = "result:" + result.identifier;
        const auto complete = removed != nullptr ? IdentificationDataAdapter::proteinHits(data, result) : vector<ProteinHit>();
        vector<ProteinHit>* extracted = nullptr;
        if (removed != nullptr)
        {
          extracted = &(*removed)[result.proteins.getIdentifier().empty() ? result.identifier : result.proteins.getIdentifier()];
        }
        vector<ProteinHit> kept;
        const auto& hits = result.proteins.getHits();
        for (Size i = 0; i < hits.size(); ++i)
        {
          if (keep(key, hits[i].getAccession()))
          {
            kept.push_back(hits[i]);
          }
          else if (extracted != nullptr)
          {
            extracted->push_back(complete[i]);
          }
        }
        if (kept.size() == hits.size()) continue;
        result.proteins.setHits(kept);
        // aliases of proteins that are neither hits nor group members are not used
        std::set<std::string> used;
        for (const auto& hit : kept) used.insert(hit.getAccession());
        for (const auto* groups : {&result.proteins.getProteinGroups(), &result.proteins.getIndistinguishableProteins()})
        {
          for (const auto& group : *groups) used.insert(group.accessions.begin(), group.accessions.end());
        }
        std::erase_if(result.qualified_accessions, [&](const auto& alias) { return !used.contains(alias.first); });
        changed = true;
      }
      for (const auto& current : data.getRuns())
      {
        if (proteinResult(data, current) != nullptr || !current.getDatabaseSequences()) continue;
        const std::string key = "run:" + current.getUuid();
        const auto& sequences = *current.getDatabaseSequences();
        vector<ID::DatabaseSequence> kept;
        for (const auto& sequence : sequences)
        {
          if (keep(key, sequence.accession)) kept.push_back(sequence);
        }
        if (kept.size() == sequences.size()) continue;
        if (removed != nullptr)
        {
          auto& extracted = (*removed)[IdentificationDataAdapter::legacyIdentifier(current)];
          for (const auto& hit : IdentificationDataAdapter::proteinHits(current))
          {
            if (!keep(key, hit.getAccession())) extracted.push_back(hit);
          }
        }
        data.getRun(current.getIdentifier()).setDatabaseSequences(std::move(kept));
      }
      if (changed)
      {
        data.clearInferenceResults();
        for (auto& result : results) data.addInferenceResult(std::move(result));
      }
    }

    /// The accessions that matches of each protein run refer to: all matches, or those that features link
    template<class MapType>
    map<std::string, unordered_set<std::string>> referencedAccessions(const MapType& map, bool include_unassigned)
    {
      const auto& data = map.getIdentificationData();
      const auto linked = include_unassigned ? set<ID::MatchReference>() : linkedMatches(map);
      std::map<std::string, unordered_set<std::string>> referenced;
      for (const auto& run : data.getRuns())
      {
        auto& accessions = referenced[proteinRunKey(data, run)];
        for (const auto& source : run.getSources())
        {
          for (const auto& query : source.identifications)
          {
            for (const auto& match : query.getMatches())
            {
              if (!include_unassigned && !linked.contains({run.getUuid(), match.getId()})) continue;
              for (const auto& evidence : match.sequence_evidence) accessions.insert(evidence.accession);
            }
          }
        }
      }
      return referenced;
    }

    /// Remove the sequence evidence of matches whose accession @p allowed (of the run) does not contain, and optionally
    /// the matches left without any
    template<class MapType>
    void filterReferences(MapType& map, const std::function<unordered_set<std::string>(const ID::Run&)>& allowed,
                          bool remove_peptides_without_reference)
    {
      auto& data = map.getIdentificationData();
      for (const auto& current : data.getRuns())
      {
        const auto accessions = allowed(current);
        data.getRun(current.getIdentifier()).transformMatches([&](ID::MatchData& match) {
          std::erase_if(match.sequence_evidence, [&](const ID::SequenceEvidence& evidence) { return !accessions.contains(evidence.accession); });
        });
      }
      if (remove_peptides_without_reference)
      {
        map.eraseMatches([](const ID::Run&, const ID::Identification&, const ID::Match& match) { return match.sequence_evidence.empty(); });
      }
    }

    /// The primary score of a match, or NaN; and whether higher is better
    double primaryScore(const ID::Run& run, const ID::Match& match)
    {
      const auto primary = run.getPrimaryScore();
      if (!primary) return std::numeric_limits<double>::quiet_NaN();
      return run.getScore(match.getId(), *primary).value_or(std::numeric_limits<double>::quiet_NaN());
    }

    bool higherBetter(const ID::Run& run)
    {
      const auto primary = run.getPrimaryScore();
      return !primary || run.getScoreDefinition(*primary).higher_better;
    }

    /// The matches of @p matches ordered best first by primary score, like PeptideIdentification::sort()
    vector<const ID::Match*> sortedMatches(const ID::Run& run, vector<const ID::Match*> matches)
    {
      const bool higher_better = higherBetter(run);
      std::stable_sort(matches.begin(), matches.end(), [&](const ID::Match* a, const ID::Match* b) {
        const double score_a = primaryScore(run, *a), score_b = primaryScore(run, *b);
        return higher_better ? score_a > score_b : score_a < score_b;
      });
      return matches;
    }

    template<class MapType>
    void keepNBest(MapType& map, Size n)
    {
      const auto& data = map.getIdentificationData();
      set<ID::MatchReference> keep;
      for (const auto& run : data.getRuns())
      {
        for (const auto& source : run.getSources())
        {
          for (const auto& query : source.identifications)
          {
            vector<const ID::Match*> matches;
            for (const auto& match : query.getMatches()) matches.push_back(&match);
            const auto sorted = sortedMatches(run, std::move(matches));
            for (Size i = 0; i < std::min(n, sorted.size()); ++i) keep.insert({run.getUuid(), sorted[i]->getId()});
          }
        }
      }
      map.eraseMatches([&](const ID::Run& run, const ID::Identification&, const ID::Match& match) {
        return !keep.contains({run.getUuid(), match.getId()});
      });
    }

    template<class MapType>
    void removeEmpty(MapType& map)
    {
      const auto& data = map.getIdentificationData();
      // a feature that links an identification, but none of its matches, has a legacy peptide identification without hits:
      for (auto& feature : map)
      {
        std::erase_if(feature.getIDQueries(), [&](const ID::QueryReference& reference) {
          const auto* run = data.findRunByUuid(reference.run_uuid);
          const auto* query = run != nullptr ? run->findIdentification(reference.query) : nullptr;
          if (query == nullptr) return false;
          return std::none_of(query->getMatches().begin(), query->getMatches().end(), [&](const ID::Match& match) {
            return feature.getIDMatches().contains({reference.run_uuid, match.getId()});
          });
        });
      }
      map.eraseIdentifications([](const ID::Run&, const ID::Identification& query) { return query.getMatches().empty(); });
    }

    /// The legacy peptide identifications of @p map: the matches that each feature links, then the unassigned ones
    template<class MapType>
    vector<ID::QueryMatches> legacyEntries(const MapType& map)
    {
      const auto& data = map.getIdentificationData();
      vector<ID::QueryMatches> entries;
      for (const auto& feature : map)
      {
        auto linked = feature.getLinkedIdentifications(data);
        entries.insert(entries.end(), std::make_move_iterator(linked.begin()), std::make_move_iterator(linked.end()));
      }
      auto unassigned = map.getUnassignedIdentifications();
      entries.insert(entries.end(), std::make_move_iterator(unassigned.begin()), std::make_move_iterator(unassigned.end()));
      return entries;
    }

    /// See IDFilter::annotateBestPerPeptideWithData(): the best matches per peptide sequence (and charge) of each
    /// protein run get "best_per_peptide" 1, the others of the top @p nr_best_spectrum of their identification 0
    template<class MapType>
    void annotateBestNative(MapType& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      auto& data = map.getIdentificationData();
      std::map<ID::MatchReference, int> annotations;
      // per protein run: sequence -> charge -> best match (with its score)
      std::map<std::string, std::map<std::string, std::map<int, std::pair<ID::MatchReference, double>>>> best;
      for (const auto& entry : legacyEntries(map))
      {
        const auto& run = *entry.run;
        auto& best_pep = best[proteinRunKey(data, run)];
        const bool higher_better = higherBetter(run);
        const auto sorted = sortedMatches(run, entry.matches);
        const Size n = nr_best_spectrum == 0 ? sorted.size() : std::min(nr_best_spectrum, sorted.size());
        for (Size i = 0; i < n; ++i)
        {
          const auto& match = *sorted[i];
          const ID::MatchReference reference {run.getUuid(), match.getId()};
          const double score = primaryScore(run, match);
          const std::string sequence = ignore_mods ? AASequence::fromString(match.representation).toUnmodifiedString() : match.representation;
          const int charge = ignore_charges ? 0 : match.charge;
          auto [found, inserted] = best_pep[sequence].emplace(charge, std::make_pair(reference, score));
          if (inserted)
          {
            annotations[reference] = 1;
          }
          else if ((higher_better && score > found->second.second) || (!higher_better && score < found->second.second))
          {
            annotations[found->second.first] = 0;
            annotations[reference] = 1;
            found->second = {reference, score};
          }
          else
          {
            annotations[reference] = 0;
          }
        }
      }
      for (const auto& current : data.getRuns())
      {
        vector<std::pair<ID::MatchId, ID::MatchData>> edits;
        for (const auto& source : current.getSources())
        {
          for (const auto& query : source.identifications)
          {
            for (const auto& match : query.getMatches())
            {
              const auto annotation = annotations.find({current.getUuid(), match.getId()});
              if (annotation == annotations.end()) continue;
              edits.emplace_back(match.getId(), match.getData());
              edits.back().second.setMetaValue("best_per_peptide", annotation->second);
            }
          }
        }
        auto& run = data.getRun(current.getIdentifier());
        for (const auto& [id, edited] : edits) run.replaceMatch(id, edited);
      }
    }

    template<class MapType>
    void keepBestNative(MapType& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      annotateBestNative(map, ignore_mods, ignore_charges, nr_best_spectrum);
      map.eraseMatches([](const ID::Run&, const ID::Identification&, const ID::Match& match) {
        return !match.metaValueExists("best_per_peptide") || int(match.getMetaValue("best_per_peptide")) != 1;
      });
    }

    template<class MapType>
    void filterByScore(MapType& map, double threshold_score)
    {
      // as HasGoodScore: a match without a primary score value is not at least as good
      map.eraseMatches([&](const ID::Run& run, const ID::Identification&, const ID::Match& match) {
        const auto primary = run.getPrimaryScore();
        if (!primary) return false;
        const auto score = run.getScore(match.getId(), *primary);
        if (!score) return true;
        return run.getScoreDefinition(*primary).higher_better ? !(*score >= threshold_score) : !(*score <= threshold_score);
      });
    }

    template<class MapType>
    void keepUnique(MapType& map)
    {
      Size n_initial = 0, n_missing = 0; // keep track of numbers of matches
      map.eraseMatches([&](const ID::Run&, const ID::Identification&, const ID::Match& match) {
        ++n_initial;
        if (!match.metaValueExists("protein_references"))
        {
          ++n_missing;
          return true;
        }
        return match.getMetaValue("protein_references") != DataValue("unique");
      });
      if (n_missing > 0)
      {
        OPENMS_LOG_WARN << "Filtering peptides by unique match to a protein removed " << n_missing << " of " << n_initial
                        << " hits (total) that were missing the required meta value ('protein_references', added by PeptideIndexer)." << endl;
      }
    }

    template<class MapType>
    void edit(MapType& map, const std::function<void(MapType&)>& operation)
    {
      IdentificationDataConverter::editAsIdentificationData(map, operation);
    }
  } // namespace


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

    const set<std::string>& mods_;

    explicit HasMatchingModification(const set<std::string>& mods) : mods_(mods)
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
          std::string mod_name = seq[i].getModification()->getFullId();
          if (mods_.contains(mod_name))
            return true;
        }
      }

      // terminal modifications:
      if (seq.hasNTerminalModification())
      {
        std::string mod_name = seq.getNTerminalModification()->getFullId();
        if (mods_.contains(mod_name))
          return true;
      }
      if (seq.hasCTerminalModification())
      {
        std::string mod_name = seq.getCTerminalModification()->getFullId();
        if (mods_.contains(mod_name))
          return true;
      }

      return false;
    }
  };


  struct IDFilter::HasMatchingSequence {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<std::string>& sequences_;
    bool ignore_mods_;

    explicit HasMatchingSequence(const set<std::string>& sequences, bool ignore_mods = false) : sequences_(sequences), ignore_mods_(ignore_mods)
    {
    }

    bool operator()(const PeptideHit& hit) const
    {
      const std::string& query = (ignore_mods_ ? hit.getSequence().toUnmodifiedString() : hit.getSequence().toString());
      return (sequences_.contains(query));
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


  void IDFilter::extractPeptideSequences(const PeptideIdentificationList& peptides, set<std::string>& sequences, bool ignore_mods)
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

  map<std::string, vector<ProteinHit>> IDFilter::extractUnassignedProteins(ConsensusMap& cmap)
  {
    map<std::string, vector<ProteinHit>> result;
    edit<ConsensusMap>(cmap, [&](ConsensusMap& map) {
      auto& data = map.getIdentificationData();
      // the accessions that peptides in features refer to, per protein run:
      const auto referenced = referencedAccessions(map, false);
      // every protein run is listed (by its legacy identifier), also without unassigned proteins:
      for (const auto& result_item : data.getInferenceResults())
      {
        result[result_item.proteins.getIdentifier().empty() ? result_item.identifier : result_item.proteins.getIdentifier()];
      }
      for (const auto& run : data.getRuns())
      {
        if (proteinResult(data, run) == nullptr) result[IdentificationDataAdapter::legacyIdentifier(run)];
      }
      keepProteins(data, [&](const std::string& key, const std::string& accession) {
        const auto found = referenced.find(key);
        return found != referenced.end() && found->second.contains(accession);
      }, &result);
    });
    return result;
  }

  void IDFilter::removeUnreferencedProteins(ConsensusMap& cmap, bool include_unassigned)
  {
    edit<ConsensusMap>(cmap, [&](ConsensusMap& map) {
      const auto referenced = referencedAccessions(map, include_unassigned);
      keepProteins(map.getIdentificationData(), [&](const std::string& key, const std::string& accession) {
        const auto found = referenced.find(key);
        return found != referenced.end() && found->second.contains(accession);
      });
    });
  }

  void IDFilter::removeUnreferencedProteins(ProteinIdentification& proteins, const PeptideIdentificationList& peptides)
  {
    auto run_to_accessions = collectReferencedAccessions(peptides);
    const unordered_set<std::string>& accessions = run_to_accessions[proteins.getIdentifier()];
    HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
    keepMatchingItems(proteins.getHits(), acc_filter);
  }

  void IDFilter::removeUnreferencedProteins(vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides)
  {
    auto run_to_accessions = collectReferencedAccessions(peptides);
    for (ProteinIdentification& prot : proteins)
    {
      const unordered_set<std::string>& accessions = run_to_accessions[prot.getIdentifier()];
      HasMatchingAccessionUnordered<ProteinHit> acc_filter(accessions);
      keepMatchingItems(prot.getHits(), acc_filter);
    }
  }

  void IDFilter::removeDanglingProteinReferences(ConsensusMap& cmap, bool remove_peptides_without_reference)
  {
    edit<ConsensusMap>(cmap, [&](ConsensusMap& map) {
      const auto& data = map.getIdentificationData();
      filterReferences(map, [&](const IdentificationData::Run& run) { return proteinAccessions(data, run); }, remove_peptides_without_reference);
    });
  }

  void IDFilter::removeDanglingProteinReferences(ConsensusMap& cmap, const ProteinIdentification& ref_run, bool remove_peptides_without_reference)
  {
    unordered_set<std::string> accessions;
    for (const ProteinHit& hit : ref_run.getHits())
    {
      accessions.insert(hit.getAccession());
    }
    edit<ConsensusMap>(cmap, [&](ConsensusMap& map) {
      filterReferences(map, [&](const IdentificationData::Run&) { return accessions; }, remove_peptides_without_reference);
    });
  }

  void IDFilter::keepNBestPeptideHits(FeatureMap& map, Size n)
  {
    edit<FeatureMap>(map, [&](FeatureMap& m) { keepNBest(m, n); });
  }

  void IDFilter::keepNBestPeptideHits(ConsensusMap& map, Size n)
  {
    edit<ConsensusMap>(map, [&](ConsensusMap& m) { keepNBest(m, n); });
  }

  void IDFilter::filterHitsByScore(FeatureMap& map, double threshold_score)
  {
    edit<FeatureMap>(map, [&](FeatureMap& m) { filterByScore(m, threshold_score); });
  }

  void IDFilter::filterHitsByScore(ConsensusMap& map, double threshold_score)
  {
    edit<ConsensusMap>(map, [&](ConsensusMap& m) { filterByScore(m, threshold_score); });
  }

  void IDFilter::keepUniquePeptidesPerProtein(FeatureMap& map)
  {
    edit<FeatureMap>(map, [](FeatureMap& m) { keepUnique(m); });
  }

  void IDFilter::keepUniquePeptidesPerProtein(ConsensusMap& map)
  {
    edit<ConsensusMap>(map, [](ConsensusMap& m) { keepUnique(m); });
  }

  void IDFilter::removeEmptyIdentifications(FeatureMap& map)
  {
    edit<FeatureMap>(map, [](FeatureMap& m) { removeEmpty(m); });
  }

  void IDFilter::removeEmptyIdentifications(ConsensusMap& map)
  {
    edit<ConsensusMap>(map, [](ConsensusMap& m) { removeEmpty(m); });
  }

  void IDFilter::annotateBestPerPeptidePerRun(FeatureMap& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
  {
    edit<FeatureMap>(map, [&](FeatureMap& m) { annotateBestNative(m, ignore_mods, ignore_charges, nr_best_spectrum); });
  }

  void IDFilter::annotateBestPerPeptidePerRun(ConsensusMap& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
  {
    edit<ConsensusMap>(map, [&](ConsensusMap& m) { annotateBestNative(m, ignore_mods, ignore_charges, nr_best_spectrum); });
  }

  void IDFilter::keepBestPerPeptidePerRun(FeatureMap& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
  {
    edit<FeatureMap>(map, [&](FeatureMap& m) { keepBestNative(m, ignore_mods, ignore_charges, nr_best_spectrum); });
  }

  void IDFilter::keepBestPerPeptidePerRun(ConsensusMap& map, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
  {
    edit<ConsensusMap>(map, [&](ConsensusMap& m) { keepBestNative(m, ignore_mods, ignore_charges, nr_best_spectrum); });
  }

  void IDFilter::removeDanglingProteinReferences(PeptideIdentificationList& peptides, const vector<ProteinIdentification>& proteins, bool remove_peptides_without_reference)
  {
    auto run_to_accessions = collectProteinAccessions(proteins);
    for (PeptideIdentification& pep : peptides)
    {
      filterPeptideReferences(pep, run_to_accessions[pep.getIdentifier()], remove_peptides_without_reference);
    }
  }


  bool IDFilter::updateProteinGroups(vector<ProteinIdentification::ProteinGroup>& groups, const vector<ProteinHit>& hits)
  {
    if (groups.empty())
      return true; // nothing to update

    // we'll do lots of look-ups, so use a suitable data structure:
    unordered_set<std::string> valid_accessions;
    for (const ProteinHit& hit : hits)
    {
      valid_accessions.insert(hit.getAccession());
    }

    bool valid = true;
    vector<ProteinIdentification::ProteinGroup> filtered_groups;
    for (ProteinIdentification::ProteinGroup& group : groups)
    {
      ProteinIdentification::ProteinGroup filtered;
      for (const std::string& acc : group.accessions)
      {
        if (valid_accessions.contains(acc))
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
        // Carry over the quantities attached by PeptideAndProteinQuant. They are indexed by
        // (fraction group, label) assay resp. by (file, channel), not by group member, so removing a
        // protein from the group does not invalidate them. Dropping them here silently discarded
        // protein abundances in every tool that filters after quantification.
        filtered.setFloatDataArrays(group.getFloatDataArrays());
        filtered.setStringDataArrays(group.getStringDataArrays());
        filtered.setIntegerDataArrays(group.getIntegerDataArrays());
        filtered_groups.push_back(std::move(filtered));
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
    unordered_set<std::string> valid_accessions;
    for (const auto& grp : groups)
    {
      valid_accessions.insert(grp.accessions.begin(), grp.accessions.end());
    }

    hits.erase(std::remove_if(hits.begin(), hits.end(), std::not_fn(HasMatchingAccessionUnordered<ProteinHit>(valid_accessions))), hits.end());
  }

  void IDFilter::keepBestPeptideHits(PeptideIdentificationList& peptides, bool strict)
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

  void IDFilter::filterPeptidesByLength(PeptideIdentificationList& peptides, Size min_length, Size max_length)
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

  void IDFilter::filterPeptidesByCharge(PeptideIdentificationList& peptides, Int min_charge, Int max_charge)
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


  void IDFilter::filterPeptidesByRT(PeptideIdentificationList& peptides, double min_rt, double max_rt)
  {
    struct HasRTInRange rt_filter(min_rt, max_rt);
    keepMatchingItems(peptides, rt_filter);
  }


  void IDFilter::filterPeptidesByMZ(PeptideIdentificationList& peptides, double min_mz, double max_mz)
  {
    struct HasMZInRange mz_filter(min_mz, max_mz);
    keepMatchingItems(peptides, mz_filter);
  }


  void IDFilter::filterPeptidesByMZError(PeptideIdentificationList& peptides, double mass_error, bool unit_ppm)
  {
    for (PeptideIdentification& pep : peptides)
    {
      struct HasLowMZError error_filter(pep.getMZ(), mass_error, unit_ppm);
      keepMatchingItems(pep.getHits(), error_filter);
    }
  }


  void IDFilter::filterPeptidesByRTPredictPValue(PeptideIdentificationList& peptides, const std::string& metavalue_key, double threshold)
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


  void IDFilter::removePeptidesWithMatchingModifications(PeptideIdentificationList& peptides, const set<std::string>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (PeptideIdentification& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), mod_filter);
    }
  }

  void IDFilter::removePeptidesWithMatchingRegEx(PeptideIdentificationList& peptides, const std::string& regex)
  {
    const std::regex re(regex);

    // true if regex matches to parts or entire unmodified sequence
    auto regex_matches = [&re](const PeptideHit& ph) -> bool { return std::regex_search(ph.getSequence().toUnmodifiedString(), re); };

    for (auto& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), regex_matches);
    }
  }

  void IDFilter::keepPeptidesWithMatchingModifications(PeptideIdentificationList& peptides, const set<std::string>& modifications)
  {
    struct HasMatchingModification mod_filter(modifications);
    for (PeptideIdentification& pep : peptides)
    {
      keepMatchingItems(pep.getHits(), mod_filter);
    }
  }


  void IDFilter::removePeptidesWithMatchingSequences(PeptideIdentificationList& peptides, const PeptideIdentificationList& bad_peptides, bool ignore_mods)
  {
    set<std::string> bad_seqs;
    extractPeptideSequences(bad_peptides, bad_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(bad_seqs, ignore_mods);
    for (PeptideIdentification& pep : peptides)
    {
      removeMatchingItems(pep.getHits(), seq_filter);
    }
  }


  void IDFilter::keepPeptidesWithMatchingSequences(PeptideIdentificationList& peptides, const PeptideIdentificationList& good_peptides, bool ignore_mods)
  {
    set<std::string> good_seqs;
    extractPeptideSequences(good_peptides, good_seqs, ignore_mods);
    struct HasMatchingSequence seq_filter(good_seqs, ignore_mods);
    for (PeptideIdentification& pep : peptides)
    {
      keepMatchingItems(pep.getHits(), seq_filter);
    }
  }


  void IDFilter::keepUniquePeptidesPerProtein(PeptideIdentificationList& peptides)
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
  void IDFilter::removeDuplicatePeptideHits(PeptideIdentificationList& peptides, bool seq_only)
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

  void IDFilter::keepNBestSpectra(PeptideIdentificationList& peptides, Size n)
  {
    std::string score_type;
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
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,std::string("PSM score types must be identical to allow proper filtering."));
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

  void IDFilter::keepBestMatchPerObservation(IdentificationData& data, const IdentificationData::ScoreDefinition& score)
  {
    data.validate();
    IdentificationData replacement(data);
    for (const auto& current : replacement.getRuns())
    {
      auto& run = replacement.getRun(current.getIdentifier());
      if (! run.getScoreDefinitions().empty()) run.retainBest(run.findScore(score), false);
    }
    data.swap(replacement);
  }

  void IDFilter::filterObservationMatchesByScore(IdentificationData& data, const IdentificationData::ScoreDefinition& score, double cutoff)
  {
    data.validate();
    IdentificationData replacement(data);
    for (const auto& current : replacement.getRuns())
    {
      auto& run = replacement.getRun(current.getIdentifier());
      if (run.getScoreDefinitions().empty()) continue;
      auto view = run.bindScore(run.findScore(score));
      run.filterMatches([&](const IdentificationData::Match& match) {
        auto value = view(match);
        return value && (score.higher_better ? *value >= cutoff : *value <= cutoff);
      });
    }
    data.swap(replacement);
  }

  void IDFilter::removeDecoys(IdentificationData& data)
  {
    using ID = IdentificationData;
    data.validate();
    ID replacement(data);
    for (const auto& current : replacement.getRuns())
    {
      auto& run = replacement.getRun(current.getIdentifier());
      std::set<std::pair<UInt32, std::string>> decoys;
      if (run.getDatabaseSequences())
      {
        auto sequences = *run.getDatabaseSequences();
        for (const auto& sequence : sequences)
          if (sequence.target_decoy == ID::TargetDecoy::DECOY) decoys.emplace(sequence.database.value, sequence.accession);
        std::erase_if(sequences, [](const auto& sequence) { return sequence.target_decoy == ID::TargetDecoy::DECOY; });
        run.setDatabaseSequences(std::move(sequences));
      }
      run.eraseMatches([](const auto& match) { return match.target_decoy == ID::TargetDecoy::DECOY; });
      run.transformMatches([&](ID::MatchData& match) {
        const auto removed = std::erase_if(match.sequence_evidence,
                                           [&](const auto& evidence) { return decoys.contains({evidence.database.value, evidence.accession}); });
        if (removed && match.target_decoy == ID::TargetDecoy::BOTH) match.target_decoy = ID::TargetDecoy::TARGET;
      });
    }
    data.swap(replacement);
  }

} // namespace OpenMS
