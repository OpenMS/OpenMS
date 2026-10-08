// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <algorithm>
#include <limits>    // for std::numeric_limits
#include <map>
#include <optional>
#include <set>
#include <type_traits>

using namespace std;

namespace OpenMS
{
  namespace
  {
    using ID = IdentificationData;

    /// An identification of a feature with its linked matches, best first
    struct Entry
    {
      ID::Run* run = nullptr;
      ID::QueryReference query;
      bool higher_better = true;
      /// Linked matches, sorted by primary score (best first, stable); the scores in the same order
      std::vector<const ID::Match*> matches;
      std::vector<double> scores;
      /// Matches that stay with the identification (default: all)
      std::vector<const ID::Match*> retained;
      /// Whether the identification stays with the feature
      bool keep = false;
    };

    std::vector<Entry> entriesOf(const BaseFeature& feature, ID& data)
    {
      std::vector<Entry> entries;
      for (const auto& linked : feature.getLinkedIdentifications(data))
      {
        Entry entry;
        entry.run = data.findRunByUuid(linked.run->getUuid());
        entry.query = {linked.run->getUuid(), linked.query->getId()};
        const auto primary = linked.run->getPrimaryScore();
        entry.higher_better = primary ? linked.run->getScoreDefinition(*primary).higher_better : true;
        std::vector<std::pair<double, const ID::Match*>> scored;
        for (const auto* match : linked.matches)
        {
          const auto score = primary ? linked.run->getScore(match->getId(), *primary) : std::nullopt;
          scored.emplace_back(score.value_or(std::numeric_limits<double>::quiet_NaN()), match);
        }
        // like PeptideIdentification::sort()
        std::stable_sort(scored.begin(), scored.end(), [&](const auto& a, const auto& b) {
          return entry.higher_better ? a.first > b.first : a.first < b.first;
        });
        for (const auto& [score, match] : scored)
        {
          entry.matches.push_back(match);
          entry.scores.push_back(score);
        }
        entry.retained = entry.matches;
        entries.push_back(std::move(entry));
      }
      return entries;
    }

    /// The best identification: the one whose best match scores best, the first of equal ones; identifications without
    /// matches only if no identification has matches
    Size bestEntry(const std::vector<Entry>& entries)
    {
      const bool higher_better = entries[0].higher_better;
      Size best = 0;
      for (Size i = 1; i < entries.size(); ++i)
      {
        const auto& candidate = entries[i];
        const auto& current = entries[best];
        if (candidate.matches.empty()) continue;
        if (current.matches.empty() || (higher_better ? candidate.scores[0] > current.scores[0] : candidate.scores[0] < current.scores[0]))
        {
          best = i;
        }
      }
      return best;
    }

    /// The edits of a resolution, applied when all features are resolved
    struct Edits
    {
      /// meta value "feature_id" of identifications
      std::map<ID::QueryReference, std::string> feature_ids;
      /// matches removed from identifications (removed from the data unless another feature links them)
      std::set<ID::MatchReference> dropped;
    };

    /// Change the links of @p feature as decided in @p entries (Entry::keep, Entry::retained)
    void relink(BaseFeature& feature, const std::vector<Entry>& entries, Edits& edits)
    {
      for (const auto& entry : entries)
      {
        const auto& uuid = entry.query.run_uuid;
        for (const auto* match : entry.matches)
        {
          const ID::MatchReference reference {uuid, match->getId()};
          const bool retained = std::find(entry.retained.begin(), entry.retained.end(), match) != entry.retained.end();
          if (! retained) edits.dropped.insert(reference);
          if (! retained || ! entry.keep) feature.getIDMatches().erase(reference);
        }
        if (! entry.keep) feature.getIDQueries().erase(entry.query);
      }
    }

    enum class Mode
    {
      BEST,
      KEEP_MATCHING,
      RANK_AGGREGATION
    };

    /// Every identification keeps only its best match; the best identification stays, the others become unassigned
    void resolveBest(std::vector<Entry>& entries)
    {
      for (auto& entry : entries)
      {
        if (! entry.matches.empty()) entry.retained = {entry.matches[0]};
      }
      entries[bestEntry(entries)].keep = true;
    }

    /// The best identification stays with all its matches; another one stays (with just that match) if it has a match
    /// with the best identification's top sequence, else it becomes unassigned. Returns the identifications that leave.
    void resolveKeepMatching(std::vector<Entry>& entries, const std::string& uid, Edits& edits)
    {
      const Size best = bestEntry(entries);
      entries[best].keep = true;
      const ID::Match* top = entries[best].matches.empty() ? nullptr : entries[best].matches[0];
      for (Size i = 0; i < entries.size(); ++i)
      {
        if (i == best) continue;
        auto& entry = entries[i];
        const auto same = top ? std::find_if(entry.matches.begin(), entry.matches.end(),
                                             [&](const ID::Match* match) { return match->representation == top->representation; })
                              : entry.matches.end();
        if (same != entry.matches.end())
        {
          entry.retained = {*same};
          entry.keep = true;
        }
        else
        {
          // annotate feature_id for later reference
          edits.feature_ids[entry.query] = uid;
        }
      }
    }

    /// Aggregate the ranks of the sequences over the identifications (see resolveAllHitRankAggregation()); the
    /// identification with the best match of the winning sequence stays with that match, the others become unassigned
    /// with their best match.
    void resolveRankAggregation(std::vector<Entry>& entries)
    {
      if (entries.size() == 1)
      {
        resolveBest(entries);
        return;
      }
      Size max_hits = 0;
      for (const auto& entry : entries)
      {
        max_hits = std::max(max_hits, entry.matches.size());
      }
      if (max_hits == 0)
      {
        entries[0].keep = true;
        return;
      }
      const Size n_runs = entries.size();
      std::map<AASequence, double> rank_sums;
      std::map<AASequence, Size> run_counts;
      std::vector<std::vector<AASequence>> sequences(entries.size());
      for (Size e = 0; e < entries.size(); ++e)
      {
        std::set<AASequence> seen_in_run;
        for (Size j = 0; j < entries[e].matches.size(); ++j)
        {
          sequences[e].push_back(AASequence::fromString(entries[e].matches[j]->representation));
          const auto& seq = sequences[e].back();
          if (seen_in_run.insert(seq).second)
          {
            rank_sums[seq] += static_cast<double>(j);
            run_counts[seq]++;
          }
        }
      }
      AASequence best_seq;
      double best_agg_score = -1.0;
      for (const auto& [seq, rank_sum] : rank_sums)
      {
        const Size runs_with_seq = run_counts.at(seq);
        const double total_rank = rank_sum + static_cast<double>((n_runs - runs_with_seq) * max_hits);
        const double agg_score = 1.0 - total_rank / (static_cast<double>(max_hits) * static_cast<double>(n_runs));
        if (agg_score > best_agg_score)
        {
          best_agg_score = agg_score;
          best_seq = seq;
        }
      }
      // among the identifications with best_seq, the one with the best score for it (its first occurrence)
      const bool higher_better = entries[0].higher_better;
      std::optional<Size> best_entry;
      Size best_match = 0;
      double best_original_score = 0.0;
      for (Size e = 0; e < entries.size(); ++e)
      {
        const auto it = std::find(sequences[e].begin(), sequences[e].end(), best_seq);
        if (it == sequences[e].end()) continue;
        const Size j = static_cast<Size>(it - sequences[e].begin());
        const double score = entries[e].scores[j];
        if (! best_entry || (higher_better && score > best_original_score) || (! higher_better && score < best_original_score))
        {
          best_entry = e;
          best_match = j;
          best_original_score = score;
        }
      }
      const Size best = best_entry.value_or(0);
      for (Size e = 0; e < entries.size(); ++e)
      {
        if (e == best) continue;
        if (! entries[e].matches.empty()) entries[e].retained = {entries[e].matches[0]};
      }
      entries[best].keep = true;
      if (best_entry) entries[best].retained = {entries[best].matches[best_match]};
    }

    /// Set the meta values and remove the dropped matches that no feature links
    template<class Map>
    void applyEdits(Map& map, const Edits& edits)
    {
      auto& data = map.getIdentificationData();
      for (const auto& [query, value] : edits.feature_ids)
      {
        auto* run = data.findRunByUuid(query.run_uuid);
        ID::Observation observation = run->getIdentification(query.query);
        observation.setMetaValue("feature_id", value);
        run->replaceObservation(query.query, observation);
      }
      if (edits.dropped.empty()) return;
      std::set<ID::MatchReference> linked;
      const auto collect = [&](const auto& self, const auto& feature) -> void {
        linked.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
          for (const auto& subordinate : feature.getSubordinates())
            self(self, subordinate);
      };
      for (const auto& feature : map)
      {
        collect(collect, feature);
      }
      for (const auto& current : data.getRuns())
      {
        auto& run = data.getRun(current.getIdentifier());
        run.eraseMatches([&](const ID::Match& match) {
          const ID::MatchReference reference {run.getUuid(), match.getId()};
          return edits.dropped.contains(reference) && ! linked.contains(reference);
        }, true);
      }
    }

    void exportIDs(FeatureMap& map)
    {
      IdentificationDataConverter::exportFeatureIDs(map);
    }
    void exportIDs(ConsensusMap& map)
    {
      IdentificationDataConverter::exportConsensusIDs(map);
    }

    template<class Map>
    void resolveMap(Map& map, Mode mode)
    {
      const bool converted = IdentificationDataConverter::moveToIdentificationData(map);
      auto& data = map.getIdentificationData();
      Edits edits;
      // annotate as not part of the resolution
      for (const auto& entry : map.getUnassignedIdentifications())
      {
        edits.feature_ids[{entry.run->getUuid(), entry.query->getId()}] = "not mapped"; // not mapped to a feature
      }
      for (auto& feature : map)
      {
        const std::string uid = StringUtils::toStr(feature.getUniqueId());
        feature.setMetaValue("feature_id", uid);
        auto entries = entriesOf(feature, data);
        if (entries.empty()) continue;
        if (mode != Mode::KEEP_MATCHING)
        {
          for (const auto& entry : entries)
          {
            edits.feature_ids[entry.query] = uid;
          }
        }
        switch (mode)
        {
          case Mode::BEST:
            resolveBest(entries);
            break;
          case Mode::KEEP_MATCHING:
            resolveKeepMatching(entries, uid, edits);
            break;
          case Mode::RANK_AGGREGATION:
            resolveRankAggregation(entries);
            break;
        }
        relink(feature, entries, edits);
      }
      applyEdits(map, edits);
      if (converted) exportIDs(map);
    }

    template<class Map>
    void resolveBetweenFeaturesMap(Map& map)
    {
      const bool converted = IdentificationDataConverter::moveToIdentificationData(map);
      auto& data = map.getIdentificationData();
      // the feature with the highest intensity for each (charge, sequence of the best match)
      std::map<std::pair<Int, AASequence>, typename Map::value_type*> feature_set;
      const auto unassign = [&](typename Map::value_type& feature) {
        feature.getIDQueries().clear();
        feature.getIDMatches().clear();
      };
      for (auto& element : map)
      {
        const auto entries = entriesOf(element, data);
        if (entries.empty()) continue;
        if (entries.size() != 1)
        {
          // Should never happen. In IDConflictResolverAlgorithm TOPP tool
          // IDConflictResolverAlgorithm::resolve() is called before IDConflictResolverAlgorithm::resolveBetweenFeatures().
          throw Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Feature does contain multiple identifications.");
        }
        if (entries[0].matches.empty()) continue;
        const auto key = std::make_pair(element.getCharge(), AASequence::fromString(entries[0].matches[0]->representation));
        const auto feature_in_set = feature_set.find(key);
        if (feature_in_set == feature_set.end())
        {
          feature_set[key] = &element;
        }
        else if (feature_in_set->second->getIntensity() < element.getIntensity())
        {
          // the annotations of the old, less intense feature become unassigned
          unassign(*feature_in_set->second);
          feature_in_set->second = &element;
        }
        else
        {
          unassign(element);
        }
      }
      if (converted) exportIDs(map);
    }
  } // namespace

  void IDConflictResolverAlgorithm::resolve(FeatureMap & features, bool keep_matching)
  {
    resolveMap(features, keep_matching ? Mode::KEEP_MATCHING : Mode::BEST);
  }

  void IDConflictResolverAlgorithm::resolve(ConsensusMap & features, bool keep_matching)
  {
    resolveMap(features, keep_matching ? Mode::KEEP_MATCHING : Mode::BEST);
  }

  void IDConflictResolverAlgorithm::resolveAllHitRankAggregation(FeatureMap& features)
  {
    resolveMap(features, Mode::RANK_AGGREGATION);
  }

  void IDConflictResolverAlgorithm::resolveAllHitRankAggregation(ConsensusMap& features)
  {
    resolveMap(features, Mode::RANK_AGGREGATION);
  }

  void IDConflictResolverAlgorithm::resolveBetweenFeatures(FeatureMap & features)
  {
    resolveBetweenFeaturesMap(features);
  }

  void IDConflictResolverAlgorithm::resolveBetweenFeatures(ConsensusMap & features)
  {
    resolveBetweenFeaturesMap(features);
  }

  IDConflictResolverAlgorithm::UnresolvedIdentifications
  IDConflictResolverAlgorithm::reduceToOnePerSpectrum(PeptideIdentificationList& ids)
  {
    UnresolvedIdentifications report;

    // One keyable entry per identification. The spectrum reference is a metavalue, so
    // getSpectrumReference() builds a string on every call - read it once here rather than
    // from inside the comparator. The sequence is referenced, never copied.
    struct Entry
    {
      const std::string* reference;
      const AASequence* sequence;
      Int charge;
      Size index;
    };
    std::vector<std::string> references(ids.size());
    std::vector<Entry> entries;
    entries.reserve(ids.size());

    for (Size i = 0; i < ids.size(); ++i)
    {
      const PeptideIdentification& id = ids[i];
      if (id.getHits().empty())
      {
        continue; // nothing to key on, and nothing that could be double-counted
      }
      references[i] = id.getSpectrumReference();
      if (references[i].empty())
      {
        ++report.without_spectrum_reference;
        continue;
      }
      const PeptideHit& hit = id.getHits().front();
      entries.push_back({&references[i], &hit.getSequence(), hit.getCharge(), i});
    }

    // Sort by the key, then by input position so the surviving member of a group is
    // deterministic when scores tie.
    std::sort(entries.begin(), entries.end(), [](const Entry& left, const Entry& right)
    {
      if (*left.reference != *right.reference) { return *left.reference < *right.reference; }
      if (*left.sequence != *right.sequence)   { return *left.sequence < *right.sequence; }
      if (left.charge != right.charge)         { return left.charge < right.charge; }
      return left.index < right.index;
    });

    std::vector<bool> remove(ids.size(), false);

    for (Size begin = 0; begin < entries.size();)
    {
      // One spectrum's entries: [begin, spectrum_end).
      Size spectrum_end = begin;
      while (spectrum_end < entries.size()
             && *entries[spectrum_end].reference == *entries[begin].reference)
      {
        ++spectrum_end;
      }

      Size surviving_here = 0;
      for (Size group = begin; group < spectrum_end;)
      {
        // One (reference, sequence, charge) group: [group, group_end).
        // Sequence membership is tested with the ordering used to sort, not with operator==:
        // the two disagree on what identifies a residue (operator< compares the one-letter code,
        // operator== the interned Residue*), so a pair that operator< calls equivalent could be
        // split apart here and its duplicate silently missed. Since the range is sorted
        // ascending, "not less than its group leader" is exactly "sorts equal to it".
        Size group_end = group;
        while (group_end < spectrum_end
               && !(*entries[group].sequence < *entries[group_end].sequence)
               && entries[group_end].charge == entries[group].charge)
        {
          ++group_end;
        }

        if (group_end - group == 1) { ++surviving_here; group = group_end; continue; }

        // "Best" is only defined if the group agrees on which direction is better. It normally
        // does - a caller that got here has one score type - but a disagreement would otherwise
        // be resolved by a coin flip, silently discarding a measurement.
        const bool higher_better = ids[entries[group].index].isHigherScoreBetter();
        bool consistent = true;
        for (Size i = group + 1; i < group_end; ++i)
        {
          if (ids[entries[i].index].isHigherScoreBetter() != higher_better) { consistent = false; break; }
        }
        if (!consistent)
        {
          ++report.inconsistent_score_direction;
          // Counted as ONE, like the reduced branch: this group is a single peptidoform that
          // could not be reduced, not evidence that the spectrum carries several. Counting its
          // members individually would push multiply_identified_spectra above 1 and report the
          // spectrum as chimeric, which is the opposite of what happened.
          ++surviving_here;
          group = group_end;
          continue;
        }

        Size best = group;
        for (Size i = group + 1; i < group_end; ++i)
        {
          const double score = ids[entries[i].index].getHits().front().getScore();
          const double best_score = ids[entries[best].index].getHits().front().getScore();
          if (higher_better ? (score > best_score) : (score < best_score)) { best = i; }
        }

        for (Size i = group; i < group_end; ++i)
        {
          if (i != best) { remove[entries[i].index] = true; }
        }
        report.removed += group_end - group - 1;
        ++surviving_here;

        if (report.example.empty())
        {
          report.example = *entries[group].reference + " / " + entries[group].sequence->toString()
                         + " / charge " + StringUtils::toStr(entries[group].charge);
        }

        group = group_end;
      }

      if (surviving_here > 1) { ++report.multiply_identified_spectra; }
      begin = spectrum_end;
    }

    if (report.removed > 0)
    {
      auto& data = ids.getData();
      Size out = 0;
      for (Size i = 0; i < data.size(); ++i)
      {
        if (remove[i]) { continue; }
        if (out != i) { data[out] = std::move(data[i]); }
        ++out;
      }
      data.resize(out);
    }

    return report;
  }

}

/// @endcond
