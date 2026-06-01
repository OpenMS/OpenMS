// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

#include <limits>    // for std::numeric_limits
#include <set>

using namespace std;

namespace OpenMS
{
  void IDConflictResolverAlgorithm::resolve(FeatureMap & features, bool keep_matching)
  {
    resolveConflict_(features, keep_matching);
  }
  
  void IDConflictResolverAlgorithm::resolve(ConsensusMap & features, bool keep_matching)
  {
    resolveConflict_(features, keep_matching);
  }

  void IDConflictResolverAlgorithm::resolveAllHitRankAggregation(FeatureMap& features)
  {
    rankAggregation_(features);
  }

  void IDConflictResolverAlgorithm::resolveAllHitRankAggregation(ConsensusMap& features)
  {
    rankAggregation_(features);
  }
  
  void IDConflictResolverAlgorithm::resolveBetweenFeatures(FeatureMap & features)
  {
    resolveBetweenFeatures_(features);
  }
  
  void IDConflictResolverAlgorithm::resolveBetweenFeatures(ConsensusMap & features)
  {
    resolveBetweenFeatures_(features);
  }

  // static
  void IDConflictResolverAlgorithm::resolveConflictKeepMatching_(
      PeptideIdentificationList & peptides,
      PeptideIdentificationList & removed,
      UInt64 uid)
  {
    if (peptides.empty()) { return; }

    for (PeptideIdentification & pep : peptides)
    {
      // sort hits
      pep.sort();
    }

    PeptideIdentificationList::iterator pos;
    if (peptides[0].isHigherScoreBetter())     // find highest-scoring ID
    {
      pos = max_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }
    else  // find lowest-scoring ID
    {
      pos = min_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }

    const AASequence& best = (*pos).getHits()[0].getSequence();
    std::swap(*peptides.begin(), *pos); // put best on first position

    // filter for matching PEP Sequence and move to unassigned/removed
    for (auto it = ++peptides.begin(); it != peptides.end();)
    {
      auto& hits = it->getHits();
      auto hit = hits.begin();
      for (; hit != hits.end(); ++hit)
      {
        if (hit->getSequence() == best)
        {
          break;
        }
      }
      if (hit != hits.end()) // found sequence
      {
        hits[0] = *hit; // put the match on first place
        hits.resize(1); // remove rest
        ++it;
      }
      else // not found
      {
        // annotate feature_id for later reference
        it->setMetaValue("feature_id", String(uid));
        // move to "removed" vector
        removed.push_back(std::move(*it));
        // erase and update iterator
        it = peptides.erase(it);
      }
    }
  }

  // static
  void IDConflictResolverAlgorithm::resolveConflict_(
    PeptideIdentificationList & peptides, 
    PeptideIdentificationList & removed,
    UInt64 uid)
  {
    if (peptides.empty()) { return; }

    for (PeptideIdentification & pep : peptides)
    {
      // sort hits
      pep.sort();

      // remove all but the best hit
      if (!pep.getHits().empty())
      {
        vector<PeptideHit> best_hit(1, pep.getHits()[0]);
        pep.setHits(best_hit);
      }
      // annotate feature id
      pep.setMetaValue("feature_id", String(uid));
    }

    PeptideIdentificationList::iterator pos;
    if (peptides[0].isHigherScoreBetter())     // find highest-scoring ID
    {
      pos = max_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }
    else  // find lowest-scoring ID
    {
      pos = min_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }

    // copy conflicting ones left to best one
    for (auto it = peptides.begin(); it != pos; ++it)
    {
      removed.push_back(*it);
    }
     
    // copy conflicting ones right of best one
    PeptideIdentificationList::iterator pos1p = pos + 1;
    for (auto it = pos1p; it != peptides.end(); ++it) // OMS_CODING_TEST_EXCLUDE
    {
      removed.push_back(*it);
    }

    // set best one to first position and shrink vector
    peptides[0] = *pos;
    peptides.resize(1);
  }

  // static
  void IDConflictResolverAlgorithm::resolveAggregateConflict_(
      PeptideIdentificationList& peptides,
      PeptideIdentificationList& removed,
      UInt64 uid)
  {
    if (peptides.empty()) { return; }

    // Annotate all IDs with feature_id and sort hits (best first)
    for (PeptideIdentification& pep : peptides)
    {
      pep.setMetaValue("feature_id", String(uid));
      pep.sort();
    }

    // If only one ID, keep only the best hit and return
    if (peptides.size() == 1)
    {
      if (!peptides[0].getHits().empty())
      {
        vector<PeptideHit> best(1, peptides[0].getHits()[0]);
        peptides[0].setHits(best);
      }
      return;
    }

    // Find the maximum number of hits across all IDs (used as penalty)
    Size max_hits = 0;
    for (const PeptideIdentification& pep : peptides)
    {
      if (pep.getHits().size() > max_hits)
      {
        max_hits = pep.getHits().size();
      }
    }

    if (max_hits == 0)
    {
      // All IDs are empty - move all but the first to removed
      for (auto it = ++peptides.begin(); it != peptides.end(); ++it)
      {
        removed.push_back(*it);
      }
      peptides.resize(1);
      return;
    }

    Size n_runs = peptides.size();

    // For each unique sequence, accumulate ranks across all IDs
    // rank = 0-based position in sorted hits; if missing from an ID, add max_hits as penalty
    map<AASequence, double> rank_sums;
    map<AASequence, Size> run_counts;

    for (const PeptideIdentification& pep : peptides)
    {
      const vector<PeptideHit>& hits = pep.getHits();
      set<AASequence> seen_in_run;

      for (Size j = 0; j < hits.size(); ++j)
      {
        const AASequence& seq = hits[j].getSequence();
        if (seen_in_run.count(seq) == 0)
        {
          // First occurrence of this sequence in this ID: use its rank
          rank_sums[seq] += static_cast<double>(j);
          run_counts[seq]++;
          seen_in_run.insert(seq);
        }
      }
    }

    // Compute aggregate score for each sequence and find the best one
    // score = 1.0 - (rank_sum + penalty_for_missing) / (max_hits * n_runs)
    // Higher score is better.
    AASequence best_seq;
    double best_agg_score = -1.0;

    for (const auto& [seq, rank_sum] : rank_sums)
    {
      Size runs_with_seq = run_counts.at(seq);
      double total_rank = rank_sum + static_cast<double>((n_runs - runs_with_seq) * max_hits);
      double agg_score = 1.0 - total_rank / (static_cast<double>(max_hits) * static_cast<double>(n_runs));

      if (agg_score > best_agg_score)
      {
        best_agg_score = agg_score;
        best_seq = seq;
      }
    }

    // Among the IDs that contain best_seq, find the one with the best original score for that hit
    bool higher_better = peptides[0].isHigherScoreBetter();
    PeptideIdentificationList::iterator best_id_it = peptides.end();
    double best_original_score = higher_better ? -numeric_limits<double>::infinity()
                                               :  numeric_limits<double>::infinity();

    for (auto it = peptides.begin(); it != peptides.end(); ++it)
    {
      for (const PeptideHit& hit : it->getHits())
      {
        if (hit.getSequence() == best_seq)
        {
          if (best_id_it == peptides.end() ||
              (higher_better  && hit.getScore() > best_original_score) ||
              (!higher_better && hit.getScore() < best_original_score))
          {
            best_id_it = it;
            best_original_score = hit.getScore();
          }
          break; // only consider the first (best-ranked) occurrence in this ID
        }
      }
    }

    if (best_id_it == peptides.end())
    {
      best_id_it = peptides.begin(); // fallback (should never happen)
    }

    // Move all IDs except the best one to removed (reduced to 1 hit each for consistency)
    for (auto it = peptides.begin(); it != peptides.end(); ++it)
    {
      if (it != best_id_it)
      {
        // Reduce to best hit before moving to removed
        if (!it->getHits().empty())
        {
          vector<PeptideHit> best_hit_only(1, it->getHits()[0]);
          it->setHits(best_hit_only);
        }
        removed.push_back(*it);
      }
    }

    // Reduce the best ID to just the winning hit and shrink the list to one element
    const vector<PeptideHit>& best_hits = best_id_it->getHits();
    for (const PeptideHit& hit : best_hits)
    {
      if (hit.getSequence() == best_seq)
      {
        best_id_it->setHits({hit});
        break;
      }
    }
    PeptideIdentificationList result;
    result.push_back(*best_id_it);
    peptides = result;
  }

  // static
  bool IDConflictResolverAlgorithm::compareIDsSmallerScores_(const PeptideIdentification & left, const PeptideIdentification & right)
  {
    // if any of them is empty, the other is considered "greater"
    // independent of the score in the first hit
    if (left.getHits().empty() || right.getHits().empty()) 
    { // also: for strict weak ordering, comp(x,x) needs to be false
      return left.getHits().size() < right.getHits().size();
    }

    return left.getHits()[0].getScore() < right.getHits()[0].getScore();
  }
}

/// @endcond
