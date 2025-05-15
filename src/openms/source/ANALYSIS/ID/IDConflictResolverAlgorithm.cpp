// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

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
      vector<PeptideIdentification> & peptides,
      vector<PeptideIdentification> & removed,
      UInt64 uid)
  {
    if (peptides.empty()) { return; }

    for (PeptideIdentification & pep : peptides)
    {
      // sort hits
      pep.sort();
    }

    vector<PeptideIdentification>::iterator pos;
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
    vector<PeptideIdentification> & peptides, 
    vector<PeptideIdentification> & removed,
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

    vector<PeptideIdentification>::iterator pos;
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
    vector<PeptideIdentification>::iterator pos1p = pos + 1;
    for (auto it = pos1p; it != peptides.end(); ++it) // OMS_CODING_TEST_EXCLUDE
    {
      removed.push_back(*it);
    }

    // set best one to first position and shrink vector
    peptides[0] = *pos;
    peptides.resize(1);
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
