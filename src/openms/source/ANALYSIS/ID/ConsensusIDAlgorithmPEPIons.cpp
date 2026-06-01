// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmPEPIons::ConsensusIDAlgorithmPEPIons()
  {
    setName("ConsensusIDAlgorithmPEPIons"); // DefaultParamHandler

    defaults_.setValue("mass_tolerance", 0.5, "Maximum difference between fragment masses (in Da) for fragments to be considered 'shared' between peptides .");
    defaults_.setMinFloat("mass_tolerance", 0.0);    
    // is the "min_shared" parameter necessary/useful?
    defaults_.setValue("min_shared", 2, "The minimal number of 'shared' fragments (between two suggested peptides) that is necessary to evaluate the similarity based on shared peak count (SPC).");
    defaults_.setMinInt("min_shared", 1);

    defaultsToParam_();
  }


   void ConsensusIDAlgorithmPEPIons::updateMembers_()
  {
    ConsensusIDAlgorithmSimilarity::updateMembers_();

    // similarity scoring based on shared peak count:
    mass_tolerance_ = param_.getValue("mass_tolerance");
    min_shared_ = param_.getValue("min_shared");

    // new parameters may affect the similarity calculation, so clear cache:
    similarities_.clear();
  }


  double ConsensusIDAlgorithmPEPIons::getSimilarity_(AASequence seq1,
                                                     AASequence seq2)
  {
    if (seq1 == seq2) return 1.0;
    // order of sequences matters for cache look-up:
    if (seq2 < seq1) std::swap(seq1, seq2); // "operator>" not defined
    pair<AASequence, AASequence> seq_pair = make_pair(seq1, seq2);
    SimilarityCache::iterator pos = similarities_.find(seq_pair);
    if (pos != similarities_.end()) return pos->second; // score found in cache

    // compare b and y ion series of seq. 1 and seq. 2:
    vector<double> ions1(2 * seq1.size()), ions2(2 * seq2.size());
    // b ions, seq. 1:
    ions1[0] = seq1.getPrefix(1).getMonoWeight(); // includes N-terminal mods
    // y ions, seq. 1:
    ions1[seq1.size()] = seq1.getSuffix(1).getMonoWeight(); // inc. C-term. mods
    for (Size i = 1; i < seq1.size(); ++i)
    {
      ions1[i] = ions1[i - 1] + seq1[i].getMonoWeight();
      ions1[seq1.size() + i] = (ions1[seq1.size() + i - 1] + 
                                seq1[seq1.size() - i - 1].getMonoWeight());
    }
    // b ions, seq. 2:
    ions2[0] = seq2.getPrefix(1).getMonoWeight(); // includes N-terminal mods
    // y ions, seq. 2:
    ions2[seq2.size()] = seq2.getSuffix(1).getMonoWeight(); // inc. C-term. mods
    for (Size i = 1; i < seq2.size(); ++i)
    {
      ions2[i] = ions2[i - 1] + seq2[i].getMonoWeight();
      ions2[seq2.size() + i] = (ions2[seq2.size() + i - 1] + 
                                seq2[seq2.size() - i - 1].getMonoWeight());
    }

    // now compare fragment masses from both sequences to find best matches
    // within the allowed tolerance; note that:
    // 1. we can be more efficient than comparing "all against all"
    // 2. an ion from seq. 2 may be the best match for two (similar) ions from
    // seq. 1 - then we want to count that ion only once, not twice 
    sort(ions1.begin(), ions1.end());
    sort(ions2.begin(), ions2.end());
    set<double> matches; // each best-matching ion counts only once
    vector<double>::iterator start = ions2.begin();
    // for each fragment in seq. 1...
    for (vector<double>::iterator it1 = ions1.begin(); it1 != ions1.end();
         ++it1)
    {
      // ...find fragments from seq. 2 that are within the mass tolerance:
      vector<double>::iterator lower = lower_bound(start, ions2.end(),
                                                   *it1 - mass_tolerance_);
      if (lower == ions2.end()) break; // all values are too low
      vector<double>::iterator upper = upper_bound(lower, ions2.end(),
                                                   *it1 + mass_tolerance_);
      double best_match = 0.0, best_diff = mass_tolerance_ + 1.0;
      // find ion from seq. 2 (*it2) that is closest to ion from seq. 1 (*it1):
      for (vector<double>::iterator it2 = lower; it2 != upper; ++it2)
      {
        double diff = fabs(*it1 - *it2);
        if (diff < best_diff)
        {
          best_diff = diff;
          best_match = *it2;
        }
      }
      if (best_diff <= mass_tolerance_) matches.insert(best_match);
      
      start = lower; // "*it1" is increasing, so lower bounds can't get lower
    }

    double score_sim = 0.0;
    if (matches.size() >= min_shared_)
    {
      score_sim = matches.size() / float(min(ions1.size(), ions2.size()));
    }
    similarities_[seq_pair] = score_sim; // cache the similarity score

    return score_sim;
  }

} // namespace OpenMS
