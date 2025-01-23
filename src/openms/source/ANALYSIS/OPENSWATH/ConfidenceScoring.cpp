// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h> 

#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <numeric> // for "accumulate"
#include <ctime> // for "time" (random number seed)
#include <random>
#include <map>

using namespace std;

namespace OpenMS
{
    /// Mapping: Q3 m/z <-> transition intensity (maybe not unique!)
    typedef boost::bimap<double, boost::bimaps::multiset_of<double> > 
    BimapType;

    ConfidenceScoring::ConfidenceScoring(bool test_mode_)
    {
      if (!test_mode_) shuffler_ = Math::RandomShuffler(0);
      else shuffler_ = Math::RandomShuffler(time(nullptr));// seed with current time
    }

    /// Randomize the list of decoy indexes
    void ConfidenceScoring::chooseDecoys_()
    {
      if (n_decoys_ == 0) return; // list is already initialized
      // somewhat inefficient to shuffle the whole list when we only need a random
      // sample, but easy to do...
      shuffler_.portable_random_shuffle(decoy_index_.begin(), decoy_index_.end());
    }

    // double rmsd_(DoubleList x, DoubleList y)
    // {
    //   double sum_of_squares = 0;
    //   for (Size i = 0; i < x.size(); i++)
    //   {
    //     double diff = x[i] - y[i];
    //     sum_of_squares += diff * diff;
    //   }
    //   return sqrt(sum_of_squares / x.size());
    // }

    /// Manhattan distance
    double ConfidenceScoring::manhattanDist_(DoubleList x, DoubleList y)
    {
      double sum = 0;
      for (Size i = 0; i < x.size(); ++i)
      {
        sum += fabs(x[i] - y[i]);
      }
      return sum;
    }

    /// Get the retention time of an assay
    double ConfidenceScoring::getAssayRT_(const TargetedExperiment::Peptide& assay)
    {
      OPENMS_PRECONDITION(assay.hasRetentionTime(), "More than zero RTs needed")
      return assay.getRetentionTime();
    }


    /// Extract the @p n_transitions highest intensities from @p intensity_map,
    /// store them in @p intensities
    void extractIntensities_(BimapType& intensity_map, Size n_transitions,
                             DoubleList& intensities)
    {
      // keep only as many transitions as needed, remove those with lowest
      // intensities:
      if (n_transitions > 0)
      {
        // use "Int" instead of "Size" to prevent overflows:
        Int diff = intensity_map.size() - n_transitions;
        for (Size i = 0; Int(i) < diff; ++i)
        {
          intensity_map.right.erase(intensity_map.right.begin());
        }
      }
      // fill output list ordered by m/z:
      intensities.clear();
      for (BimapType::left_map::iterator int_it = intensity_map.left.begin();
           int_it != intensity_map.left.end(); ++int_it)
      {
        // missing values might be "-1"
        intensities.push_back(max(0.0, int_it->second));
      }
    }

    /// Score the assay @p assay against feature data (@p feature_rt,
    /// @p feature_intensities), optionally using only the specified transitions
    /// (@p transition_ids)
    double ConfidenceScoring::scoreAssay_(const TargetedExperiment::Peptide& assay, 
                           double feature_rt, DoubleList& feature_intensities,
                           const std::set<String>& transition_ids)
    {
      // compute RT difference:
      double assay_rt = rt_norm_(getAssayRT_(assay));
      double diff_rt = assay_rt - feature_rt;

      // collect transition intensities:
      BimapType intensity_map;
      for (IntList::iterator trans_it = transition_map_[assay.id].begin();
           trans_it != transition_map_[assay.id].end(); ++trans_it)
      {
        const ReactionMonitoringTransition& transition = 
          library_.getTransitions()[*trans_it];
        // for the "true" assay, we need to choose the same transitions as for the
        // feature:
        if (!transition_ids.empty() && 
            (transition_ids.count(transition.getNativeID()) == 0)) continue;
        // seems like Boost's Bimap doesn't support "operator[]"...
        intensity_map.left.insert(make_pair(transition.getProductMZ(), 
                                            transition.getLibraryIntensity()));
      }
      DoubleList assay_intensities;
      extractIntensities_(intensity_map, feature_intensities.size(), 
                          assay_intensities);

      if (feature_intensities.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Feature intensities were empty - please provide feature subordinate with intensities");
      }
      if (feature_intensities.size()!=assay_intensities.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Did not find a feature for each assay provided - each feature needs "
          "to have n subordinates with the meta-value 'native_id' set to the corresponding transition.");
      }

      // compute intensity distance:
      OpenSwath::Scoring::normalize_sum(&feature_intensities[0],
                                        boost::numeric_cast<int>(feature_intensities.size()));
      OpenSwath::Scoring::normalize_sum(&assay_intensities[0],
                                        boost::numeric_cast<int>(assay_intensities.size()));
      double dist_int = manhattanDist_(feature_intensities, 
                                           assay_intensities);

      double score = glm_(diff_rt, dist_int);

      OPENMS_LOG_DEBUG << "\ndelta_RT:  " << fabs(diff_rt)
                << "\ndist_int:  " << dist_int
                << "\nGLM_score: " << score << endl;

      return score;
    }

    /// Score a feature
    void ConfidenceScoring::scoreFeature_(Feature& feature)
    {
      // extract predictors from feature:
      double feature_rt = rt_norm_(rt_trafo_.apply(feature.getRT()));
      BimapType intensity_map;
      // for the "true" assay, we need to make sure we compare based on the same
      // transitions, so keep track of them:
      std::map<double, String> trans_id_map; // Q3 m/z -> transition ID
      for (vector<Feature>::iterator sub_it = feature.getSubordinates().begin();
           sub_it != feature.getSubordinates().end(); ++sub_it)
      {
        // seems like Boost's Bimap doesn't support "operator[]"...
        intensity_map.left.insert(make_pair(sub_it->getMZ(), 
                                            sub_it->getIntensity()));
        trans_id_map[sub_it->getMZ()] = sub_it->getMetaValue("native_id");
      }
      DoubleList feature_intensities;
      extractIntensities_(intensity_map, n_transitions_, feature_intensities);
      if ((n_transitions_ > 0) && (feature_intensities.size() < n_transitions_))
      {
        OPENMS_LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
                 << "' contains only " << feature_intensities.size()
                 << " transitions." << endl;
      }
      // "intensity_map" now only contains the transitions we need later:
      std::set<String> transition_ids;
      for (BimapType::left_map::iterator int_it = intensity_map.left.begin();
           int_it != intensity_map.left.end(); ++int_it)
      {
        transition_ids.insert(trans_id_map[int_it->first]);
      }

      DoubleList scores; // "true" score is in "scores[0]", decoy scores follow

      // compare to "true" assay:
      String true_id = feature.getMetaValue("PeptideRef");
      OPENMS_LOG_DEBUG << "True assay (ID '" << true_id << "')" << endl;
      if (true_id.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Feature does not contain meta value 'PeptideRef' (reference to assay)");
      }
      scores.push_back(scoreAssay_(library_.getPeptideByRef(true_id), feature_rt,
                            feature_intensities, transition_ids));

      // compare to decoy assays:
      chooseDecoys_();
      Size counter = 0;
      for (IntList::iterator decoy_it = decoy_index_.begin(); 
           decoy_it != decoy_index_.end(); ++decoy_it)
      {
        const TargetedExperiment::Peptide& decoy_assay = 
          library_.getPeptides()[*decoy_it];

        // skip the "true" assay and assays with too few transitions:
        // TODO: maybe add an option to include assays with too few transitions?
        if ((decoy_assay.id == true_id) || 
            (transition_map_[decoy_assay.id].size() < feature_intensities.size()))
        {
          continue;
        }
        OPENMS_LOG_DEBUG << "Decoy assay " << scores.size() << " (ID '" << decoy_assay.id
                  << "')" << endl;

        scores.push_back(scoreAssay_(decoy_assay, feature_rt, feature_intensities));

        if ((n_decoys_ > 0) && (++counter >= n_decoys_)) break; // enough decoys
      }
      
      Size n_scores = scores.size();
      if (n_scores - 1 < n_decoys_)
      {
        OPENMS_LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
                 << "': Couldn't find enough decoy assays with at least "
                 << feature_intensities.size() << " transitions. "
                 << "Scoring based on " << n_scores - 1 << " decoys." << endl;
      }
      // TODO: this warning may trigger for every feature and get annoying
      if ((n_decoys_ == 0) && (n_scores < library_.getPeptides().size()))
      {
        OPENMS_LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
                 << "': Skipped some decoy assays with fewer than " 
                 << feature_intensities.size() << " transitions. "
                 << "Scoring based on " << n_scores - 1 << " decoys." << endl;
      }

      // count decoy scores that are greater than the "true" score:
      counter = 0;
      for (DoubleList::iterator it = ++scores.begin(); it != scores.end(); ++it)
      {
        if (*it > scores[0]) counter++;
      }

      // annotate feature:
      feature.setMetaValue("GLM_score", scores[0]);
      double local_fdr = counter / (n_scores - 1.0);
      feature.setMetaValue("local_FDR", local_fdr);
      feature.setOverallQuality(1.0 - local_fdr);
    }

}
