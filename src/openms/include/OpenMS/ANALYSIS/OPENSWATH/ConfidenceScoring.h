// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <cmath> // for "exp"
#include <limits> // for "infinity"
#include <map>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <OpenMS/MATH/MathFunctions.h>

namespace OpenMS
{

  class OPENMS_DLLAPI ConfidenceScoring :
      public ProgressLogger
  {
  public:

      /// Constructor
      explicit ConfidenceScoring(bool test_mode_ = false);

      ~ConfidenceScoring() override {}

  protected:

      /// Binomial GLM
      struct GLM_
      {
        double intercept;
        double rt_coef;
        double int_coef;

        double operator()(double diff_rt, double dist_int) const
        {
          double lm = intercept + rt_coef * diff_rt * diff_rt + 
            int_coef * dist_int;
          return 1.0 / (1.0 + exp(-lm));
        }
      } glm_;

      /// Helper for RT normalization (range 0-100)
      struct RTNorm_
      {
        double min_rt;
        double max_rt;
        
        double operator()(double rt) const
        {
          return (rt - min_rt) / (max_rt - min_rt) * 100;
        }
      } rt_norm_;

      TargetedExperiment library_; ///< assay library

      IntList decoy_index_; ///< indexes of assays to use as decoys

      Size n_decoys_; ///< number of decoys to use (per feature/true assay)

      std::map<String, IntList> transition_map_; ///< assay (ID) -> transitions (indexes)

      Size n_transitions_; ///< number of transitions to consider

      /// RT transformation to map measured RTs to assay RTs
      TransformationDescription rt_trafo_;

      Math::RandomShuffler shuffler_; ///< random shuffler for container

      /// Randomize the list of decoy indexes
      void chooseDecoys_();

      /// Manhattan distance
      double manhattanDist_(DoubleList x, DoubleList y);

      /// Get the retention time of an assay
      double getAssayRT_(const TargetedExperiment::Peptide& assay);

      /// Score the assay @p assay against feature data (@p feature_rt,
      /// @p feature_intensities), optionally using only the specified transitions
      /// (@p transition_ids)
      double scoreAssay_(const TargetedExperiment::Peptide& assay, 
                         double feature_rt, DoubleList& feature_intensities,
                         const std::set<String>& transition_ids = std::set<String>());

      /// Score a feature
      void scoreFeature_(Feature& feature);

  public:

      void initialize(const TargetedExperiment& library, const Size n_decoys, const Size n_transitions, const TransformationDescription& rt_trafo)
      {
        library_ = library; 
        n_decoys_ = n_decoys;
        n_transitions_ = n_transitions;
        rt_trafo_ = rt_trafo;
      }

      void initializeGlm(double intercept, double rt_coef, double int_coef)
      {
        glm_.intercept = intercept;
        glm_.rt_coef = rt_coef;
        glm_.int_coef = int_coef;
      }

      /**
        @brief Score a feature map -> make sure the class is properly initialized

        both functions initializeGlm and initialize need to be called first.

        The input to the program is 
        - a transition library which contains peptides with corresponding assays.
        - a feature map where each feature corresponds to an assay (mapped with
          MetaValue "PeptideRef") and each feature has as many subordinates as the
          assay has transitions (mapped with MetaValue "native_id").

      */
      void scoreMap(FeatureMap & features)
      {
        // are there enough assays in the library?
        Size n_assays = library_.getPeptides().size();
        if (n_assays < 2)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "There need to be at least 2 assays in the library for ConfidenceScoring.");

        }
        if (n_assays - 1 < n_decoys_)
        {
          OPENMS_LOG_WARN << "Warning: Parameter 'decoys' (" << n_decoys_ 
                   << ") is higher than the number of unrelated assays in the "
                   << "library (" << n_assays - 1 << "). "
                   << "Using all unrelated assays as decoys." << std::endl;
        }
        if (n_assays - 1 <= n_decoys_) n_decoys_ = 0; // use all available assays

        decoy_index_.resize(n_assays);
        for (Size i = 0; i < n_assays; ++i) decoy_index_[i] = boost::numeric_cast<Int>(i);

        // build mapping between assays and transitions:
        OPENMS_LOG_DEBUG << "Building transition map..." << std::endl;
        for (Size i = 0; i < library_.getTransitions().size(); ++i)
        {
          const String& ref = library_.getTransitions()[i].getPeptideRef();
          transition_map_[ref].push_back(boost::numeric_cast<Int>(i));
        }
        // find min./max. RT in the library:
        OPENMS_LOG_DEBUG << "Determining retention time range..." << std::endl;
        rt_norm_.min_rt = std::numeric_limits<double>::infinity();
        rt_norm_.max_rt = -std::numeric_limits<double>::infinity();
        for (std::vector<TargetedExperiment::Peptide>::const_iterator it = 
               library_.getPeptides().begin(); it != library_.getPeptides().end();
             ++it)
        {
          double current_rt = getAssayRT_(*it);
          if (current_rt == -1.0) continue; // indicates a missing value
          rt_norm_.min_rt = std::min(rt_norm_.min_rt, current_rt);
          rt_norm_.max_rt = std::max(rt_norm_.max_rt, current_rt);
        }

        // log scoring progress:
        OPENMS_LOG_DEBUG << "Scoring features..." << std::endl;
        startProgress(0, features.size(), "scoring features");

        for (FeatureMap::Iterator feat_it = features.begin(); 
             feat_it != features.end(); ++feat_it)
        {
          OPENMS_LOG_DEBUG << "Feature " << feat_it - features.begin() + 1 
                    << " (ID '" << feat_it->getUniqueId() << "')"<< std::endl;
          scoreFeature_(*feat_it);
          setProgress(feat_it - features.begin());
        }
        endProgress();

      }

  };

}

