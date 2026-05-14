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

  /**
    @brief Confidence scoring for SRM/MRM/PRM features against a targeted assay library.

    For each detected feature the algorithm scores the feature against
      - the matching @em true assay (referenced by the feature's @c PeptideRef MetaValue), and
      - a random sample of @em decoy assays (randomly chosen non-matching assays from the
        same library) used to estimate the score distribution under the null hypothesis.

    The per-assay score is the output of a binomial GLM (see @c GLM_) on two features
    derived from the feature/assay pair:
      - the difference between the feature's RT and the assay's library RT, after
        normalisation to the 0–100 range (see @c RTNorm_); the GLM itself applies the
        square (its rt_coef multiplies @c diff_rt*diff_rt).
      - the Manhattan distance between the feature's transition-intensity vector and
        the assay's expected transition-intensity vector (intensities matched by Q3 m/z).

    Initialisation flow:
      -# initialize() — install the library, decoy count, transition count, and RT transform.
      -# initializeGlm() — install the pre-fitted GLM coefficients (intercept / RT² / intensity).
      -# scoreMap() — iterate all features in the input @ref FeatureMap and write the per-feature
         scoring result back (see scoreMap() for the exact MetaValue set + overall-quality update).

    Feature → assay matching is by MetaValue: each feature must carry a @c "PeptideRef" string
    (the assay id), and each subordinate (one per transition) must carry a @c "native_id" string
    matching a transition id from the library.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI ConfidenceScoring :
      public ProgressLogger
  {
  public:

      /**
        @brief Construct an empty scorer.

        @param[in] test_mode_ If @c true, seed the internal random shuffler with @c time(NULL)
                              so that decoy selection is deterministic only within a single
                              process run. If @c false (default), seed with 0 for fully
                              reproducible decoy draws across runs and platforms.
      */
      explicit ConfidenceScoring(bool test_mode_ = false);

      /// Destructor
      ~ConfidenceScoring() override {}

  protected:

      /**
        @brief Binomial GLM used to map (squared-normalised-RT-diff, intensity-distance) -> [0, 1] confidence.

        Implements @f$\mathrm{sigmoid}(\mathrm{intercept} + \mathrm{rt\_coef} \cdot \Delta\mathrm{RT}^2 + \mathrm{int\_coef} \cdot d_\mathrm{int})@f$
        with coefficients fitted externally (see @ref initializeGlm).
      */
      struct GLM_
      {
        double intercept; ///< GLM intercept term
        double rt_coef;   ///< GLM coefficient on the squared RT difference (units: 1/RT²)
        double int_coef;  ///< GLM coefficient on the Manhattan intensity distance

        /// Evaluate the GLM at (@p diff_rt, @p dist_int); returns a probability in [0, 1].
        double operator()(double diff_rt, double dist_int) const
        {
          double lm = intercept + rt_coef * diff_rt * diff_rt +
            int_coef * dist_int;
          return 1.0 / (1.0 + exp(-lm));
        }
      } glm_;

      /**
        @brief Map RT values into the [0, 100] interval using min/max RT of the assay library.

        Populated by scoreMap() from the library RTs and applied before computing the GLM's
        @c diff_rt term so the GLM coefficients are unit-free across different assay-library
        RT ranges.
      */
      struct RTNorm_
      {
        double min_rt; ///< Smallest assay RT in the library; set by scoreMap()
        double max_rt; ///< Largest assay RT in the library; set by scoreMap()

        /// Map @p rt into [0, 100] using the cached min/max library RTs.
        double operator()(double rt) const
        {
          return (rt - min_rt) / (max_rt - min_rt) * 100;
        }
      } rt_norm_;

      TargetedExperiment library_; ///< Targeted-assay library: one peptide per assay, each with its transitions

      IntList decoy_index_; ///< Indexes into @c library_.getPeptides() used as decoys for the current feature

      Size n_decoys_; ///< Number of decoy assays to sample per feature (0 = use all unrelated assays as decoys)

      std::map<String, IntList> transition_map_; ///< Lookup assay-id -> indexes into @c library_.getTransitions()

      Size n_transitions_; ///< Number of top-intensity transitions to keep when computing the intensity-distance term (0 = keep all)

      TransformationDescription rt_trafo_; ///< Optional RT transformation applied to measured feature RTs before comparison with library RTs

      Math::RandomShuffler shuffler_; ///< Random shuffler used to draw decoy samples (seed depends on test mode — see ctor)

      /// Permute @c decoy_index_ in place to pick a fresh random decoy sample for the next feature
      void chooseDecoys_();

      /// Manhattan (L1) distance between two equal-length vectors
      double manhattanDist_(DoubleList x, DoubleList y);

      /// Read the (single) retention time from an assay's @ref TargetedExperiment::Peptide; the assay is required to have an RT set (checked in debug builds)
      double getAssayRT_(const TargetedExperiment::Peptide& assay);

      /**
        @brief Score one feature against one candidate assay.

        Reduces the assay's transitions to those listed in @p transition_ids (or all if the
        set is empty), matches them up with the feature's transition intensities by Q3 m/z,
        computes the Manhattan distance on the resulting paired intensity vectors, evaluates
        the GLM on (RT² diff, intensity distance) and returns the GLM output in [0, 1].

        @param[in]  assay                Candidate assay (true or decoy).
        @param[in]  feature_rt           Feature retention time after @c rt_trafo_ and @c rt_norm_.
        @param[in,out] feature_intensities Feature's transition intensities; reordered to match the assay's transition list.
        @param[in]  transition_ids       Optional subset of transition ids to consider; empty means "use all transitions of the assay".
        @return GLM-derived score in [0, 1].
      */
      double scoreAssay_(const TargetedExperiment::Peptide& assay,
                         double feature_rt, DoubleList& feature_intensities,
                         const std::set<String>& transition_ids = std::set<String>());

      /// Score one feature against its matching assay plus a random decoy sample; writes the per-assay scores into @p feature's MetaValues
      void scoreFeature_(Feature& feature);

  public:

      /**
        @brief Install the configuration needed before scoreMap() can run.

        @param[in] library        Targeted-assay library; must have at least 2 peptides for scoring to be defined.
        @param[in] n_decoys       Number of decoy assays to sample per feature (clamped to all-available if larger, with warning).
        @param[in] n_transitions  Number of top-intensity transitions to keep when computing the intensity-distance term; 0 keeps all.
        @param[in] rt_trafo       RT transformation applied to feature RTs before comparison against library RTs.
      */
      void initialize(const TargetedExperiment& library, const Size n_decoys, const Size n_transitions, const TransformationDescription& rt_trafo)
      {
        library_ = library;
        n_decoys_ = n_decoys;
        n_transitions_ = n_transitions;
        rt_trafo_ = rt_trafo;
      }

      /**
        @brief Install the GLM coefficients fitted externally on a training set.

        Coefficients correspond to the @ref GLM_ formula
        @f$\mathrm{sigmoid}(\mathrm{intercept} + \mathrm{rt\_coef} \cdot \Delta\mathrm{RT}^2 + \mathrm{int\_coef} \cdot d_\mathrm{int})@f$.

        @param[in] intercept GLM intercept term.
        @param[in] rt_coef   Coefficient on @f$\Delta\mathrm{RT}^2@f$ (using normalised RT in [0, 100]).
        @param[in] int_coef  Coefficient on the Manhattan intensity distance.
      */
      void initializeGlm(double intercept, double rt_coef, double int_coef)
      {
        glm_.intercept = intercept;
        glm_.rt_coef = rt_coef;
        glm_.int_coef = int_coef;
      }

      /**
        @brief Score every feature in @p features in place, writing per-feature scores and an updated overall-quality value.

        Prerequisite: both @ref initialize and @ref initializeGlm must have been called first.

        The input contract is:
          - @p features each carry the assay id in MetaValue @c "PeptideRef".
          - Each feature has one @ref Feature::getSubordinates entry per assay transition,
            each carrying the transition id in MetaValue @c "native_id".
          - The targeted-assay library @ref initialize() was called with must contain at
            least 2 peptides; the function throws @ref Exception::IllegalArgument otherwise.

        Per-feature side effects:
          - MetaValue @c "GLM_score" is set to the GLM score against the true assay.
          - MetaValue @c "local_FDR" is set to the rank-based local FDR estimated against
            the per-feature decoy scores.
          - The feature's overall quality is set to @c 1.0 - local_FDR.

        Other side effects:
          - The library's RT range is computed here and cached in @c rt_norm_.
          - If @c n_decoys_ exceeds the number of unrelated assays in the library, it is
            clamped to 0 (= use all available assays as decoys) and a warning is logged.

        @param[in,out] features Feature map to score; per-feature MetaValues and overall quality are written in place.
        @throws Exception::IllegalArgument If @c library_ has fewer than 2 peptides.
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
