// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief Helper utilities for matching MS2 spectra to LC-MS features by precursor (m/z, RT).

    The class is a thin namespace holding two nested record types and one static algorithm:

      - @ref FeatureMappingInfo   – input bundle (feature maps + a @ref KDTreeFeatureMaps view)
      - @ref FeatureToMs2Indices  – output bundle (feature → spectrum indices map + unassigned list)
      - @ref assignMS2IndexToFeature – the actual assignment routine

    A typical caller builds a @ref FeatureMappingInfo once (loading the maps and constructing the
    kd-tree), then calls @ref assignMS2IndexToFeature per MS run. The kd-tree must be built from
    the same maps stored in @c feature_maps so that the indices returned by range queries refer
    back into the caller-visible feature pointers.

    Workflows that use this include SIRIUS export and other "annotate MS2 with quantified feature"
    pipelines.

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI FeatureMapping
  {
    public:

    /**
      @brief Input bundle: feature maps and a spatial index over them.

      @c feature_maps owns the raw @ref FeatureMap objects (one per input run); @c kd_tree is a
      @ref KDTreeFeatureMaps pre-populated with the same features and provides fast
      (RT, m/z) range queries. The caller is responsible for keeping the two consistent
      (i.e. building the kd-tree from the same maps stored here).
    */
    class FeatureMappingInfo
    {
    public:
      std::vector<FeatureMap> feature_maps; ///< feature data, one map per input run
      KDTreeFeatureMaps kd_tree;            ///< (RT, m/z) kd-tree referencing the features in @c feature_maps
    };

    /**
      @brief Output bundle: per-feature MS2 spectrum indices + the list of unassigned MS2 spectra.

      Keys of @c assignedMS2 are non-owning pointers into the @ref FeatureMap data held by the
      @ref FeatureMappingInfo passed to @ref assignMS2IndexToFeature; the values are spectrum
      indices into the @c spectra input. @c unassignedMS2 holds the indices of MS2 spectra that
      had at least one precursor but no feature inside the tolerance window (MS2 spectra without
      a precursor are silently dropped and appear in neither container).
    */
    class FeatureToMs2Indices
    {
    public:
       /// MS2 spectrum indices grouped by the feature they were assigned to
       std::map<const BaseFeature*, std::vector<size_t>> assignedMS2;
       /// Indices of MS2 spectra that had a precursor but no feature inside the tolerance window
       std::vector<size_t> unassignedMS2;
    };

    /**
      @brief Assign each MS2 spectrum to the feature whose precursor matches it best.

      For every MS2 spectrum in @p spectra:
        - Spectra at MS level != 2 are skipped.
        - Spectra with no precursor are silently dropped (not added to @c unassignedMS2).
        - Spectra with at least one precursor are matched against features in
          @c fm_info.kd_tree using a rectangular tolerance window:
            * RT:  @p precursor_rt_tolerance  (absolute, around the spectrum's RT)
            * m/z: @p precursor_mz_tolerance  (ppm if @p ppm is true, otherwise Th, around the first precursor's m/z)
          Features from any map are eligible (the kd-tree query is run with
          @c include_features_from_same_map = true).
        - If no feature falls in the window, the spectrum index is appended to
          @c unassignedMS2.
        - If several features fall in the window, the one with the smallest absolute m/z
          difference to the precursor wins (RT is not used as a tie-breaker).
        - Only the first precursor of the spectrum is considered (chimeric / multi-precursor
          spectra are matched against @c precursors()[0] only).

      The result groups spectrum indices by their winning feature; each spectrum index appears
      at most once in either the assigned or the unassigned collection.

      @param[in] spectra                  Run-level spectrum container (MS1 and MS2 mixed; MS2 spectra are processed).
      @param[in] fm_info                  Feature maps + kd-tree bundle (the kd-tree must reference @c fm_info.feature_maps).
      @param[in] precursor_mz_tolerance   Half-width of the m/z tolerance window (ppm if @p ppm, else Th).
      @param[in] precursor_rt_tolerance   Half-width of the RT tolerance window in seconds (absolute).
      @param[in] ppm                      If true, interpret @p precursor_mz_tolerance as ppm; otherwise as Th.
      @return Mapping result containing one entry per matched feature plus the list of unassigned MS2 spectra.
    */
    static FeatureToMs2Indices assignMS2IndexToFeature(const MSExperiment& spectra,
                                                       const FeatureMappingInfo& fm_info,
                                                       const double& precursor_mz_tolerance,
                                                       const double& precursor_rt_tolerance,
                                                       bool ppm);

  };
} // namespace OpenMS
