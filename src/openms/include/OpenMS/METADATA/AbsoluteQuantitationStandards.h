// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**
    @brief Bridges sample-level concentration tables with the @ref FeatureMap
           output of a feature finder, for building calibration curves used in
           absolute quantitation.

    The user supplies a list of @ref runConcentration entries -- one per
    (sample, target component, internal-standard) triple with the known
    concentrations -- together with the @ref FeatureMap instances produced
    for each sample. The class then locates each component (and its
    internal standard, when given) in the matching @ref FeatureMap and
    pairs it with the entry's concentration, yielding @ref featureConcentration
    records that downstream calibration code can consume.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI AbsoluteQuantitationStandards
  {

public:
    /// Default constructor.
    AbsoluteQuantitationStandards() = default;

    /// Destructor.
    ~AbsoluteQuantitationStandards() = default;

    /**
      @brief One (sample, target, internal-standard, concentration) entry as
             provided by the user.

      Typically populated from a calibration-table file before the standards
      are mapped to feature data.
    */
    struct runConcentration
    {
      String sample_name;              ///< Sample identifier; matched against @ref FeatureMap::getPrimaryMSRunPath after stripping a @c .mzML or @c .txt suffix.
      String component_name;           ///< Identifier of the target component (matched against subordinate-feature meta value @c "native_id").
      String IS_component_name;        ///< Identifier of the paired internal standard, or empty if there is none.
      double actual_concentration;     ///< Known concentration of the target component in this sample.
      double IS_actual_concentration;  ///< Known concentration of the internal standard in this sample.
      String concentration_units;      ///< Units the concentration values are expressed in.
      double dilution_factor;          ///< Sample dilution factor.
    };

    /**
      @brief One target component matched to a feature in the corresponding
             @ref FeatureMap, paired with its known concentration.

      Produced by @ref mapComponentsToConcentrations and
      @ref getComponentFeatureConcentrations.
    */
    struct featureConcentration
    {
      Feature feature;                 ///< Subordinate feature found in the matching @ref FeatureMap (see @ref mapComponentsToConcentrations for the match rule).
      Feature IS_feature;              ///< Subordinate feature for the internal standard. Default-constructed when @ref runConcentration::IS_component_name was empty or no IS subordinate was found.
      double actual_concentration;     ///< Copied from @ref runConcentration::actual_concentration.
      double IS_actual_concentration;  ///< Copied from @ref runConcentration::IS_actual_concentration.
      String concentration_units;      ///< Copied from @ref runConcentration::concentration_units.
      double dilution_factor;          ///< Copied from @ref runConcentration::dilution_factor.
    };

    /**
      @brief Match each @ref runConcentration to a feature in @p feature_maps and
             produce one @ref featureConcentration record per match, keyed by
             component name.

      For every entry in @p run_concentrations:
        - Entries with an empty @c sample_name or @c component_name are
          skipped.
        - A @ref FeatureMap from @p feature_maps is considered a candidate
          when its first @ref FeatureMap::getPrimaryMSRunPath entry, with a
          trailing @c .mzML or @c .txt suffix removed, equals the entry's
          @c sample_name. When the primary-run-path list is empty the
          candidate is accepted unconditionally.
        - Within the chosen @ref FeatureMap, the target component is located
          as a subordinate feature whose @c "native_id" meta value equals
          @c component_name. If no such subordinate exists, the entry is
          skipped.
        - When @c IS_component_name is non-empty, the same subordinate
          lookup is performed for the internal standard, but the entry is
          kept even when the IS lookup fails (the resulting
          @ref featureConcentration::IS_feature is then default-constructed).
        - At most one @ref FeatureMap is matched per entry; matching stops
          at the first hit.

      @param[in]  run_concentrations           Calibration entries to consider.
      @param[in]  feature_maps                 Feature maps to search.
      @param[out] components_to_concentrations Result, keyed by
                                               @c component_name. Cleared at
                                               the start of the call.
    */
    void mapComponentsToConcentrations(
      const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
      const std::vector<FeatureMap>& feature_maps,
      std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>>& components_to_concentrations
    ) const;

    /**
      @brief Return the @ref featureConcentration records for a single
             component name.

      Equivalent to filtering @p run_concentrations down to entries whose
      @c component_name matches @p component_name and running
      @ref mapComponentsToConcentrations on the result; the records for
      @p component_name are then written to @p feature_concentrations.
      If no matching record is produced, @p feature_concentrations is left
      unchanged.

      @param[in]  run_concentrations    Calibration entries to filter.
      @param[in]  feature_maps          Feature maps to search.
      @param[in]  component_name        Only entries with this component name are considered.
      @param[out] feature_concentrations Result vector; only assigned when at least one match is produced.
    */
    void getComponentFeatureConcentrations(
      const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
      const std::vector<FeatureMap>& feature_maps,
      const String& component_name,
      std::vector<AbsoluteQuantitationStandards::featureConcentration>& feature_concentrations
    ) const;

private:
    /**
      @brief Locate a subordinate feature whose @c "native_id" meta value equals @p component_name.

      Walks every @ref Feature in @p feature_map and, for each, its
      subordinates; the first subordinate whose @c "native_id" meta value
      equals @p component_name is copied into @p feature_found.

      @param[in]  feature_map    Features to search (top-level features and their subordinates).
      @param[in]  component_name Required value of the subordinate's @c "native_id" meta value.
      @param[out] feature_found  The matching subordinate feature; only written on a successful return.
      @return @c true if a matching subordinate was found, @c false otherwise.
    */
    bool findComponentFeature_(
      const FeatureMap& feature_map,
      const String& component_name,
      Feature& feature_found
    ) const;
  };
}

