// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  /**

    @brief The MRMFeatureQC is a class to handle the parameters and options for
     MRMFeatureFilter.

     The format is based loosely on the TraML format and can be stored and loaded to disk using MRMFeatureQCFile.

     Quality control parameters are available on multiple levels:
       - the level of a single component (or transition) representing a single transition
       - the level of a component group (or transition group) representing a single chemical entity
       - the level of component group pairs (e.g. isotopic pairs) representing multiple chemical entities that may be related by isotopic pairing
  */
  class OPENMS_DLLAPI MRMFeatureQC
  {

public:

    //@{
    /// Constructor
    MRMFeatureQC() = default;

    /// Destructor
    ~MRMFeatureQC() = default;
    //@}

    // Members
    //
    /**@brief Quality Controls (QCs) for individual components

      A component is analogous to a transition or subfeature.
      A component is a transition that corresponds to a single peptide or metabolite

    */
    struct ComponentQCs
    {
      bool operator==(const ComponentQCs& other) const {
        bool members_eq = 
          std::tie(
            component_name,
            retention_time_l,
            retention_time_u,
            intensity_l,
            intensity_u,
            overall_quality_l,
            overall_quality_u
          ) == std::tie(
            other.component_name,
            other.retention_time_l,
            other.retention_time_u,
            other.intensity_l,
            other.intensity_u,
            other.overall_quality_l,
            other.overall_quality_u
          );
        auto compare_maps = [](std::pair<String, std::pair<double, double>> lhs, std::pair<String, std::pair<double, double>> rhs) {return (lhs.first == rhs.first && lhs.second.first == rhs.second.first && lhs.second.second == rhs.second.second); };
        bool meta_values_eq = std::equal(meta_value_qc.begin(), meta_value_qc.end(), other.meta_value_qc.begin(), compare_maps);
        return members_eq && meta_values_eq;
      }
      bool operator!=(const ComponentQCs& other) const
      {
        return !(*this == other);
      }

      /// name of the component
      String component_name;

      // Feature members
      /// retention time lower bound
      double retention_time_l { 0.0 };
      /// retention time upper bound
      double retention_time_u { 1e12 };
      /// intensity lower bound
      double intensity_l { 0.0 };
      /// intensity upper bound
      double intensity_u { 1e12 };
      /// overall quality lower bound
      double overall_quality_l { 0.0 };
      /// overall quality upper bound
      double overall_quality_u { 1e12 };

      /// Feature MetaValues
      std::map<String,std::pair<double,double>> meta_value_qc;

    };

    /**@brief Quality Controls (QCs) within a component group

      A component group is analogous to a transition group or feature.
      A component group includes all transitions that correspond to a given component (i.e., peptide or metabolite)

    */
    struct ComponentGroupQCs
    {
      bool operator==(const ComponentGroupQCs& other) const {
        bool members_eq =
          std::tie(
            component_group_name,
            retention_time_l,
            retention_time_u,
            intensity_l,
            intensity_u,
            overall_quality_l,
            overall_quality_u,
            n_heavy_l,
            n_heavy_u,
            n_light_l,
            n_light_u,
            n_detecting_l,
            n_detecting_u,
            n_quantifying_l,
            n_quantifying_u,
            n_identifying_l,
            n_identifying_u,
            n_transitions_l,
            n_transitions_u,
            ion_ratio_pair_name_1,
            ion_ratio_pair_name_2,
            ion_ratio_l,
            ion_ratio_u,
            ion_ratio_feature_name
          ) == std::tie(
            other.component_group_name,
            other.retention_time_l,
            other.retention_time_u,
            other.intensity_l,
            other.intensity_u,
            other.overall_quality_l,
            other.overall_quality_u,
            other.n_heavy_l,
            other.n_heavy_u,
            other.n_light_l,
            other.n_light_u,
            other.n_detecting_l,
            other.n_detecting_u,
            other.n_quantifying_l,
            other.n_quantifying_u,
            other.n_identifying_l,
            other.n_identifying_u,
            other.n_transitions_l,
            other.n_transitions_u,
            other.ion_ratio_pair_name_1,
            other.ion_ratio_pair_name_2,
            other.ion_ratio_l,
            other.ion_ratio_u,
            other.ion_ratio_feature_name
          );
        auto compare_maps = [](std::pair<String, std::pair<double, double>> lhs, std::pair<String, std::pair<double, double>> rhs) {return (lhs.first == rhs.first && lhs.second.first == rhs.second.first && lhs.second.second == rhs.second.second); };
        bool meta_values_eq = std::equal(meta_value_qc.begin(), meta_value_qc.end(), other.meta_value_qc.begin(), compare_maps);
        return members_eq && meta_values_eq;
      }
      bool operator!=(const ComponentGroupQCs& other) const
      {
        return !(*this == other);
      }
      /// name of the component group
      String component_group_name;

      /// retention time lower bound
      double retention_time_l { 0.0 };
      /// retention time upper bound
      double retention_time_u { 1e12 };
      /// intensity lower bound
      double intensity_l { 0.0 };
      /// intensity upper bound
      double intensity_u { 1e12 };
      /// overall quality lower bound
      double overall_quality_l { 0.0 };
      /// overall quality upper bound
      double overall_quality_u { 1e12 };

      // number of transitions and labels
      /// number of heavy ion lower bound
      Int n_heavy_l { 0 };
      /// number of heavy ion upper bound
      Int n_heavy_u { 100 };
      Int n_light_l { 0 };
      Int n_light_u { 100 };
      Int n_detecting_l { 0 };
      Int n_detecting_u { 100 };
      Int n_quantifying_l { 0 };
      Int n_quantifying_u { 100 };
      Int n_identifying_l { 0 };
      Int n_identifying_u { 100 };
      Int n_transitions_l { 0 };
      Int n_transitions_u { 100 };

      // Ion Ratio QCs
      String ion_ratio_pair_name_1;
      String ion_ratio_pair_name_2;
      double ion_ratio_l { 0.0 };
      double ion_ratio_u { 1e12 };
      String ion_ratio_feature_name;
      std::map<String,std::pair<double,double>> meta_value_qc;

    };

    /**@brief Quality Controls (QCs) for multiple components (between or within component_groups)

      This structure contains upper and lower bounds for parameters that involve two or more
      component groups.  For example, a quality control that is based on a minimum retention
      time difference between two components would be suitable for this struct.

    */
    struct ComponentGroupPairQCs
    {

      /// name of the component
      String component_group_name;
      /// name of the component to calculate the resolution or retention time
      String resolution_pair_name;
      /// resolution lower bound
      double resolution_l;
      /// resolution upper bound
      double resolution_u;
      /// retention time lower bound
      double rt_diff_l;
      /// retention time upper bound
      double rt_diff_u;
    };

    //members
    /// list of all component QCs
    std::vector<ComponentQCs> component_qcs;
    /// list of all component group QCs
    std::vector<ComponentGroupQCs> component_group_qcs;
    /// list of all component group pair QCs
    std::vector<ComponentGroupPairQCs> component_group_pair_qcs;
  };
}


