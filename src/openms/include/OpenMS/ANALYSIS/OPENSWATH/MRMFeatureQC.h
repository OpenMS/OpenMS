// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/LogStream.h>

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
    MRMFeatureQC();

    /// Destructor
    ~MRMFeatureQC();
    //@}

    // Members
    //
    /**@brief Quality Controls (QCs) for individual components

      A component is analogous to a transition or subfeature.
      A component is a transition that corresponds to a single peptide or metabolite

    */
    struct ComponentQCs
    {
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


