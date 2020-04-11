// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  /**

    @brief The MRMFeatureFilter either flags components and/or transitions that do not pass the QC criteria or filters out
      components and/or transitions that do not pass the QC criteria.

    @htmlinclude OpenMS_MRMFeatureFilter.parameters

  */
  class OPENMS_DLLAPI MRMFeatureFilter :
    public DefaultParamHandler
  {

public:

    //@{
    /// Constructor
    MRMFeatureFilter();

    /// Destructor
    ~MRMFeatureFilter() override;
    //@}

    /**
      @brief Get the class' default parameters

      @param[out] params Output parameters
    */
    void getDefaultParameters(Param& params) const;

    /// Synchronize members with param class
    void updateMembers_() override;

    /**
      @brief Flags or filters features and subordinates in a FeatureMap

      @param features FeatureMap to flag or filter
      @param filter_criteria MRMFeatureQC class defining QC parameters
      @param transitions transitions from a TargetedExperiment

    */
    void FilterFeatureMap(FeatureMap& features, const MRMFeatureQC& filter_criteria,
      const TargetedExperiment & transitions);

    /**
      @brief Estimate the lower and upper bound values for the MRMFeatureQC class based on a
        user supplied template.

      @param[in] samples Samples (typically Standards) from which to estimate the lower and upper bound values for the MRMFeatureQC members
      @param[in, out] filter_template A MRMFeatureQC class that will be used as a template to fill in the estimated lower and upper values.
        A "template" is needed so that the MRMFeatureQC::meta_value_qc parameters that the FeatureMap::MetaValues that user would like estimated are known.
      @param transitions transitions from a TargetedExperiment
    */
    void EstimateDefaultMRMFeatureQCValues(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions);

    /**
      @brief Transfer the lower and upper bound values for the calculated concentrations
        based off of the AbsoluteQuantitationMethod

      @param[in] quantitation_methods The absolute quantitation methods that has been determined for each component
      @param[in, out] filter_template A MRMFeatureQC class that will be used as a template to fill in the 
        MRMFeatureQC::ComponentQCs.calculated_concetration bounds based on the LLOQ and ULOQ values given in the quantitation_method.
    */
    void TransferLLOQAndULOQToCalculatedConcentrationBounds(const std::vector<AbsoluteQuantitationMethod>& quantitation_method, MRMFeatureQC& filter_template);

    /**
      @brief Estimate the feature variability from multiple pooled QC samples
        or replicate Unknown samples.  The returned map can then be used by
        `FilterFeatureMap` in order to filter on either the 
        `MRMFeatureQC::ComponentQCs.perc_rsd_qc` and `MRMFeatureQC::ComponentGroupQCs.perc_rsd_qc`
        or `MRMFeatureQC::ComponentQCs.perc_rsd_rep` members based on the
        position of the returned map in the arguments list to `FilterFeatureMap`

      @param[in] qc_or_reps multiple pooled QC samples or replicate Unknown samples FeatureMaps
      @param[in, out] perc_rsd A consensus FeatureMap of %RSD values.
      @param[in, out] filter_template A MRMFeatureQC class that will be used as a template to determine what FeatureMap values
        to estimate the %RSD for
      @param transitions transitions from a TargetedExperiment
    */
    void EstimatePercRSD(const std::vector<FeatureMap>& qc_or_reps, FeatureMap& perc_rsd, MRMFeatureQC& filter_template, const TargetedExperiment& transitions);

    /**
      @brief Estimate the background interference level from Blank samples.
        The returned map can then be used by `FilterFeatureMap` in order to filter on the `MRMFeatureQC::ComponentQCs.perc_background`
        or the `MRMFeatureQC::ComponentGroupQCs.perc_background` members based on the
        position of the returned map in the arguments list to `FilterFeatureMap`

      @param[in] qc_or_reps multiple pooled QC samples or replicate Unknown samples FeatureMaps
      @param[in, out] background A consensus feature map of background intensity values.
      @param[in, out] filter_template A MRMFeatureQC class that will be used as a template to determine what FeatureMap values
        to estimate the %RSD for
      @param transitions transitions from a TargetedExperiment
    */
    void EstimateBackgroundInterferences(const std::vector<FeatureMap>& blanks, FeatureMap& background, MRMFeatureQC& filter_template, const TargetedExperiment& transitions);

    /**
      @brief Calculates the ion ratio between two transitions

      @param component_1 component of the numerator
      @param component_2 component of the denominator
      @param feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.
    */
    double calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name);

    /**
      @brief Calculates the retention time difference between two features

      @param component_1 First eluting component
      @param component_2 Second eluting component

      @return The difference.
    */
    double calculateRTDifference(Feature & component_1, Feature & component_2);

    /**
      @brief Calculates the resolution between two features

      @param component_1 component 1
      @param component_2 component 2

      @return The difference.
    */
    double calculateResolution(Feature & component_1, Feature & component_2);

    /**
      @brief Checks if the metaValue is within the user specified range

      @param[in] component component of the numerator
      @param[in] meta_value_key Name of the metaValue
      @param[in] meta_value_l Lower bound (inclusive) for the metaValue range
      @param[in] meta_value_u Upper bound (inclusive) for the metaValue range
      @param[out] key_exists true if the given key is found, false otherwise

      @return True if the metaValue is within the bounds, and False otherwise.
    */
    bool checkMetaValue(
      const Feature & component,
      const String & meta_value_key,
      const double & meta_value_l,
      const double & meta_value_u,
      bool & key_exists
    ) const;

    /**
      @brief Checks if the metaValue is within the user specified range

      @param[in] component component of the numerator
      @param[in] meta_value_key Name of the metaValue
      @param[in, out] meta_value_l Lower bound (inclusive) for the metaValue range
      @param[in, out] meta_value_u Upper bound (inclusive) for the metaValue range
      @param[out] key_exists true if the given key is found, false otherwise
    */
    void updateMetaValue(
      const Feature & component,
      const String & meta_value_key,
      double & meta_value_l,
      double & meta_value_u,
      bool & key_exists
    ) const;

    /**
      @brief Count the number of heavy/light labels and quantifying/detecting/identifying transitions

      @param component component_group with subordinates
      @param transitions transitions from a TargetedExperiment

      @return Map of labels/transition types and their corresponding number.
    */
    std::map<String,int> countLabelsAndTransitionTypes(const Feature & component_group,
      const TargetedExperiment & transitions);

    /**
      @brief Sorts the messages and returns a copy without duplicates

      @param[in] messages A StringList containing the failure messages

      @return A copy of the input, without duplicates
    */
    StringList getUniqueSorted(const StringList& messages) const;

    template <typename T>
    bool checkRange(const T& value, const T& value_l, const T& value_u) const;
    template <typename T>
    void updateRange(const T& value, T& value_l, T& value_u) const;

private:
    // Members
    /// flag or filter (i.e., remove) features that do not pass the QC
    String flag_or_filter_;
    /// whether to use intensity or calculated concentration for background interference estimation and filtering
    bool use_calculated_concentration_background_;
    /// whether to use intensity or calculated concentration for %RSD estimation and filtering
    bool use_calculated_concentration_rsd_;
  };
}

