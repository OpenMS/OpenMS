// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  class AbsoluteQuantitationMethod;

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
      const TargetedExperiment& transitions);

    /**
      @brief Flags or filters features and subordinates in a FeatureMap based on a user
        defined set of filter values derived from calling `EstimatePercRSD`.  The user supplied
        filter_criteria represents the bounds on acceptable PercentRSD values.
        NOTE that flagging nor filtering will be done on the labels and transitions type counts.

      @param features FeatureMap to flag or filter
      @param filter_criteria MRMFeatureQC class defining QC parameters defining the acceptable limits of the PercentRSD
        where PercentRSD = (value std dev)/(value mean)*100Percent
      @param filter_values MRMFeatureQC class filled with bounds representing the PercentRSD found in e.g., pooled QC samples or replicate Unknown samples

    */
    void FilterFeatureMapPercRSD(FeatureMap& features, const MRMFeatureQC& filter_criteria, const MRMFeatureQC& filter_values);

    /**
      @brief Flags or filters features and subordinates in a FeatureMap based on a user
        defined set of filter values derived from calling `EstimateBackgroundInterferences`.  The user supplied
        filter_criteria represents the bounds on acceptable PercentBackgroundInterference values.
        NOTE that filtering is only done on the `Intensity` member.


      @param features FeatureMap to flag or filter
      @param filter_criteria MRMFeatureQC class defining QC parameters defining the acceptable limits of PercentBackgroundInterference
        where PercentBackgroundInterference = (value Sample)/(value Blank)*100Percent
      @param filter_values MRMFeatureQC class filled with bounds representing the average values found in e.g., pooled QC samples or replicate Unknown samples

    */
    void FilterFeatureMapBackgroundInterference(FeatureMap& features, const MRMFeatureQC& filter_criteria, const MRMFeatureQC& filter_values);

    /**
      @brief Estimate the lower and upper bound values for the MRMFeatureQC class based on a
        user supplied template.  The template can either be initialized from the first sample 
        (meaning all initial template values are over written) or not (meaning the initial template
        values will be updated if the ranges are found to be too narrow)

      @param[in] samples Samples (typically Standards) from which to estimate the lower and upper bound values for the MRMFeatureQC members
      @param[in,out] filter_template A MRMFeatureQC class that will be used as a template to fill in the estimated lower and upper values.
        A "template" is needed so that the MRMFeatureQC::meta_value_qc parameters that the FeatureMap::MetaValues that user would like estimated are known.
      @param[in] transitions transitions from a TargetedExperiment
      @param[in] init_template_values Boolean indicating whether to initialize the template values based on the first sample
    */
    void EstimateDefaultMRMFeatureQCValues(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions, const bool& init_template_values) const;

    /**
      @brief Transfer the lower and upper bound values for the calculated concentrations
        based off of the AbsoluteQuantitationMethod

      @param[in] quantitation_method The absolute quantitation methods that has been determined for each component
      @param[in,out] filter_template A MRMFeatureQC class that will be used as a template to fill in the 
        MRMFeatureQC::ComponentQCs.calculated_concentration bounds based on the LLOQ and ULOQ values given in the quantitation_method.
    */
    void TransferLLOQAndULOQToCalculatedConcentrationBounds(const std::vector<AbsoluteQuantitationMethod>& quantitation_method, MRMFeatureQC& filter_template);

    /**
      @brief Estimate the feature variability as measured by PercentRSD from multiple pooled QC samples
        or replicate Unknown samples.  The returned filter_template can then be used by
        `FilterFeatureMapPercRSD` in order to filter based on the PercentRSD user defined limits.

      @param[in] samples multiple pooled QC samples or replicate Unknown samples FeatureMaps
      @param[in,out] filter_template A MRMFeatureQC class that will be used as a template to determine what FeatureMap values
        to estimate the PercentRSD for.  The PercentRSD values will be stored in the upper bound parameter of the filter_template
      @param[in] transitions transitions from a TargetedExperiment
    */
    void EstimatePercRSD(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const;

    /**
      @brief Estimate the background interference level based on the average values from Blank samples.
        The returned filter_template can then be used by `FilterFeatureMapBackgroundInterference` 
        in order to filter on the `Intensity` members of MRMFeatureQC::ComponentGroupQCs and MRMFeatureQC::ComponentQCs.

      @param[in] samples multiple Blank samples to estimate the background intensity values FeatureMaps
      @param[in,out] filter_template A MRMFeatureQC class that will be used as a template to determine what FeatureMap values
        to estimate the PercentInterference.  The average values will be stored in the upper bound parameter of the filter_template
      @param[in] transitions transitions from a TargetedExperiment
    */
    void EstimateBackgroundInterferences(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const;

    /**
      @brief Calculates the ion ratio between two transitions

      @param component_1 component of the numerator
      @param component_2 component of the denominator
      @param feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.
    */
    double calculateIonRatio(const Feature& component_1, const Feature& component_2, const String& feature_name) const;

    /**
      @brief Calculates the retention time difference between two features

      @param component_1 First eluting component
      @param component_2 Second eluting component

      @return The difference.
    */
    double calculateRTDifference(Feature& component_1, Feature& component_2) const;

    /**
      @brief Calculates the resolution between two features

      @param component_1 component 1
      @param component_2 component 2

      @return The difference.
    */
    double calculateResolution(Feature& component_1, Feature& component_2) const;

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
      const Feature& component,
      const String& meta_value_key,
      const double& meta_value_l,
      const double& meta_value_u,
      bool& key_exists
    ) const;

    /**
      @brief Updates the metaValue ranges based on the value given

      @param[in] component component of the numerator
      @param[in] meta_value_key Name of the metaValue
      @param[in,out] meta_value_l Lower bound (inclusive) for the metaValue range
      @param[in,out] meta_value_u Upper bound (inclusive) for the metaValue range
      @param[out] key_exists true if the given key is found, false otherwise
    */
    void updateMetaValue(
      const Feature& component,
      const String& meta_value_key,
      double& meta_value_l,
      double& meta_value_u,
      bool& key_exists
    ) const;

    /**
      @brief Uses the supplied value to set the metaValue ranges

      @param[in] component component of the numerator
      @param[in] meta_value_key Name of the metaValue
      @param[in,out] meta_value_l Lower bound (inclusive) for the metaValue range
      @param[in,out] meta_value_u Upper bound (inclusive) for the metaValue range
      @param[out] key_exists true if the given key is found, false otherwise
    */
    void setMetaValue(
      const Feature& component,
      const String& meta_value_key,
      double& meta_value_l,
      double& meta_value_u,
      bool& key_exists
    ) const;

    /**
      @brief Uses the supplied value to initialize the metaValue ranges to the same value

      @param[in] component component of the numerator
      @param[in] meta_value_key Name of the metaValue
      @param[in,out] meta_value_l Lower bound (inclusive) for the metaValue range
      @param[in,out] meta_value_u Upper bound (inclusive) for the metaValue range
      @param[out] key_exists true if the given key is found, false otherwise
    */
    void initMetaValue(
      const Feature& component,
      const String& meta_value_key,
      double& meta_value_l,
      double& meta_value_u,
      bool& key_exists
    ) const;

    /**
      @brief Count the number of heavy/light labels and quantifying/detecting/identifying transitions

      @param component_group Component group with subordinates
      @param transitions Transitions from a TargetedExperiment

      @return Map of labels/transition types and their corresponding number.
    */
    std::map<String,int> countLabelsAndTransitionTypes(const Feature& component_group,
      const TargetedExperiment& transitions) const;

    /**
      @brief Sorts the messages and returns a copy without duplicates

      @param[in] messages A StringList containing the failure messages

      @return A copy of the input, without duplicates
    */
    StringList getUniqueSorted(const StringList& messages) const;

    /**
      @brief Accumulate feature values from a list of FeatureMaps

      @param[out] filter_values A list of MRMFeatureQC objects filled with the values determined by the filter_template from the samples
      @param[in] samples A list of feature maps
      @param[in] filter_template A MRMFeatureQC object that will be used as a template to fill in values derived from the sample feature maps
      @param[in] transitions transitions from a TargetedExperiment
    */
    void accumulateFilterValues(std::vector<MRMFeatureQC>& filter_values, const std::vector<FeatureMap>& samples, const MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const;

    /**
      @brief Set all members in MRMFeatureQC to zero 

      @param[out] filter_zeros A MRMFeatureQC object whose members have been set to 0
      @param[in] filter_template A MRMFeatureQC object that will be used as a template to fill in values
    */
    void zeroFilterValues(MRMFeatureQC& filter_zeros, const MRMFeatureQC& filter_template) const;

    /**
      @brief Calculate the mean of each MRMFeatureQC parameter from a list of MRMFeatureQC classes

      @param[out] filter_mean A MRMFeatureQC object whose members will be replaced by the mean values
      @param[in] filter_values A list of MRMFeatureQC objects
      @param[in] filter_template A MRMFeatureQC object that will be used as a template to fill in values
    */
    void calculateFilterValuesMean(MRMFeatureQC& filter_mean, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_template) const;

    /**
      @brief Calculate the var of each MRMFeatureQC parameter from a list of MRMFeatureQC classes

      @param[out] filter_var A MRMFeatureQC object whose members will be replaced by the variance of the values
      @param[in] filter_values A list of MRMFeatureQC objects
      @param[in] filter_mean A MRMFeatureQC object with the mean values of the filter_values
      @param[in] filter_template A MRMFeatureQC object that will be used as a template to fill in values
    */
    void calculateFilterValuesVar(MRMFeatureQC& filter_var, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_mean, const MRMFeatureQC& filter_template) const;

    /**
      @brief Calculate the relative standard deviation (PercentRSD) of each MRMFeatureQC parameter from pre-computed mean and variance values

      @param[out] filter_rsd A MRMFeatureQC object whose members will be replaced by PercentRSD
      @param[in] filter_var A MRMFeatureQC object with the variance
      @param[in] filter_mean A MRMFeatureQC object with the mean
    */
    void calculateFilterValuesPercRSD(MRMFeatureQC& filter_rsd, const MRMFeatureQC& filter_mean, const MRMFeatureQC& filter_var) const;

    /// Checks that the range of value is bracketed by value_l and value_u
    template <typename T>
    bool checkRange(const T& value, const T& value_l, const T& value_u) const;

    /// Updates value_l and value_u according to whether value is greater than value_u or less than value_l
    template <typename T>
    void updateRange(const T& value, T& value_l, T& value_u) const;

    /// Sets value_l and value_u to bracket the range 0 to value or value to 0 depending on if value is >0
    template <typename T>
    void setRange(const T& value, T& value_l, T& value_u) const;

    /// Sets value_l and value_u to value
    template <typename T>
    void initRange(const T& value, T& value_l, T& value_u) const;

private:
    // Members
    /// flag or filter (i.e., remove) features that do not pass the QC
    String flag_or_filter_;
  };
}

