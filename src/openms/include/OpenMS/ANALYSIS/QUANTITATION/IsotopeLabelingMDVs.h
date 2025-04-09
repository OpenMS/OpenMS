// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ahmed Khalil $
// $Authors: Ahmed Khalil $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

namespace OpenMS
{
  /**
    @brief IsotopeLabelingMDVs is a class to support and analyze isotopic labeling experiments
            (i.e. MDVs : Mass Distribution Vectors, also known as Mass Isotopomer Distribution (MID))
  */
  class OPENMS_DLLAPI IsotopeLabelingMDVs :
    public DefaultParamHandler
  {
  public:
    //@{
    /// Constructor
    IsotopeLabelingMDVs();

    /// Destructor
    ~IsotopeLabelingMDVs() override;
    //@}
    
    enum class DerivatizationAgent
    {
      NOT_SELECTED,
      TBDMS,
      SIZE_OF_DERIVATIZATIONAGENT
    };
    
    enum class MassIntensityType
    {
      NORM_MAX,
      NORM_SUM,
      SIZE_OF_MASSINTENSITYTYPE
    };
    
    static const std::string NamesOfDerivatizationAgent[static_cast<int>(DerivatizationAgent::SIZE_OF_DERIVATIZATIONAGENT)];
    
    static const std::string NamesOfMassIntensityType[static_cast<int>(MassIntensityType::SIZE_OF_MASSINTENSITYTYPE)];
 
    /**
      @brief This function performs an isotopic correction to account for unlabeled abundances coming from
      the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
      MDV_corrected = CM^-1 * MDV_normalized_feature
      The formula is obtained from "The importance of accurately correcting for the natural abundance of stable isotopes",
      Midani et al, doi:10.1016/j.ab.2016.12.011
     
      @param  normalized_feature Feature with normalized values for each component and unlabeled chemical formula for each component group.
      @param[out] corrected_feature Feature with corrected values for each component.
      @param  correction_matrix  Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of
      naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent.
      @param  correction_matrix_agent name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied,
      only "TBDMS" is supported for now.
    */
    void isotopicCorrection(
      const Feature& normalized_feature,
      Feature& corrected_feature,
      const Matrix<double>& correction_matrix,
      const DerivatizationAgent& correction_matrix_agent);
    
    /**
      @brief This function performs an isotopic correction to account for unlabeled abundances coming from
      the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
      MDV_corrected = CM^-1 * MDV_normalized_feature
      The formula is obtained from "The importance of accurately correcting for the natural abundance of stable isotopes",
      Midani et al, doi:10.1016/j.ab.2016.12.011
     
      @param  measured_fm FeatureMap with normalized values for each component and unlabeled chemical formula for each component group.
      @param[out] corrected_fm FeatureMap with corrected values for each component.
      @param  correction_matrix Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of
      naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent.
      @param  correction_matrix_agent name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied,
      only "TBDMS" is supported for now.
    */
    void isotopicCorrections(
      const FeatureMap& measured_fm,
      FeatureMap& corrected_fm,
      const Matrix<double>& correction_matrix,
      const DerivatizationAgent& correction_matrix_agent);

    /**
      @brief This function calculates the isotopic purity of the MDV using the following formula:
      isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
      where n in M+n is represented as the index of the result.
      The formula is extracted from "High-resolution 13C metabolic flux analysis",
      Long et al, doi:10.1038/s41596-019-0204-0

      @param[in,out]  normalized_feature Feature with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a Feature with the calculated isotopic purity for the component group.
      @param[in]      experiment_data Vector of experiment data in percent.
      @param[in]      isotopic_purity_name Name of the isotopic purity tracer to be saved as a meta value.
    */
    void calculateIsotopicPurity(
      Feature& normalized_feature,
      const std::vector<double>& experiment_data, 
      const std::string& isotopic_purity_name);
    
    /**
      @brief This function calculates the isotopic purity of the MDVs using the following formula:
      isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
      where n in M+n is represented as the index of the result.
      The formula is extracted from "High-resolution 13C metabolic flux analysis",
      Long et al, doi:10.1038/s41596-019-0204-0

      @param[in,out]  measured_fm FeatureMap with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a FeatureMap with the calculated isotopic purity for the component group.
      @param[in]  experiment_data Vector of vectors of Experiment data in percent.
      @param[in]  isotopic_purity_name Vector of names of the isotopic purity tracer to be saved as a meta value.
    */
    void calculateIsotopicPurities(
      FeatureMap& measured_fm,
      const std::vector<std::vector<double>>& experiment_data,
      const std::vector<std::string>& isotopic_purity_name);

    /**
      @brief This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
      using average deviation to the mean. The result is mapped to the meta value "average_accuracy" in the updated Feature.
   
      @param[in,out]  normalized_feature Feature with normalized values for each component and the chemical formula of the component group. Out is a Feature with the component group accuracy and accuracy for the error for each component.
      @param[in]      feature_name Feature name to determine which features are needed to apply calculations on.
      @param[in]      fragment_isotopomer_theoretical_formula Empirical formula from which the theoretical values will be generated.
    */
    void calculateMDVAccuracy(
      Feature& normalized_feature,
      const std::string& feature_name, 
      const std::string& fragment_isotopomer_theoretical_formula);
    
    /**
       @brief This function calculates the accuracy of the MDVs as compared to the theoretical MDVs (only for 12C quality control experiments)
       using average deviation to the mean.
    
       @param[in]  normalized_fm FeatureMap with normalized values for each component and the chemical formula of the component group. Out is a FeatureMap with the component group accuracy and accuracy for the error for each component.
       @param[in]  feature_name Feature name to determine which features are needed to apply calculations on.
       @param[in]  fragment_isotopomer_theoretical_formulas A map of ProteinName/peptideRef to Empirical formula from which the theoretical values will be generated.
    */
    void calculateMDVAccuracies(
      FeatureMap& normalized_fm,
      const std::string& feature_name,
      const std::map<std::string, std::string>& fragment_isotopomer_theoretical_formulas);
 
    /**
      @brief This function calculates the mass distribution vector (MDV)
      either normalized to the highest mass intensity (norm_max) or normalized
      to the sum of all mass intensities (norm_sum)
     
      @param[in]   measured_feature Feature with measured intensity for each component.
      @param[out]  normalized_feature Feature with normalized values for each component.
      @param[in]   mass_intensity_type Mass intensity type (either norm_max or norm_sum).
      @param[in]   feature_name Feature name to determine which features are needed to apply calculations on.
    */
    void calculateMDV(
      const Feature& measured_feature, 
      Feature& normalized_feature,
      const MassIntensityType& mass_intensity_type, 
      const std::string& feature_name);
    
    /**
      @brief This function calculates the mass distribution vector (MDV)
      either normalized to the highest mass intensity (norm_max) or normalized
      to the sum of all mass intensities (norm_sum)
     
      @param[in]   measured_fm FeatureMap with measured intensity for each component.
      @param[out]  normalized_fm FeatureMap with normalized values for each component.
      @param[in]   mass_intensity_type Mass intensity type (either norm_max or norm_sum).
      @param[in]   feature_name Feature name to determine which features are needed to apply calculations on.
    */
    void calculateMDVs(
      const FeatureMap& measured_fm, 
      FeatureMap& normalized_fm,
      const MassIntensityType& mass_intensity_type, 
      const std::string& feature_name);
    
  protected:
    /// Synchronize members with param class
    void updateMembers_() override;
  };
}
