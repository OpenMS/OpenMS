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
// $Maintainer: Ahmed Khalil $
// $Authors: Ahmed Khalil $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

//Kernal classes
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
//#include <unordered_map>
#include <algorithm>
#include <Eigen/Dense>

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
    ~IsotopeLabelingMDVs();
    //@}
    
    enum class DerivatizationAgent
    {
      NOT_SELECTED,
      TBDMS,
      SIZE_OF_DERIVATIZATIONAGENT
    };
    
    enum class FeatureName
    {
      INTENSITY,
      PEAK_APEX_INT,
      SIZE_OF_FEATURENAME
    };
    
    enum class MassIntensityType
    {
      NORM_MAX,
      NORM_SUM,
      SIZE_OF_MASSINTENSITYTYPE
    };
    
    static const std::string NamesOfFeatureName[static_cast<int>(FeatureName::SIZE_OF_FEATURENAME)];
 
    /**
      @brief This function performs an isotopic correction to account for unlabeled abundances coming from
      the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
      MDV_corrected = CM^-1 * MDV_normalized_feature
      The formula is obtained from "The importance of accurately correcting for the natural abundance of stable isotopes",
      Midani et al, doi:10.1016/j.ab.2016.12.011
     
      @param[in]  normalized_feature Feature with normalized values for each component and unlabeled chemical formula for each component group.
      @param[in]  correction_matrix  Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of
      naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent.
      @param[in]  correction_matrix_agent name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied,
      only "TBDMS" is supported for now.
      @param[out] corrected_feature Feature with corrected values for each component.
    */
    void isotopicCorrection(
      const Feature& normalized_feature, Feature& corrected_feature,
      const Matrix<double>& correction_matrix, const DerivatizationAgent& correction_matrix_agent);
    
    /**
      @brief This function performs an isotopic correction to account for unlabeled abundances coming from
      the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
      MDV_corrected = CM^-1 * MDV_normalized_feature
      The formula is obtained from "The importance of accurately correcting for the natural abundance of stable isotopes",
      Midani et al, doi:10.1016/j.ab.2016.12.011
     
      @param[in]  normalized_featuremap FeatureMap with normalized values for each component and unlabeled chemical formula for each component group.
      @param[in]  correction_matrix Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of
      naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent.
      @param[in]  correction_matrix_agent name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied,
      only "TBDMS" is supported for now.
      @param[out] corrected_featuremap FeatureMap with corrected values for each component.
    */
    void isotopicCorrections(
      const FeatureMap& normalized_featureMap, FeatureMap& corrected_featureMap,
      const Matrix<double>& correction_matrix, const DerivatizationAgent& correction_matrix_agent);

    /**
      @brief This function calculates the isotopic purity of the MDV using the following formula:
      isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
      where n in M+n is represented as the index of the result.
      The formula is extracted from "High-resolution 13C metabolic flux analysis",
      Long et al, doi:10.1038/s41596-019-0204-0

      @param[in]  normalized_feature Feature with normalized values for each component and the number of heavy labeled e.g., carbons.
      @param[in]  experiment_data Experiment data in percent.
      @param[in]  isotopic_purity_name Name of the isotopic purity tracer to be saved as a meta value.
      @param[out] feature_with_isotopic_purity Feature with the calculated isotopic purity for the component group.
    */
    void calculateIsotopicPurity(
      const Feature& normalized_feature, Feature& feature_with_isotopic_purity,
      const std::vector<double>& experiment_data, const std::string& isotopic_purity_name);
    
    /**
      @brief This function calculates the isotopic purity of the MDVs using the following formula:
      isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
      where n in M+n is represented as the index of the result.
      The formula is extracted from "High-resolution 13C metabolic flux analysis",
      Long et al, doi:10.1038/s41596-019-0204-0

      @param[in]  normalized_featureMap FeatureMap with normalized values for each component and the number of heavy labeled e.g., carbons.
      @param[in]  experiment_data Experiment data in percent.
      @param[in]  isotopic_purity_name Name of the isotopic purity tracer to be saved as a meta value.
      @param[out] featureMap_with_isotopic_purity FeatureMap with the calculated isotopic purity for the component group.
    */
    void calculateIsotopicPurities(
      const FeatureMap& normalized_featureMap, FeatureMap& featureMap_with_isotopic_purity,
      const std::vector<double>& experiment_data, const std::string& isotopic_purity_name);

    /**
      @brief This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
      using average deviation to the mean. The result is mapped to the meta value "average_accuracy" in the updated Feature.
   
      @param[in]  normalized_feature Feature with normalized values for each component and the chemical formula of the component group.
      @param[in]  fragment_isotopomer_measured Measured scan values.
      @param[in]  fragment_isotopomer_theoretical_formula Empirical formula from which the theoretical values will be generated.
      @param[out] feature_with_accuracy_info Feature with the component group accuracy and accuracy for the error for each component.
    */
    void calculateMDVAccuracy(
      const Feature& normalized_feature, Feature& feature_with_accuracy_info,
      const std::vector<double>& fragment_isotopomer_measured, const std::string& fragment_isotopomer_theoretical_formula);
    
    /**
       @brief This function calculates the accuracy of the MDVs as compared to the theoretical MDVs (only for 12C quality control experiments)
       using average deviation to the mean.
    
       @param[in]  normalized_featuremap FeatureMap with normalized values for each component and the chemical formula of the component group.
       @param[in]  fragment_isotopomer_measured Measured scan values.
       @param[in]  fragment_isotopomer_theoretical_formula Empirical formula from which the theoretical values will be generated.
       @param[out] featuremap_with_accuracy_info FeatureMap with the component group accuracy and accuracy for the error for each component.
    */
    void calculateMDVAccuracies(
      const FeatureMap& normalized_featureMap, FeatureMap& featureMap_with_accuracy_info,
      const std::vector<double>& fragment_isotopomer_measured, const std::string& fragment_isotopomer_theoretical_formula);
 
    /**
      @brief This function calculates the mass distribution vector (MDV)
      either normalized to the highest mass intensity (norm_max) or normalized
      to the sum of all mass intensities (norm_sum)
     
      @param[in]   mass_intensity_type Mass intensity type (either norm_max or norm_sum).
      @param[in]   feature_name Feature name to determine which features are needed to apply calculations on.
      @param[in]   measured_feature Feature with measured intensity for each component.
      @param[out]  normalized_feature Feature with normalized values for each component.
    */
    void calculateMDV(
      const Feature& measured_feature, Feature& normalized_feature,
      const MassIntensityType& mass_intensity_type, const FeatureName& feature_name);
    
    /**
      @brief This function calculates the mass distribution vector (MDV)
      either normalized to the highest mass intensity (norm_max) or normalized
      to the sum of all mass intensities (norm_sum)
     
      @param[in]   mass_intensity_type Mass intensity type (either norm_max or norm_sum).
      @param[in]   feature_name Feature name to determine which features are needed to apply calculations on.
      @param[in]   measured_featuremap FeatureMap with measured intensity for each component.
      @param[out]  normalized_featuremap FeatureMap with normalized values for each component.
    */
    void calculateMDVs(
      const FeatureMap& measured_featureMap, FeatureMap& normalized_featureMap,
      const MassIntensityType& mass_intensity_type, const FeatureName& feature_name);
    
  protected:
    /// Synchronize members with param class
    void updateMembers_() override;
  };
}
