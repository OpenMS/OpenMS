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

#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>

namespace OpenMS
{

  const std::string IsotopeLabelingMDVs::NamesOfFeatureName[] = {"intensity", "peak_apex_int"};

  IsotopeLabelingMDVs::IsotopeLabelingMDVs() :
    DefaultParamHandler("IsotopeLabelingMDVs")
  {
  }

  IsotopeLabelingMDVs::~IsotopeLabelingMDVs()
  {
  }

  void IsotopeLabelingMDVs::updateMembers_()
  {
  }
  
  void IsotopeLabelingMDVs::isotopicCorrection(
    const Feature& normalized_feature,
    Feature& corrected_feature,
    const Matrix<double>& correction_matrix,
    const DerivatizationAgent& correction_matrix_agent)
  {
    // MDV_corrected = correction_matrix_inversed * MDV_observed (normalized_features)
    
    /// Correction Matrices for various derivatization agents
    const std::map<DerivatizationAgent, std::vector<std::vector<double>> > correction_matrices =
    {
      { DerivatizationAgent::TBDMS, {{0.8213, 0.1053, 0.0734, 0.0000},
                                     {0.8420, 0.0963, 0.0617, 0.0000},
                                     {0.8466, 0.0957, 0.0343, 0.0233},
                                     {0.8484, 0.0954, 0.0337, 0.0225}}
      }
    };
    
    Eigen::MatrixXd correction_matrix_eigen;
    auto correction_matrix_search = correction_matrices.find(correction_matrix_agent);
    
    // use the internally saved correction matrix if the derivatization agent name is supplied and found
    if (correction_matrix_agent != DerivatizationAgent::NOT_SELECTED && correction_matrix_search != correction_matrices.end())
    {
      correction_matrix_eigen.resize(correction_matrix_search->second.size(), correction_matrix_search->second[0].size());
      for (size_t i = 0; i < correction_matrix_search->second.size(); ++i)
      {
        for (size_t j = 0; j < correction_matrix_search->second[0].size(); ++j)
        {
          correction_matrix_eigen(i,j) = correction_matrix_search->second[i][j];
        }
      }
    }
    // copy correction_matrix to an eigen matrix, when correction_matrix is supplied
    else
    {
      correction_matrix_eigen.resize(correction_matrix.rows(), correction_matrix.cols());
      for (size_t i = 0; i < correction_matrix.rows(); ++i)
      {
        for (size_t j = 0; j < correction_matrix.cols(); ++j)
        {
          correction_matrix_eigen(i,j) = correction_matrix(i,j);
        }
      }
    }
    
    // 1- inversion of correction matrix
    Eigen::MatrixXd correction_matrix_eigen_inversed = correction_matrix_eigen.inverse();
    
    // 2- element-wise expansion with MDV_observed
    std::vector<double> MDV_observed;
    for (const auto& feature : normalized_feature.getSubordinates())
    {
      MDV_observed.push_back(feature.getMetaValue("peak_apex_int"));
    }
    
    corrected_feature = normalized_feature;
    for (int i = 0; i < correction_matrix_eigen_inversed.rows(); ++i)
    {
      double corrected_value = 0.0;
      for (int j = 0; j < correction_matrix_eigen_inversed.cols(); ++j)
      {
        corrected_value += correction_matrix_eigen_inversed(i,j) * MDV_observed[j];
      }
      corrected_feature.getSubordinates().at(i).setIntensity(corrected_value);
    }
  }

  void IsotopeLabelingMDVs::isotopicCorrections(
    const FeatureMap& normalized_featureMap,
    FeatureMap& corrected_featureMap,
    const Matrix<double>& correction_matrix,
    const DerivatizationAgent& correction_matrix_agent)
  {
    for (const Feature& feature : normalized_featureMap)
    {
      Feature corrected_feature;
      isotopicCorrection(feature, corrected_feature, correction_matrix, correction_matrix_agent);
      corrected_featureMap.push_back(corrected_feature);
    }
  }

  void IsotopeLabelingMDVs::calculateIsotopicPurity(
    const Feature& normalized_featuremap,
    Feature& featuremap_with_isotopic_purity,
    const std::vector<double>& experiment_data,
    const std::string& isotopic_purity_name)
  {
    featuremap_with_isotopic_purity = normalized_featuremap;
    
    if (!experiment_data.empty())
    {
      double experiment_data_peak = 0.0;
      std::vector<double> experiment_data_ = experiment_data;
      std::vector<double>::iterator max_it = std::max_element(experiment_data_.begin(), experiment_data_.end());
      uint64_t experiment_data_peak_idx = std::distance(experiment_data_.begin(), max_it);
      experiment_data_peak = experiment_data_[experiment_data_peak_idx];
      
      if (experiment_data_peak_idx >= 1 && experiment_data_peak != 0.0)
      {
        double previous_experiment_data_peak = experiment_data[experiment_data_peak_idx - 1];
        double isotopic_purity = experiment_data_peak_idx / (experiment_data_peak_idx + (previous_experiment_data_peak / experiment_data_peak));
        featuremap_with_isotopic_purity.setMetaValue(isotopic_purity_name, isotopic_purity);
      }
    }
  }

  void IsotopeLabelingMDVs::calculateIsotopicPurities(
    const FeatureMap& normalized_featureMap,
    FeatureMap& featureMap_with_isotopic_purity,
    const std::vector<double>& experiment_data,
    const std::string& isotopic_purity_name)
  {
    for (const Feature& feature : normalized_featureMap)
    {
      Feature feature_with_isotopic_purity;
      calculateIsotopicPurity(feature, feature_with_isotopic_purity, experiment_data, isotopic_purity_name);
      featureMap_with_isotopic_purity.push_back(feature_with_isotopic_purity);
    }
  }
  
  void IsotopeLabelingMDVs::calculateMDVAccuracy(
    const Feature& normalized_feature,
    Feature& feature_with_accuracy_info,
    const std::vector<double>& fragment_isotopomer_measured,
    const std::string& fragment_isotopomer_theoretical_formula)
  {
    feature_with_accuracy_info = normalized_feature;
    
    std::vector<double> fragment_isotopomer_theoretical;
    
    IsotopeDistribution theoretical_iso(EmpiricalFormula(fragment_isotopomer_theoretical_formula).getIsotopeDistribution(CoarseIsotopePatternGenerator(fragment_isotopomer_measured.size())));
    for (IsotopeDistribution::ConstIterator it = theoretical_iso.begin(); it != theoretical_iso.end(); ++it)
    {
      fragment_isotopomer_theoretical.push_back( it->getIntensity() );
    }
    
    std::vector<double> fragment_isotopomer_abs_diff;
    for (size_t i = 0; i < fragment_isotopomer_theoretical.size(); ++i)
    {
      fragment_isotopomer_abs_diff.push_back(std::fabs(fragment_isotopomer_theoretical[i] - fragment_isotopomer_measured[i]));
    }
    
    if (!fragment_isotopomer_abs_diff.empty())
    {
      double diff_mean = OpenMS::Math::mean(fragment_isotopomer_abs_diff.begin(), fragment_isotopomer_abs_diff.end());
    
      diff_mean = OpenMS::Math::MeanAbsoluteDeviation(fragment_isotopomer_abs_diff.begin(), fragment_isotopomer_abs_diff.end(), diff_mean);
    
      feature_with_accuracy_info.setMetaValue("average_accuracy", diff_mean);
    }
  }

  void IsotopeLabelingMDVs::calculateMDVAccuracies(
    const FeatureMap& normalized_featureMap,
    FeatureMap& featureMap_with_accuracy_info,
    const std::vector<double>& fragment_isotopomer_measured,
    const std::string& fragment_isotopomer_theoretical_formula)
  {
    featureMap_with_accuracy_info.clear();
    for (const Feature& feature : normalized_featureMap)
    {
      Feature feature_with_accuracy_info;
      calculateMDVAccuracy(feature, feature_with_accuracy_info, fragment_isotopomer_measured, fragment_isotopomer_theoretical_formula);
      featureMap_with_accuracy_info.push_back(feature_with_accuracy_info);
    }
  }

  void IsotopeLabelingMDVs::calculateMDV(
    const Feature& measured_feature,
    Feature& normalized_feature,
    const MassIntensityType& mass_intensity_type,
    const FeatureName& feature_name)
  {
    std::vector<Feature> measured_feature_subordinates = measured_feature.getSubordinates();
    normalized_feature = measured_feature;
    if (mass_intensity_type == MassIntensityType::NORM_MAX)
    {
      if (feature_name == FeatureName::INTENSITY)
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getIntensity());
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
      
        if (measured_feature_max != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setIntensity(normalized_feature.getSubordinates().at(i).getIntensity() /  measured_feature_max);
          }
        }
      }
      // for every other case where feature_name isn't 'intensity', i.e. 'peak_apex_int'
      else if (feature_name == FeatureName::PEAK_APEX_INT)
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getMetaValue(NamesOfFeatureName[static_cast<int>(FeatureName::PEAK_APEX_INT)]));
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
          
        if (measured_feature_max != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(NamesOfFeatureName[static_cast<int>(FeatureName::PEAK_APEX_INT)]) / measured_feature_max);
          }
        }
      }
    }
    else if (mass_intensity_type == MassIntensityType::NORM_SUM)
    {
      if (feature_name == FeatureName::INTENSITY)
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += it->getIntensity();
        }
        
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_feature.setMetaValue((it - measured_feature_subordinates.begin()), (it->getIntensity() / feature_peak_apex_intensity_sum));
        }
      }
      // for every other case where feature_name isn't 'intensity', i.e. 'peak_apex_int'
      else if (feature_name == FeatureName::PEAK_APEX_INT)
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += (Peak2D::IntensityType)it->getMetaValue(NamesOfFeatureName[static_cast<int>(FeatureName::PEAK_APEX_INT)]);
        }

        if (feature_peak_apex_intensity_sum != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(NamesOfFeatureName[static_cast<int>(FeatureName::PEAK_APEX_INT)]) / feature_peak_apex_intensity_sum);
          }
        }
      }
    }
  }

  void IsotopeLabelingMDVs::calculateMDVs(
    const FeatureMap& measured_featureMap, FeatureMap& normalized_featureMap,
    const MassIntensityType& mass_intensity_type, const FeatureName& feature_name)
  {
    normalized_featureMap.clear();
    
    for (const Feature& feature : measured_featureMap)
    {
      Feature normalized_feature;
      calculateMDV(feature, normalized_feature, mass_intensity_type, feature_name);
      normalized_featureMap.push_back(normalized_feature);
    }
  }
} // namespace
