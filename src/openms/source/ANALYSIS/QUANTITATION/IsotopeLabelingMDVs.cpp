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
    const std::vector<std::vector<double>> correction_matrix,
    const std::string correction_matrix_agent)
  {
    // MDV_corrected = correction_matrix_inversed * MDV_observed (normalized_features)
    
    std::vector<std::vector<double>> selected_correction_matrix;
    auto correction_matrix_search = correction_matrices_.find(correction_matrix_agent);
    
    if (!correction_matrix_agent.empty() && correction_matrix_search != correction_matrices_.end())
    {
      selected_correction_matrix = correction_matrix_search->second ;
    }
    else
    {
      selected_correction_matrix = correction_matrix;
    }
    
    uint16_t correction_matrix_n = selected_correction_matrix.size();
    std::vector<std::vector<double>> correction_matrix_inversed(correction_matrix_n, std::vector<double>(correction_matrix_n,0));
    
    // 1- correction matrix inversion
    Eigen::MatrixXd correction_matrix_eigen(selected_correction_matrix.size(), selected_correction_matrix[0].size());
    for (uint i = 0; i < selected_correction_matrix.size(); i++){
      for (uint j = 0; j < selected_correction_matrix[0].size(); j++){
        correction_matrix_eigen(i,j) = selected_correction_matrix[i][j];
      }
    }
    Eigen::MatrixXd correction_matrix_eigen_inversed = correction_matrix_eigen.inverse();
    for (int i = 0; i < correction_matrix_eigen_inversed.rows(); i++){
      for (int j = 0; j < correction_matrix_eigen_inversed.cols(); j++){
        correction_matrix_inversed[i][j] = correction_matrix_eigen_inversed(i,j);
      }
    }
    
    // 2- element-wise expansion with MDV_observed
    std::vector<Feature> normalized_feature_subordinates = normalized_feature.getSubordinates();
    std::vector<double> MDV_observed;
    for (auto it = normalized_feature_subordinates.begin(); it != normalized_feature_subordinates.end(); it++)
    {
      MDV_observed.push_back(it->getMetaValue("peak_apex_int"));
    }
    
    corrected_feature = normalized_feature;
    for (size_t i = 0; i < correction_matrix_inversed.size(); ++i) {
      double corrected_value = 0.0;
      for (size_t j = 0; j < correction_matrix_inversed[0].size(); ++j) {
        corrected_value += correction_matrix_inversed[i][j] * MDV_observed[j];
      }
      corrected_feature.getSubordinates().at(i).setIntensity(corrected_value);
    }
  }

  void IsotopeLabelingMDVs::isotopicCorrections(
    const FeatureMap& normalized_featureMap,
    FeatureMap& corrected_featureMap,
    const std::vector<std::vector<double>> correction_matrix,
    const std::string correction_matrix_agent)
  {
    for (const Feature& feature : normalized_featureMap) {
      Feature corrected_feature;
      isotopicCorrection( feature, corrected_feature, correction_matrix, correction_matrix_agent);
      corrected_featureMap.push_back(corrected_feature);
    }
  }

  void IsotopeLabelingMDVs::calculateIsotopicPurity(
    const Feature& normalized_featuremap,
    Feature& featuremap_with_isotopic_purity,
    std::vector<double>& experiment_data,
    std::string& isotopic_purity_name)
  {
    featuremap_with_isotopic_purity = normalized_featuremap;
    
    if ( !experiment_data.empty() )
    {
      double experiment_data_peak = 0.0;
      std::vector<double>::iterator max_it = std::max_element(experiment_data.begin(), experiment_data.end());
      uint64_t experiment_data_peak_idx = std::distance(experiment_data.begin(), max_it);
      experiment_data_peak = experiment_data[experiment_data_peak_idx];
      
      if ( experiment_data_peak_idx >= 1 && experiment_data_peak != 0.0)
      {
        double previous_experiment_data_peak = experiment_data[experiment_data_peak_idx - 1];
        double isotopic_purity = experiment_data_peak_idx / (experiment_data_peak_idx + ( previous_experiment_data_peak / experiment_data_peak));
        featuremap_with_isotopic_purity.setMetaValue(isotopic_purity_name, isotopic_purity);
      }
    }
  }

  void IsotopeLabelingMDVs::calculateIsotopicPurities(
    const FeatureMap& normalized_featureMap,
    FeatureMap& featureMap_with_isotopic_purity,
    std::vector<double>& experiment_data,
    std::string& isotopic_purity_name)
  {
    for (const Feature& feature : normalized_featureMap){
      Feature feature_with_isotopic_purity;
      calculateIsotopicPurity( feature, feature_with_isotopic_purity, experiment_data, isotopic_purity_name);
      featureMap_with_isotopic_purity.push_back(feature_with_isotopic_purity);
    }
  }
  
  void IsotopeLabelingMDVs::calculateMDVAccuracy(
    const Feature& normalized_feature,
    Feature& feature_with_accuracy_info,
    const std::vector<double>& fragment_isotopomer_measured,
    const std::vector<double>& fragment_isotopomer_theoretical)
  {
    feature_with_accuracy_info = normalized_feature;
    
    std::vector<double> fragment_isotopomer_abs_diff;
    for (size_t i = 0; i < fragment_isotopomer_theoretical.size(); ++i) {
      fragment_isotopomer_abs_diff.push_back(std::abs(fragment_isotopomer_theoretical[i] - fragment_isotopomer_measured[i]));
    }
    
    double diff_mean = 0.0;
    for (size_t i = 0; i < fragment_isotopomer_abs_diff.size(); ++i) {
      diff_mean += fragment_isotopomer_abs_diff.at(i);
    }
    
    diff_mean /= fragment_isotopomer_abs_diff.size();
    
    std::vector<double> fragment_isotopomer_abs_diff_;
    for (size_t i = 0; i < fragment_isotopomer_abs_diff.size(); ++i) {
      fragment_isotopomer_abs_diff_.push_back(std::abs(fragment_isotopomer_abs_diff[i] - diff_mean));
    }
    
    diff_mean = 0.0;
    for (size_t i = 0; i < fragment_isotopomer_abs_diff_.size(); ++i) {
      diff_mean += fragment_isotopomer_abs_diff_.at(i);
    }
    
    diff_mean /= fragment_isotopomer_abs_diff_.size();
    
    feature_with_accuracy_info = normalized_feature;
    feature_with_accuracy_info.setMetaValue("average_accuracy", diff_mean);
  }

  void IsotopeLabelingMDVs::calculateMDVAccuracies(
    const FeatureMap& normalized_featureMap,
    FeatureMap& featureMap_with_accuracy_info,
    const std::vector<double>& fragment_isotopomer_measured,
    const std::vector<double>& fragment_isotopomer_theoretical)
  {
    for (const Feature& feature : normalized_featureMap) {
      Feature feature_with_accuracy_info;
      calculateMDVAccuracy(feature, feature_with_accuracy_info, fragment_isotopomer_measured, fragment_isotopomer_theoretical);
      featureMap_with_accuracy_info.push_back(feature_with_accuracy_info);
    }
  }

  void IsotopeLabelingMDVs::calculateMDV(
    const Feature& measured_feature,
    Feature& normalized_feature,
    const String& mass_intensity_type,
    const String& feature_name)
  {
    std::vector<Feature> measured_feature_subordinates = measured_feature.getSubordinates();
    normalized_feature = measured_feature;
    if (mass_intensity_type == "norm_max")
    {
      if (feature_name == "intensity")
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getIntensity());
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
      
        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if(measured_feature_max != 0.0) {
            normalized_feature.getSubordinates().at(i).setIntensity(normalized_feature.getSubordinates().at(i).getIntensity() /  measured_feature_max);
          }
        }
      }
      
      else
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getMetaValue(feature_name));
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
          
        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if (measured_feature_max != 0.0) {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / measured_feature_max);
          }
        }
      }
    }
    
    else if (mass_intensity_type == "norm_sum")
    {
      if (feature_name == "intensity")
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += it->getIntensity();
        }
        
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_feature.setMetaValue((it - measured_feature_subordinates.begin()), (it->getIntensity() / feature_peak_apex_intensity_sum));
        }
      }
      
      else
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += (Peak2D::IntensityType)it->getMetaValue(feature_name);
        }

        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if (feature_peak_apex_intensity_sum != 0.0)
          {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / feature_peak_apex_intensity_sum);
          }
        }
      }
    }
  }

  void IsotopeLabelingMDVs::calculateMDVs(
    const FeatureMap& measured_featureMap, FeatureMap& normalized_featureMap,
    const String& mass_intensity_type, const String& feature_name)
  {
    for (const Feature& feature : measured_featureMap) {
      Feature normalized_feature;
      calculateMDV(feature, normalized_feature, mass_intensity_type, feature_name);
      normalized_featureMap.push_back(normalized_feature);
    }
  }

  const std::unordered_map<std::string, std::vector<std::vector<double>> > IsotopeLabelingMDVs::correction_matrices_ =
  std::unordered_map<std::string, std::vector<std::vector<double>> >
  {
    std::unordered_map<std::string, std::vector<std::vector<double>> >
    {
      { "TBDMS", {{0.8213, 0.1053, 0.0734, 0.0000},
                  {0.8420, 0.0963, 0.0617, 0.0000},
                  {0.8466, 0.0957, 0.0343, 0.0233},
                  {0.8484, 0.0954, 0.0337, 0.0225}}
      }
    }
  };
  
} // namespace
