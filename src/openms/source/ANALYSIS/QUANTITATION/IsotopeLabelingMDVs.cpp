// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ahmed Khalil $
// $Authors: Ahmed Khalil $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/CONCEPT/LogStream.h>

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
  const std::string IsotopeLabelingMDVs::NamesOfDerivatizationAgent[] = {"NOT_SELECTED", "TBDMS"};
  
  const std::string IsotopeLabelingMDVs::NamesOfMassIntensityType[] = {"norm_max", "norm_sum"};

  IsotopeLabelingMDVs::IsotopeLabelingMDVs() :
    DefaultParamHandler("IsotopeLabelingMDVs")
  {
  }

  IsotopeLabelingMDVs::~IsotopeLabelingMDVs() = default;

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
    auto& em = correction_matrix.getEigenMatrix();
    if (em.isIdentity() && !(em.size() == 0))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "IsotopeLabelingMDVs: The given isotope correction matrix is an identity matrix leading to no correction."
                                        "Please provide a valid correction_matrix.");
    }
    
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
      for (long int i = 0; i < correction_matrix.rows(); ++i)
      {
        for (long int j = 0; j < correction_matrix.cols(); ++j)
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

    // Update MDV_observed to match the size of inversed correction_matrix
    if (static_cast<unsigned long>(correction_matrix_eigen_inversed.cols()) > MDV_observed.size())
    {
      size_t resize_diff = correction_matrix_eigen_inversed.cols() - MDV_observed.size();
      for (size_t i = 0; i < resize_diff; ++i)
      {
        MDV_observed.push_back(0.0);
      }
    }
    // Expand the inversed correction_matrix to be of an equivalent size to MDV_observed
    else if (MDV_observed.size() > static_cast<unsigned long>(correction_matrix_eigen_inversed.cols()))
    {
      size_t resize_diff = MDV_observed.size() - correction_matrix_eigen_inversed.cols();
      
      correction_matrix_eigen_inversed.conservativeResize(correction_matrix_eigen_inversed.rows() + resize_diff,
                                                          correction_matrix_eigen_inversed.cols() + resize_diff);
      
      for (int i = correction_matrix_eigen_inversed.rows() - resize_diff; i < correction_matrix_eigen_inversed.rows(); ++i)
      {
        for (int j = correction_matrix_eigen_inversed.cols() - resize_diff; j < correction_matrix_eigen_inversed.cols(); ++j)
        {
          correction_matrix_eigen_inversed(i,j) = 0.0;
        }
      }
    }

    for (int i = 0; i < correction_matrix_eigen_inversed.rows(); ++i)
    {
      double corrected_value = 0.0;
      for (int j = 0; j < correction_matrix_eigen_inversed.cols(); ++j)
      {
        corrected_value += correction_matrix_eigen_inversed(i,j) * MDV_observed[j];
      }
      corrected_feature.getSubordinates().at(i).setMetaValue("peak_apex_int", std::isnan(corrected_value) ? 0.0 : corrected_value);
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
    Feature& normalized_featuremap,
    const std::vector<double>& experiment_data,
    const std::string& isotopic_purity_name)
  {
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
        normalized_featuremap.setMetaValue(isotopic_purity_name, isotopic_purity);
      }
    }
  }

  void IsotopeLabelingMDVs::calculateIsotopicPurities(
    FeatureMap& normalized_featureMap,
    const std::vector<std::vector<double>>& experiment_data,
    const std::vector<std::string>& isotopic_purity_names)
  {
    for (size_t feature_idx = 0; feature_idx < normalized_featureMap.size(); ++feature_idx)
    {
      calculateIsotopicPurity(normalized_featureMap.at(feature_idx), experiment_data.at(feature_idx), isotopic_purity_names.at(feature_idx));
    }
  }
  
  void IsotopeLabelingMDVs::calculateMDVAccuracy(
    Feature& normalized_feature,
    const std::string& feature_name,
    const std::string& fragment_isotopomer_theoretical_formula)
  {
    std::vector<double> fragment_isotopomer_theoretical, fragment_isotopomer_measured;
    
    for (auto it = normalized_feature.getSubordinates().begin(); it != normalized_feature.getSubordinates().end(); it++)
    {
      if (feature_name == "intensity")
      {
        fragment_isotopomer_measured.push_back((double)(it->getIntensity()));
      }
      else if (feature_name != "intensity" && it->metaValueExists(feature_name))
      {
        fragment_isotopomer_measured.push_back(it->getMetaValue(feature_name));
      }
    }
    
    if (normalized_feature.getSubordinates().size() != fragment_isotopomer_measured.size() || fragment_isotopomer_measured.empty()) {
      OpenMS_Log_fatal << "Missing values for the Measured Isotopomer Fragment, Please make sure the Subordinates are accordingly updated." << std::endl;
    }
    
    // Generate theoretical values with the exact same length as fragment_isotopomer_measured
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
    
      normalized_feature.setMetaValue("average_accuracy", diff_mean);
      
      for (size_t feature_subordinate = 0; feature_subordinate < normalized_feature.getSubordinates().size(); ++feature_subordinate)
      {
        normalized_feature.getSubordinates().at(feature_subordinate).setMetaValue("absolute_difference", fragment_isotopomer_abs_diff.at(feature_subordinate));
      }
    }
  }

  void IsotopeLabelingMDVs::calculateMDVAccuracies(
    FeatureMap& normalized_featureMap,
    const std::string& feature_name,
    const std::map<std::string, std::string>& fragment_isotopomer_theoretical_formulas)
  {
    for (size_t feature_idx = 0; feature_idx < normalized_featureMap.size(); ++feature_idx)
    {
      if (normalized_featureMap.at(feature_idx).metaValueExists("PeptideRef"))
      {
        calculateMDVAccuracy(normalized_featureMap.at(feature_idx),
                             feature_name,
                             fragment_isotopomer_theoretical_formulas.find((std::string)normalized_featureMap.at(feature_idx).getMetaValue("PeptideRef"))->second);
      }
      else
      {
        OPENMS_LOG_ERROR << "No PeptideRef in FeatureMap (MetaValue doesn't exist)!" << std::endl;
      }
    }
  }

  void IsotopeLabelingMDVs::calculateMDV(
    const Feature& measured_feature,
    Feature& normalized_feature,
    const MassIntensityType& mass_intensity_type,
    const std::string& feature_name)
  {
    std::vector<Feature> measured_feature_subordinates = measured_feature.getSubordinates();
    normalized_feature = measured_feature;
    if (mass_intensity_type == MassIntensityType::NORM_MAX)
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
      
        if (measured_feature_max != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setIntensity(normalized_feature.getSubordinates().at(i).getIntensity() /  measured_feature_max);
          }
        }
      }
      // for every other case where feature_name isn't 'intensity', i.e. 'peak_apex_int'
      else
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getMetaValue(feature_name));
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
          
        if (measured_feature_max != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setMetaValue(feature_name, (OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / measured_feature_max);
          }
        }
      }
    }
    else if (mass_intensity_type == MassIntensityType::NORM_SUM)
    {
      if (feature_name == "intensity")
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += it->getIntensity();
        }
        
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_feature.getSubordinates().at(it - measured_feature_subordinates.begin()).setIntensity((it->getIntensity() / feature_peak_apex_intensity_sum));
        }
      }
      // for every other case where feature_name isn't 'intensity', i.e. 'peak_apex_int'
      else
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += (Peak2D::IntensityType)it->getMetaValue(feature_name);
        }

        if (feature_peak_apex_intensity_sum != 0.0)
        {
          for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
          {
            normalized_feature.getSubordinates().at(i).setMetaValue(feature_name, (OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / feature_peak_apex_intensity_sum);
          }
        }
      }
    }
  }

  void IsotopeLabelingMDVs::calculateMDVs(
    const FeatureMap& measured_featureMap, FeatureMap& normalized_featureMap,
    const MassIntensityType& mass_intensity_type, const std::string& feature_name)
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
