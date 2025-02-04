// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ahmed Khalil $
// $Authors: Ahmed Khalil $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <cassert>

using namespace OpenMS;

START_TEST(IsotopeLabelingMDVs, "$Id$")

IsotopeLabelingMDVs* ptr          = nullptr;
IsotopeLabelingMDVs* nullPointer  = nullptr;

START_SECTION((IsotopeLabelingMDVs()))
  ptr = new IsotopeLabelingMDVs();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IsotopeLabelingMDVs()))
  delete ptr;
END_SECTION

START_SECTION(( void IsotopeLabelingMDVs::calculateMDV(
                                              const Feature& measured_feature,
                                              Feature& normalized_featuremap,
                                              const String& mass_intensity_type,
                                              const String& feature_name) ))

  // case 1:  intensity with norm max and norm sum (x)  : intensity (peak area) not supplied
  // case 2:  peak apex with norm max and norm sum      : - Lactate1 & Lactate2 - peak_apex_int - norm_max
  //                                                      - Lactate1 & Lactate2 - peak_apex_int - norm_sum

  IsotopeLabelingMDVs                isotopelabelingmdvs;
  
  // From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  std::vector<double> L1_peak_apex_int {3.61e+08, 1.20e+04, 1.02e+05, 2.59e+04};
  std::vector<double> L2_peak_apex_int {2.77e+07, 5.45e+04, 6.26e+05, 7.46e+04, 2.75e+04};

  std::vector<double> L1_norm_max {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<double> L1_norm_sum {9.9961e-01, 3.3228e-05, 2.8243e-04, 7.1717e-05};

  std::vector<double> L2_norm_max {1.00e+00, 1.967e-03, 2.259e-02, 2.693e-03, 9.927e-04};
  std::vector<double> L2_norm_sum {9.7252e-01, 1.9134e-03, 2.1978e-02, 2.6191e-03, 9.655e-04};

  // Lactate1 & Lactate2 - peak_apex_int - norm_max
  OpenMS::Feature               lactate_1_normmax;
  OpenMS::Feature               lactate_1_normalized_normmax;
  std::vector<OpenMS::Feature>  L1_subordinates_normmax;

  lactate_1_normmax.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_peak_apex_int[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normmax.setSubordinates(L1_subordinates_normmax);
  
  isotopelabelingmdvs.calculateMDV(lactate_1_normmax, lactate_1_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, "peak_apex_int");

  for(size_t i = 0; i < lactate_1_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normmax.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_norm_max.at(i));
  }

  OpenMS::Feature               lactate_2_normmax;
  OpenMS::Feature               lactate_2_normalized_normmax;
  std::vector<OpenMS::Feature>  L2_subordinates_normmax;

  lactate_2_normmax.setMetaValue("PeptideRef", "Lactate2");
  for (uint16_t i = 0; i < L2_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(219+i));
    sub.setMetaValue("peak_apex_int", L2_peak_apex_int[i]);
    L2_subordinates_normmax.push_back(sub);
  }
  lactate_2_normmax.setSubordinates(L2_subordinates_normmax);

  isotopelabelingmdvs.calculateMDV(lactate_2_normmax, lactate_2_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, "peak_apex_int");

  for(size_t i = 0; i < lactate_2_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normmax.getSubordinates().at(i).getMetaValue("peak_apex_int"), L2_norm_max.at(i));
  }

  // Lactate1 & Lactate2 - peak_apex_int - norm_sum
  OpenMS::Feature               lactate_1_normsum;
  OpenMS::Feature               lactate_1_normalized_normsum;
  std::vector<OpenMS::Feature>  L1_subordinates_normsum;

  lactate_1_normsum.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_peak_apex_int[i]);
    L1_subordinates_normsum.push_back(sub);
  }
  lactate_1_normsum.setSubordinates(L1_subordinates_normsum);

  isotopelabelingmdvs.calculateMDV(lactate_1_normsum, lactate_1_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, "peak_apex_int");

  for(size_t i = 0; i < lactate_1_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normsum.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_norm_sum.at(i));
  }

  OpenMS::Feature lactate_2_normsum; OpenMS::Feature lactate_2_normalized_normsum;
  std::vector<OpenMS::Feature> L2_subordinates_normsum;

  lactate_2_normsum.setMetaValue("PeptideRef", "Lactate2");
  for (uint16_t i = 0; i < L2_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(219+i));
    sub.setMetaValue("peak_apex_int", L2_peak_apex_int[i]);
    L2_subordinates_normsum.push_back(sub);
  }
  lactate_2_normsum.setSubordinates(L2_subordinates_normsum);

//  isotopelabelingmdvs.calculateMDV(lactate_2_normsum, lactate_2_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_2_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normsum.getSubordinates().at(i).getIntensity(), L2_norm_sum.at(i));
  }

END_SECTION

START_SECTION(( void IsotopeLabelingMDVs::calculateMDVs(
                                              const FeatureMap& measured_feature,
                                              FeatureMap& normalized_featuremap,
                                              const String& mass_intensity_type,
                                              const String& feature_name) ))

  // case 1:  intensity with norm max and norm sum (x)  : intensity (peak area) not supplied
  // case 2:  peak apex with norm max and norm sum      : - Lactate1 & Lactate2 - peak_apex_int - norm_max
  //                                                      - Lactate1 & Lactate2 - peak_apex_int - norm_sum

  IsotopeLabelingMDVs                isotopelabelingmdvs;
  
  // From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  std::vector<double> L1_peak_apex_int {3.61e+08, 1.20e+04, 1.02e+05, 2.59e+04};
  std::vector<double> L2_peak_apex_int {2.77e+07, 5.45e+04, 6.26e+05, 7.46e+04, 2.75e+04};

  std::vector<double> L1_norm_max {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<double> L1_norm_sum {9.9961e-01, 3.3228e-05, 2.8243e-04, 7.1717e-05};

  std::vector<double> L2_norm_max {1.00e+00, 1.967e-03, 2.259e-02, 2.693e-03, 9.927e-04};
  std::vector<double> L2_norm_sum {9.7252e-01, 1.9134e-03, 2.1978e-02, 2.6191e-03, 9.655e-04};
  

  // Lactate1 & Lactate2 - peak_apex_int - norm_max
  OpenMS::Feature               lactate_1_normmax;
  OpenMS::Feature               lactate_1_normalized_normmax;
  std::vector<OpenMS::Feature>  L1_subordinates_normmax;

  lactate_1_normmax.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_peak_apex_int[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normmax.setSubordinates(L1_subordinates_normmax);
  
  isotopelabelingmdvs.calculateMDV(lactate_1_normmax, lactate_1_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, "peak_apex_int");

  for(size_t i = 0; i < lactate_1_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normmax.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_norm_max.at(i));
  }

  OpenMS::Feature               lactate_2_normmax;
  OpenMS::Feature               lactate_2_normalized_normmax;
  std::vector<OpenMS::Feature>  L2_subordinates_normmax;

  lactate_2_normmax.setMetaValue("PeptideRef", "Lactate2");
  for (uint16_t i = 0; i < L2_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(219+i));
    sub.setMetaValue("peak_apex_int", L2_peak_apex_int[i]);
    L2_subordinates_normmax.push_back(sub);
  }
  lactate_2_normmax.setSubordinates(L2_subordinates_normmax);

  isotopelabelingmdvs.calculateMDV(lactate_2_normmax, lactate_2_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, "peak_apex_int");

  for(size_t i = 0; i < lactate_2_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normmax.getSubordinates().at(i).getMetaValue("peak_apex_int"), L2_norm_max.at(i));
  }

  // Lactate1 & Lactate2 - peak_apex_int - norm_sum
  OpenMS::Feature               lactate_1_normsum;
  OpenMS::Feature               lactate_1_normalized_normsum;
  std::vector<OpenMS::Feature>  L1_subordinates_normsum;

  lactate_1_normsum.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_peak_apex_int[i]);
    L1_subordinates_normsum.push_back(sub);
  }
  lactate_1_normsum.setSubordinates(L1_subordinates_normsum);

  isotopelabelingmdvs.calculateMDV(lactate_1_normsum, lactate_1_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, "peak_apex_int");

  for(size_t i = 0; i < lactate_1_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normsum.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_norm_sum.at(i));
  }

  OpenMS::Feature lactate_2_normsum; OpenMS::Feature lactate_2_normalized_normsum;
  std::vector<OpenMS::Feature> L2_subordinates_normsum;

  lactate_2_normsum.setMetaValue("PeptideRef", "Lactate2");
  for (uint16_t i = 0; i < L2_peak_apex_int.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(219+i));
    sub.setMetaValue("peak_apex_int", L2_peak_apex_int[i]);
    L2_subordinates_normsum.push_back(sub);
  }
  lactate_2_normsum.setSubordinates(L2_subordinates_normsum);

  isotopelabelingmdvs.calculateMDV(lactate_2_normsum, lactate_2_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, "peak_apex_int");

  for(size_t i = 0; i < lactate_2_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normsum.getSubordinates().at(i).getMetaValue("peak_apex_int"), L2_norm_sum.at(i));
  }

  OpenMS::FeatureMap  lactate_normmax;
  OpenMS::FeatureMap  lactate_normsum;
  OpenMS::FeatureMap  lactate_normalized_normmax;
  OpenMS::FeatureMap  lactate_normalized_normsum;

  lactate_normmax.push_back(lactate_1_normmax);
  lactate_normmax.push_back(lactate_2_normmax);

  lactate_normsum.push_back(lactate_1_normsum);
  lactate_normsum.push_back(lactate_2_normsum);

  isotopelabelingmdvs.calculateMDVs(lactate_normmax, lactate_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, "peak_apex_int");
  isotopelabelingmdvs.calculateMDVs(lactate_normsum, lactate_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, "peak_apex_int");

  for(size_t i = 0; i < lactate_normalized_normmax.size(); ++i)
  {
    for(size_t j = 0; j < lactate_normalized_normmax.at(i).getSubordinates().size(); ++j)
    {
      if (i == 0) // lactate_1
      {
        TEST_REAL_SIMILAR(lactate_normalized_normmax.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L1_norm_max.at(j));
      }
      else if (i == 1) // lactate_2
      {
        TEST_REAL_SIMILAR(lactate_normalized_normmax.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L2_norm_max.at(j));
      }
    }
  }

  for(size_t i = 0; i < lactate_normalized_normsum.size(); ++i)
  {
    for(size_t j = 0; j < lactate_normalized_normsum.at(i).getSubordinates().size(); ++j)
    {
      if (i == 0) // lactate_1
      {
        TEST_REAL_SIMILAR(lactate_normalized_normsum.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L1_norm_sum.at(j));
      }
      else if (i == 1) // lactate_2
      {
        TEST_REAL_SIMILAR(lactate_normalized_normsum.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L2_norm_sum.at(j));
      }
    }
  }

END_SECTION

START_SECTION(( void IsotopeLabelingMDVs::isotopicCorrection(
                                              const Feature& normalized_feature,
                                              Feature& corrected_feature,
                                              const Matrix<double>& correction_matrix),
                                              const std::string correction_matrix_agent ))

  // case 1: validating matrix inverse (separately tested)
  // case 2: validating corrected results (corrected peak_apex_int)

  IsotopeLabelingMDVs                   isotopelabelingmdvs;
  OpenMS::Feature                       lactate_1_normalized;
  OpenMS::Feature                       lactate_1_corrected;
  std::vector<std::vector<double>>      correction_matrix_inversed(4, std::vector<double>(4,0));
  // Correction Matrix extracted from "TOOLS FOR MASS ISOTOPOMER DATA EVALUATION IN 13C FLUX ANALYSIS,
  // Wahl et al, P.263, Table I"
  Matrix<double> correction_matrix_tBDMS;
  double correction_matrix_tBDMS_[4][4] = {
    {0.8213, 0.1053, 0.0734, 0.0000},
    {0.8420, 0.0963, 0.0617, 0.0000},
    {0.8466, 0.0957, 0.0343, 0.0233},
    {0.8484, 0.0954, 0.0337, 0.0225}
  };
  correction_matrix_tBDMS.setMatrix<double,4,4>(correction_matrix_tBDMS_);

  // L1_norm_max, L1_peak_apex_int From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // L1_corrected self calculated
  std::vector<double>                   L1_norm_max       {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<double>                   L1_corrected      {-12.7699, 140.7289, -45.3788, -47.2081};
  std::vector<Peak2D::IntensityType>    L1_peak_apex_int  {3.61e+08, 1.20e+04, 1.02e+05, 2.59e+04};
  std::vector<OpenMS::Feature>          L1_subordinates_normmax;

  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  isotopelabelingmdvs.isotopicCorrection(lactate_1_normalized, lactate_1_corrected,
                                         correction_matrix_tBDMS, IsotopeLabelingMDVs::DerivatizationAgent::NOT_SELECTED);

  for(size_t i = 0; i < lactate_1_corrected.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_corrected.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_corrected[i]);
  }

  isotopelabelingmdvs.isotopicCorrection(lactate_1_normalized, lactate_1_corrected,
                                         {}, IsotopeLabelingMDVs::DerivatizationAgent::TBDMS);

  for(size_t i = 0; i < lactate_1_corrected.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_corrected.getSubordinates().at(i).getMetaValue("peak_apex_int"), L1_corrected[i]);
  }

END_SECTION

START_SECTION(( void IsotopeLabelingMDVs::isotopicCorrections(
                                              const FeatureMap& normalized_feature,
                                              FeatureMap& corrected_feature,
                                              const Matrix<double>& correction_matrix,
                                              const std::string correction_matrix_agent) ))

  // case 1: validating corrected results (corrected peak_apex_int)

  IsotopeLabelingMDVs                   isotopelabelingmdvs;
  OpenMS::Feature                       lactate_1_normalized;
  OpenMS::FeatureMap                    lactate_1_featureMap;
  OpenMS::FeatureMap                    lactate_1_corrected_featureMap;
  std::vector<std::vector<double>>      correction_matrix_inversed(4, std::vector<double>(4,0));
  // Correction Matrix extracted from "TOOLS FOR MASS ISOTOPOMER DATA EVALUATION IN 13C FLUX ANALYSIS,
  // Wahl et al, P.263, Table I"
  Matrix<double> correction_matrix_tBDMS;
  double correction_matrix_tBDMS_[4][4] = {
    {0.8213, 0.1053, 0.0734, 0.0000},
    {0.8420, 0.0963, 0.0617, 0.0000},
    {0.8466, 0.0957, 0.0343, 0.0233},
    {0.8484, 0.0954, 0.0337, 0.0225}
  };
  correction_matrix_tBDMS.setMatrix<double, 4, 4>(correction_matrix_tBDMS_);

  // L1_norm_max, L1_peak_apex_int From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // L1_corrected self calculated
  std::vector<double>                   L1_norm_max       {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<double>                   L1_corrected      {-12.7699, 140.7289, -45.3788, -47.2081};
  std::vector<Peak2D::IntensityType>    L1_peak_apex_int  {3.61e+08, 1.20e+04, 1.02e+05, 2.59e+04};
  std::vector<OpenMS::Feature>          L1_subordinates_normmax;

  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  for(uint8_t i = 0; i < 3; ++i)
  {
    lactate_1_featureMap.push_back(lactate_1_normalized);
  }

  isotopelabelingmdvs.isotopicCorrections(lactate_1_featureMap, lactate_1_corrected_featureMap, correction_matrix_tBDMS, IsotopeLabelingMDVs::DerivatizationAgent::NOT_SELECTED);
  for(uint8_t i = 0; i < lactate_1_corrected_featureMap.size(); ++i)
  {
    for(uint8_t j = 0; j < lactate_1_corrected_featureMap.at(i).getSubordinates().size(); ++j)
    {
      TEST_REAL_SIMILAR(lactate_1_corrected_featureMap.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L1_corrected[j]);
    }
  }

  lactate_1_corrected_featureMap.clear();

  isotopelabelingmdvs.isotopicCorrections(lactate_1_featureMap, lactate_1_corrected_featureMap, {}, IsotopeLabelingMDVs::DerivatizationAgent::TBDMS);
  for(uint8_t i = 0; i < lactate_1_corrected_featureMap.size(); ++i)
  {
    for(uint8_t j = 0; j <lactate_1_corrected_featureMap.at(i).getSubordinates().size(); ++j)
    {
      TEST_REAL_SIMILAR(lactate_1_corrected_featureMap.at(i).getSubordinates().at(j).getMetaValue("peak_apex_int"), L1_corrected[j]);
    }
  }

END_SECTION


START_SECTION(( void IsotopeLabelingMDVs::calculateIsotopicPurity(
                                              Feature& normalized_featuremap,
                                              const std::vector<double>& experiment_data,
                                              const std::string& isotopic_purity_name) ))

  // case 1: calculating isotopic purity on 1_2_13C, U_13C sample experiment data

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // L1_1_2_13C_glucose_experiment, L1_U_13C_glucose_experiment & L1_isotopic_purity_ground_truth
  // from "High-resolution 13C metabolic flux analysis",Long et al, doi:10.1038/s41596-019-0204-0,
  // P.2869, Box 4
  std::vector<double>           L1_norm_max                     {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<double>           L1_1_2_13C_glucose_experiment   {0.5, 0.7, 98.8, 0.0, 0.0, 0.0};
  std::vector<double>           L1_U_13C_glucose_experiment     {0.5, 0.0, 0.1, 0.2, 3.6, 95.5};
  std::vector<double>           L1_isotopic_purity_ground_truth {99.6469, 99.2517};  // [1_2_13C, U_13C]

  std::string                   L1_1_2_13C_glucose = "1_2-13C_glucose_experiment";
  std::string                   L1_U_13C_glucose = "U-13C_glucose_experiment";

  std::vector<OpenMS::Feature>  L1_subordinates_normmax;

  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  isotopelabelingmdvs.calculateIsotopicPurity(lactate_1_normalized, L1_1_2_13C_glucose_experiment, L1_1_2_13C_glucose);
  TEST_REAL_SIMILAR( (double)(lactate_1_normalized.getMetaValue(L1_1_2_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[0]);

  isotopelabelingmdvs.calculateIsotopicPurity(lactate_1_normalized, L1_U_13C_glucose_experiment, L1_U_13C_glucose);
  TEST_REAL_SIMILAR( (double)(lactate_1_normalized.getMetaValue(L1_U_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[1]);

END_SECTION


START_SECTION(( void IsotopeLabelingMDVs::calculateIsotopicPurities(
                                              Feature& normalized_featuremap,
                                              const std::vector<std::vector<double>>& experiment_data,
                                              const std::vector<std::string>& isotopic_purity_names) ))

  // case 1: calculating isotopic purity on 1_2_13C, U_13C sample experiment data

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  OpenMS::FeatureMap            lactate_1_featureMap;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // L1_1_2_13C_glucose_experiment, L1_U_13C_glucose_experiment & L1_isotopic_purity_ground_truth
  // from "High-resolution 13C metabolic flux analysis",Long et al, doi:10.1038/s41596-019-0204-0,
  // P.2869, Box 4
  std::vector<double>               L1_norm_max                     {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};
  std::vector<std::vector<double>>  L1_1_2_13C_glucose_experiment   {{0.5, 0.7, 98.8, 0.0, 0.0, 0.0},{0.5, 0.7, 98.8, 0.0, 0.0, 0.0},{0.5, 0.7, 98.8, 0.0, 0.0, 0.0}};
  std::vector<std::vector<double>>  L1_U_13C_glucose_experiment     {{0.5, 0.0, 0.1, 0.2, 3.6, 95.5},{0.5, 0.0, 0.1, 0.2, 3.6, 95.5},{0.5, 0.0, 0.1, 0.2, 3.6, 95.5}};
  std::vector<double>               L1_isotopic_purity_ground_truth {99.6469, 99.2517};  // [1_2_13C, U_13C]

  std::vector<std::string>          L1_1_2_13C_glucose = {"1_2-13C_glucose_experiment","1_2-13C_glucose_experiment","1_2-13C_glucose_experiment"};
  std::vector<std::string>          L1_U_13C_glucose = {"U-13C_glucose_experiment","U-13C_glucose_experiment","U-13C_glucose_experiment"};

  std::vector<OpenMS::Feature>  L1_subordinates_normmax;
  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  for(uint8_t i = 0; i < 3; ++i)
  {
    lactate_1_featureMap.push_back(lactate_1_normalized);
  }

  isotopelabelingmdvs.calculateIsotopicPurities(lactate_1_featureMap, L1_1_2_13C_glucose_experiment, L1_1_2_13C_glucose);
  for(uint8_t i = 0; i < lactate_1_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( (double)(lactate_1_featureMap.at(i).getMetaValue(L1_1_2_13C_glucose.at(i))) * 100, L1_isotopic_purity_ground_truth[0]);
  }

  isotopelabelingmdvs.calculateIsotopicPurities(lactate_1_featureMap, L1_U_13C_glucose_experiment, L1_U_13C_glucose);
  for(uint8_t i = 0; i < lactate_1_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( (double)(lactate_1_featureMap.at(i).getMetaValue(L1_U_13C_glucose.at(i))) * 100, L1_isotopic_purity_ground_truth[1]);
  }

END_SECTION


START_SECTION(( IsotopeLabelingMDVs::calculateMDVAccuracy(
                                              Feature& normalized_feature,
                                              const std::string& feature_name,
                                              const std::string& fragment_isotopomer_theoretical_formula) ))

  // case 1: calculating accuracy given measured values of 2 seperate features

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               feature_1, feature_2;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // accoa_C23H37N7O17P3S_MRM_measured_13 & fad_C27H32N9O15P2_EPI_measured_48 are extracted from
  // "MID Max: LC–MS/MS Method for Measuring the Precursor and Product Mass Isotopomer Distributions
  // of Metabolic Intermediates and Cofactors for Metabolic Flux Analysis Applications, McCloskey et al",
  // DOI: 10.1021/acs.analchem.5b03887, Supporting Information: Table S-2
  std::vector<double>           L1_norm_max                             {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};

  std::vector<double>           accoa_C23H37N7O17P3S_MRM_measured_13    {0.627, 0.253, 0.096, 0.02, 0.004, 0.001};
  std::vector<double>           accoa_C23H37N7O17P3S_abs_diff           {0.0632108, 0.0505238, 0.0119821, 0.0014131, 0.0000315, 0.0003232};
  std::vector<double>           fad_C27H32N9O15P2_EPI_measured_48       {0.638, 0.355, 0.1, 0.0, 0.0, 0.0};
  std::vector<double>           fad_C27H32N9O15P2_abs_diff              {0.0570446, 0.1223954, 0.0407946, 0.0111298, 0.0017729, 0.0002426};
  std::vector<double>           Average_accuracy_groundtruth            {0.02374, 0.03451}; // [accoa_13, fad_48]

  std::vector<OpenMS::Feature>  L1_subordinates, L2_subordinates;
  
  feature_1.setMetaValue("PeptideRef", "accoa");
  for (uint16_t i = 0; i < accoa_C23H37N7O17P3S_MRM_measured_13.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", accoa_C23H37N7O17P3S_MRM_measured_13[i]);
    L1_subordinates.push_back(sub);
  }
  feature_1.setSubordinates(L1_subordinates);

  feature_2.setMetaValue("PeptideRef", "fad");
  for (uint16_t i = 0; i < fad_C27H32N9O15P2_EPI_measured_48.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", fad_C27H32N9O15P2_EPI_measured_48[i]);
    L2_subordinates.push_back(sub);
  }
  feature_2.setSubordinates(L2_subordinates);

  isotopelabelingmdvs.calculateMDVAccuracy(feature_1, "peak_apex_int", "C23H37N7O17P3S");
  TEST_REAL_SIMILAR( feature_1.getMetaValue("average_accuracy"), Average_accuracy_groundtruth[0] );
  for (size_t feature_subordinate = 0; feature_subordinate < feature_1.getSubordinates().size(); ++feature_subordinate)
  {
    TEST_REAL_SIMILAR( feature_1.getSubordinates().at(feature_subordinate).getMetaValue("absolute_difference"), accoa_C23H37N7O17P3S_abs_diff[feature_subordinate] );
  }

  isotopelabelingmdvs.calculateMDVAccuracy(feature_2, "peak_apex_int", "C27H32N9O15P2");
  TEST_REAL_SIMILAR( feature_2.getMetaValue("average_accuracy"), Average_accuracy_groundtruth[1] );
  for (size_t feature_subordinate = 0; feature_subordinate < feature_2.getSubordinates().size(); ++feature_subordinate)
  {
    TEST_REAL_SIMILAR( feature_2.getSubordinates().at(feature_subordinate).getMetaValue("absolute_difference"), fad_C27H32N9O15P2_abs_diff[feature_subordinate] );
  }

END_SECTION

START_SECTION(( IsotopeLabelingMDVs::calculateMDVAccuracies(
                                              FeatureMap& normalized_featureMap,
                                              const std::string& feature_name,
                                              const std::map<std::string, std::string>& fragment_isotopomer_theoretical_formulas) ))

  // case 1: calculating accuracy given theoretical and measured values

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               feature_1, feature_2;
  OpenMS::FeatureMap            featureMap_1, featureMap_2;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // accoa_C23H37N7O17P3S_MRM_measured_13 & fad_C27H32N9O15P2_EPI_measured_48 are extracted from
  // "MID Max: LC–MS/MS Method for Measuring the Precursor and Product Mass Isotopomer Distributions
  // of Metabolic Intermediates and Cofactors for Metabolic Flux Analysis Applications, McCloskey et al",
  // DOI: 10.1021/acs.analchem.5b03887, Supporting Information: Table S-2
  std::vector<double>               L1_norm_max                             {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};

  std::vector<double>               accoa_C23H37N7O17P3S_MRM_measured_13    {0.627, 0.253, 0.096, 0.02, 0.004, 0.001};
  std::vector<double>               accoa_C23H37N7O17P3S_abs_diff           {0.0632108, 0.0505238, 0.0119821, 0.0014131, 0.0000315, 0.0003232};
  std::vector<double>               fad_C27H32N9O15P2_EPI_measured_48       {0.638, 0.355, 0.1, 0.0, 0.0, 0.0};
  std::vector<double>               fad_C27H32N9O15P2_abs_diff              {0.0570446, 0.1223954, 0.0407946, 0.0111298, 0.0017729, 0.0002426};
  std::vector<double>               Average_accuracy_groundtruth            {0.02374, 0.03451}; // [accoa_13, fad_48]

  std::map<std::string,std::string> theoretical_formulas                    {{"accoa","C23H37N7O17P3S"},{"fad","C27H32N9O15P2"}};

  std::vector<OpenMS::Feature>      L1_subordinates, L2_subordinates;
  
  feature_1.setMetaValue("PeptideRef", "accoa");
  for (uint16_t i = 0; i < accoa_C23H37N7O17P3S_MRM_measured_13.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", accoa_C23H37N7O17P3S_MRM_measured_13[i]);
    L1_subordinates.push_back(sub);
  }
  feature_1.setSubordinates(L1_subordinates);

  for (uint8_t i = 0; i < 3; ++i)
  {
    featureMap_1.push_back(feature_1);
  }

  feature_2.setMetaValue("PeptideRef", "fad");
  for (uint16_t i = 0; i < fad_C27H32N9O15P2_EPI_measured_48.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate2_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", fad_C27H32N9O15P2_EPI_measured_48[i]); // TODO:test for a missing val
    L2_subordinates.push_back(sub);
  }
  feature_2.setSubordinates(L2_subordinates);

  for (uint8_t i = 0; i < 3; ++i)
  {
    featureMap_2.push_back(feature_2);
  }

  isotopelabelingmdvs.calculateMDVAccuracies(featureMap_1, "peak_apex_int", theoretical_formulas);
  for (size_t i = 0; i < featureMap_1.size(); ++i)
  {
    TEST_REAL_SIMILAR(featureMap_1.at(i).getMetaValue("average_accuracy"), Average_accuracy_groundtruth[0]);
    
    for (size_t feature_subordinate = 0; feature_subordinate < featureMap_1.at(i).getSubordinates().size(); ++feature_subordinate)
    {
      TEST_REAL_SIMILAR(featureMap_1.at(i).getSubordinates().at(feature_subordinate).getMetaValue("absolute_difference"),
                        accoa_C23H37N7O17P3S_abs_diff[feature_subordinate]);
    }
  }

  isotopelabelingmdvs.calculateMDVAccuracies(featureMap_2, "peak_apex_int", theoretical_formulas);
  for (size_t i = 0; i < featureMap_2.size(); ++i)
  {
    TEST_REAL_SIMILAR(featureMap_2.at(i).getMetaValue("average_accuracy"), Average_accuracy_groundtruth[1]);
    
    for (size_t feature_subordinate = 0; feature_subordinate < featureMap_2.at(i).getSubordinates().size(); ++feature_subordinate)
    {
      TEST_REAL_SIMILAR(featureMap_2.at(i).getSubordinates().at(feature_subordinate).getMetaValue("absolute_difference"),
                        fad_C27H32N9O15P2_abs_diff[feature_subordinate]);
    }
  }

END_SECTION

END_TEST
