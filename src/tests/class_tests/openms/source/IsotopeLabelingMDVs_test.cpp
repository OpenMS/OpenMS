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
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <assert.h>

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
  
  isotopelabelingmdvs.calculateMDV(lactate_1_normmax, lactate_1_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_1_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normmax.getSubordinates().at(i).getIntensity(), L1_norm_max.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_2_normmax, lactate_2_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_2_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normmax.getSubordinates().at(i).getIntensity(), L2_norm_max.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_1_normsum, lactate_1_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_1_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normsum.getSubordinates().at(i).getIntensity(), L1_norm_sum.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_2_normsum, lactate_2_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

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
  
  isotopelabelingmdvs.calculateMDV(lactate_1_normmax, lactate_1_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_1_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normmax.getSubordinates().at(i).getIntensity(), L1_norm_max.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_2_normmax, lactate_2_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_2_normalized_normmax.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normmax.getSubordinates().at(i).getIntensity(), L2_norm_max.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_1_normsum, lactate_1_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_1_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_normalized_normsum.getSubordinates().at(i).getIntensity(), L1_norm_sum.at(i));
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

  isotopelabelingmdvs.calculateMDV(lactate_2_normsum, lactate_2_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_2_normalized_normsum.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_2_normalized_normsum.getSubordinates().at(i).getIntensity(), L2_norm_sum.at(i));
  }

  OpenMS::FeatureMap  lactate_normmax;
  OpenMS::FeatureMap  lactate_normsum;
  OpenMS::FeatureMap  lactate_normalized_normmax;
  OpenMS::FeatureMap  lactate_normalized_normsum;

  lactate_normmax.push_back(lactate_1_normmax);
  lactate_normmax.push_back(lactate_2_normmax);

  lactate_normsum.push_back(lactate_1_normsum);
  lactate_normsum.push_back(lactate_2_normsum);

  isotopelabelingmdvs.calculateMDVs(lactate_normmax, lactate_normalized_normmax, IsotopeLabelingMDVs::MassIntensityType::NORM_MAX, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);
  isotopelabelingmdvs.calculateMDVs(lactate_normsum, lactate_normalized_normsum, IsotopeLabelingMDVs::MassIntensityType::NORM_SUM, IsotopeLabelingMDVs::FeatureName::PEAK_APEX_INT);

  for(size_t i = 0; i < lactate_normalized_normmax.size(); ++i)
  {
    for(size_t j = 0; j < lactate_normalized_normmax.at(i).getSubordinates().size(); ++j)
    {
      if (i == 0) // lactate_1
      {
        TEST_REAL_SIMILAR(lactate_normalized_normmax.at(i).getSubordinates().at(j).getIntensity(), L1_norm_max.at(j));
      }
      else if (i == 1) // lactate_2
      {
        TEST_REAL_SIMILAR(lactate_normalized_normmax.at(i).getSubordinates().at(j).getIntensity(), L2_norm_max.at(j));
      }
    }
  }

  for(size_t i = 0; i < lactate_normalized_normsum.size(); ++i)
  {
    for(size_t j = 0; j < lactate_normalized_normsum.at(i).getSubordinates().size(); ++j)
    {
      if (i == 0) // lactate_1
      {
        TEST_REAL_SIMILAR(lactate_normalized_normsum.at(i).getSubordinates().at(j).getIntensity(), L1_norm_sum.at(j));
      }
      else if (i == 1) // lactate_2
      {
        TEST_REAL_SIMILAR(lactate_normalized_normsum.at(i).getSubordinates().at(j).getIntensity(), L2_norm_sum.at(j));
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
  correction_matrix_tBDMS.setMatrix<4,4>(correction_matrix_tBDMS_);

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

isotopelabelingmdvs.isotopicCorrection(lactate_1_normalized, lactate_1_corrected, correction_matrix_tBDMS, IsotopeLabelingMDVs::DerivatizationAgent::NOT_SELECTED);
  for(size_t i = 0; i < lactate_1_corrected.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_corrected.getSubordinates().at(i).getIntensity(), L1_corrected[i]);
  }

isotopelabelingmdvs.isotopicCorrection(lactate_1_normalized, lactate_1_corrected, {}, IsotopeLabelingMDVs::DerivatizationAgent::TBDMS);
  for(size_t i = 0; i < lactate_1_corrected.getSubordinates().size(); ++i)
  {
    TEST_REAL_SIMILAR(lactate_1_corrected.getSubordinates().at(i).getIntensity(), L1_corrected[i]);
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
  correction_matrix_tBDMS.setMatrix<4,4>(correction_matrix_tBDMS_);

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
      TEST_REAL_SIMILAR(lactate_1_corrected_featureMap.at(i).getSubordinates().at(j).getIntensity(), L1_corrected[j]);
    }
  }

  lactate_1_corrected_featureMap.clear();

  isotopelabelingmdvs.isotopicCorrections(lactate_1_featureMap, lactate_1_corrected_featureMap, {}, IsotopeLabelingMDVs::DerivatizationAgent::TBDMS);
  for(uint8_t i = 0; i < lactate_1_corrected_featureMap.size(); ++i)
  {
//    for(uint8_t j = 0; j < lactate_1_corrected.getSubordinates().size(); ++j)
    {
//      TEST_REAL_SIMILAR(lactate_1_corrected_featureMap.at(i).getSubordinates().at(j).getIntensity(), L1_corrected[j]);
    }
  }

END_SECTION


START_SECTION(( void IsotopeLabelingMDVs::calculateIsotopicPurity(
                                              const Feature& normalized_featuremap,
                                              Feature& featuremap_with_isotopic_purity,
                                              const std::vector<double>& experiment_data,
                                              const std::string& isotopic_purity_name) ))

  // case 1: calculating isotopic purity on 1_2_13C, U_13C sample experiment data

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  OpenMS::Feature               lactate_1_with_isotopic_purity;
  
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

  isotopelabelingmdvs.calculateIsotopicPurity(lactate_1_normalized, lactate_1_with_isotopic_purity, L1_1_2_13C_glucose_experiment, L1_1_2_13C_glucose);
  TEST_REAL_SIMILAR( (double)(lactate_1_with_isotopic_purity.getMetaValue(L1_1_2_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[0]);

  isotopelabelingmdvs.calculateIsotopicPurity(lactate_1_normalized, lactate_1_with_isotopic_purity, L1_U_13C_glucose_experiment, L1_U_13C_glucose);
  TEST_REAL_SIMILAR( (double)(lactate_1_with_isotopic_purity.getMetaValue(L1_U_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[1]);

END_SECTION


START_SECTION(( void IsotopeLabelingMDVs::calculateIsotopicPurities(
                                              const Feature& normalized_featuremap,
                                              const Feature& featuremap_with_isotopic_purity,
                                              std::vector<double>& experiment_data,
                                              const std::string& isotopic_purity_name) ))

  // case 1: calculating isotopic purity on 1_2_13C, U_13C sample experiment data

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  OpenMS::Feature               lactate_1_with_isotopic_purity;
  OpenMS::FeatureMap            lactate_1_featureMap;
  OpenMS::FeatureMap            lactate_1_with_isotopic_purity_featureMap;
  
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

  for(uint8_t i = 0; i < 3; ++i)
  {
    lactate_1_featureMap.push_back(lactate_1_normalized);
  }

  isotopelabelingmdvs.calculateIsotopicPurities(lactate_1_featureMap, lactate_1_with_isotopic_purity_featureMap, L1_1_2_13C_glucose_experiment, L1_1_2_13C_glucose);
  for(uint8_t i = 0; i < lactate_1_with_isotopic_purity_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( (double)(lactate_1_with_isotopic_purity_featureMap.at(i).getMetaValue(L1_1_2_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[0]);
  }

  lactate_1_with_isotopic_purity_featureMap.clear();
  isotopelabelingmdvs.calculateIsotopicPurities(lactate_1_featureMap, lactate_1_with_isotopic_purity_featureMap, L1_U_13C_glucose_experiment, L1_U_13C_glucose);
  for(uint8_t i = 0; i < lactate_1_with_isotopic_purity_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( (double)(lactate_1_with_isotopic_purity_featureMap.at(i).getMetaValue(L1_U_13C_glucose)) * 100, L1_isotopic_purity_ground_truth[1]);
  }

END_SECTION


START_SECTION(( IsotopeLabelingMDVs::calculateMDVAccuracy(
                                              const Feature& normalized_feature,
                                              Feature& feature_with_accuracy_info,
                                              const std::vector<double>& fragment_isotopomer_measured,
                                              const std::string& fragment_isotopomer_theoretical_formula) ))

  // case 1: calculating accuracy given theoretical and measured values

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  OpenMS::Feature               lactate_1_with_accuracy_info;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // accoa_C23H37N7O17P3S_MRM_measured_13 & fad_C27H32N9O15P2_EPI_measured_48 are extracted from
  // "MID Max: LC–MS/MS Method for Measuring the Precursor and Product Mass Isotopomer Distributions
  // of Metabolic Intermediates and Cofactors for Metabolic Flux Analysis Applications, McCloskey et al",
  // DOI: 10.1021/acs.analchem.5b03887, Supporting Information: Table S-2
  std::vector<double>           L1_norm_max                             {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};

  std::vector<double>           accoa_C23H37N7O17P3S_MRM_theoretical_13 ;
  std::vector<double>           fad_C27H32N9O15P2_EPI_theoretical_48    ;

  std::vector<double>           accoa_C23H37N7O17P3S_MRM_measured_13    {0.627, 0.253, 0.096, 0.02, 0.004, 0.001};
  std::vector<double>           fad_C27H32N9O15P2_EPI_measured_48       {0.638, 0.355, 0.1, 0.0, 0.0, 0.0};
  std::vector<double>           Average_accuracy_groundtruth            {0.02374, 0.03451}; // [accoa_13, fad_48]

  std::vector<OpenMS::Feature>  L1_subordinates_normmax;

  IsotopeDistribution accoa_C23H37N7O17P3S_MRM_theoretical_iso(EmpiricalFormula("C23H37N7O17P3S").getIsotopeDistribution(CoarseIsotopePatternGenerator(6)));
  for (IsotopeDistribution::ConstIterator it = accoa_C23H37N7O17P3S_MRM_theoretical_iso.begin(); it != accoa_C23H37N7O17P3S_MRM_theoretical_iso.end(); ++it)
  {
    accoa_C23H37N7O17P3S_MRM_theoretical_13.push_back( it->getIntensity() );
  }

  IsotopeDistribution fad_C27H32N9O15P2_EPI_theoretical_iso(EmpiricalFormula("C27H32N9O15P2").getIsotopeDistribution(CoarseIsotopePatternGenerator(6)));
  for (IsotopeDistribution::ConstIterator it = fad_C27H32N9O15P2_EPI_theoretical_iso.begin(); it != fad_C27H32N9O15P2_EPI_theoretical_iso.end(); ++it)
  {
    fad_C27H32N9O15P2_EPI_theoretical_48.push_back( it->getIntensity() );
  }
  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  isotopelabelingmdvs.calculateMDVAccuracy(lactate_1_normalized, lactate_1_with_accuracy_info, accoa_C23H37N7O17P3S_MRM_measured_13, "C23H37N7O17P3S");
  TEST_REAL_SIMILAR( lactate_1_with_accuracy_info.getMetaValue("average_accuracy"), Average_accuracy_groundtruth[0] );
  lactate_1_with_accuracy_info.clearMetaInfo();

  isotopelabelingmdvs.calculateMDVAccuracy(lactate_1_normalized, lactate_1_with_accuracy_info, fad_C27H32N9O15P2_EPI_measured_48, "C27H32N9O15P2");
  TEST_REAL_SIMILAR( lactate_1_with_accuracy_info.getMetaValue("average_accuracy"), Average_accuracy_groundtruth[1] );

END_SECTION

START_SECTION(( IsotopeLabelingMDVs::calculateMDVAccuracies(
                                              const FeatureMap& normalized_featureMap,
                                              FeatureMap& featureMap_with_accuracy_info,
                                              const std::vector<double>& fragment_isotopomer_measured,
                                              const std::string& fragment_isotopomer_theoretical_formula) ))

  // case 1: calculating accuracy given theoretical and measured values

  IsotopeLabelingMDVs           isotopelabelingmdvs;
  OpenMS::Feature               lactate_1_normalized;
  OpenMS::Feature               lactate_1_with_accuracy_info;
  OpenMS::FeatureMap            lactate_1_featureMap;
  OpenMS::FeatureMap            lactate_1_with_accuracy_info_featureMap;
  
  // L1_norm_max From CHO_190316_Flux.xlsx provided by Douglas McCloskey
  // accoa_C23H37N7O17P3S_MRM_measured_13 & fad_C27H32N9O15P2_EPI_measured_48 are extracted from
  // "MID Max: LC–MS/MS Method for Measuring the Precursor and Product Mass Isotopomer Distributions
  // of Metabolic Intermediates and Cofactors for Metabolic Flux Analysis Applications, McCloskey et al",
  // DOI: 10.1021/acs.analchem.5b03887, Supporting Information: Table S-2
  std::vector<double>           L1_norm_max                             {1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05};

  std::vector<double>           accoa_C23H37N7O17P3S_MRM_theoretical_13 ;
  std::vector<double>           fad_C27H32N9O15P2_EPI_theoretical_48    ;

  std::vector<double>           accoa_C23H37N7O17P3S_MRM_measured_13    {0.627, 0.253, 0.096, 0.02, 0.004, 0.001};
  std::vector<double>           fad_C27H32N9O15P2_EPI_measured_48       {0.638, 0.355, 0.1, 0.0, 0.0, 0.0};
  std::vector<double>           Average_accuracy_groundtruth            {0.02374, 0.03451}; // [accoa_13, fad_48]

  std::vector<OpenMS::Feature>  L1_subordinates_normmax;

  IsotopeDistribution accoa_C23H37N7O17P3S_MRM_theoretical_iso(EmpiricalFormula("C23H37N7O17P3S").getIsotopeDistribution(CoarseIsotopePatternGenerator(6)));
  for (IsotopeDistribution::ConstIterator it = accoa_C23H37N7O17P3S_MRM_theoretical_iso.begin(); it != accoa_C23H37N7O17P3S_MRM_theoretical_iso.end(); ++it)
  {
    accoa_C23H37N7O17P3S_MRM_theoretical_13.push_back( it->getIntensity() );
  }

  IsotopeDistribution fad_C27H32N9O15P2_EPI_theoretical_iso(EmpiricalFormula("C27H32N9O15P2").getIsotopeDistribution(CoarseIsotopePatternGenerator(6)));
  for (IsotopeDistribution::ConstIterator it = fad_C27H32N9O15P2_EPI_theoretical_iso.begin(); it != fad_C27H32N9O15P2_EPI_theoretical_iso.end(); ++it)
  {
    fad_C27H32N9O15P2_EPI_theoretical_48.push_back( it->getIntensity() );
  }
  
  lactate_1_normalized.setMetaValue("PeptideRef", "Lactate1");
  for (uint16_t i = 0; i < L1_norm_max.size(); ++i)
  {
    OpenMS::Feature sub;
    sub.setMetaValue("native_id", "Lactate1_"+std::to_string(117+i));
    sub.setMetaValue("peak_apex_int", L1_norm_max[i]);
    L1_subordinates_normmax.push_back(sub);
  }
  lactate_1_normalized.setSubordinates(L1_subordinates_normmax);

  for (uint8_t i = 0; i < 3; ++i)
  {
    lactate_1_featureMap.push_back(lactate_1_normalized);
  }

  isotopelabelingmdvs.calculateMDVAccuracies(lactate_1_featureMap, lactate_1_with_accuracy_info_featureMap, accoa_C23H37N7O17P3S_MRM_measured_13, "C23H37N7O17P3S");
  for (size_t i = 0; i < lactate_1_with_accuracy_info_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( lactate_1_with_accuracy_info_featureMap.at(i).getMetaValue("average_accuracy"), Average_accuracy_groundtruth[0] );
  }
  lactate_1_with_accuracy_info.clearMetaInfo();
  lactate_1_with_accuracy_info_featureMap.clear();

  isotopelabelingmdvs.calculateMDVAccuracies(lactate_1_featureMap, lactate_1_with_accuracy_info_featureMap, fad_C27H32N9O15P2_EPI_measured_48, "C27H32N9O15P2");
  for (size_t i = 0; i < lactate_1_with_accuracy_info_featureMap.size(); ++i)
  {
    TEST_REAL_SIMILAR( lactate_1_with_accuracy_info_featureMap.at(i).getMetaValue("average_accuracy"), Average_accuracy_groundtruth[1] );
  }

END_SECTION

END_TEST
