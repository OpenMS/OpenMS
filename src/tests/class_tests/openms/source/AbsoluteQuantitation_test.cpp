// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>

using namespace OpenMS;
using namespace std;
///////////////////////////

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_serL_standards()
{
  // TEST 1: ser-L
  static const double arrx1[] = {2.32e4,2.45e4,1.78e4,2.11e4,1.91e4,
    2.06e4,1.85e4,1.53e4,1.40e4,1.03e4,1.07e4,6.68e3,5.27e3,2.83e3};
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = {4.94e3,6.55e3,7.37e3,1.54e4,2.87e4,
    5.41e4,1.16e5,1.85e5,3.41e5,7.54e5,9.76e5,1.42e6,1.93e6,2.23e6};
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );
  static const double arrz1[] = {1.00e-2,2.00e-2,4.00e-2,1.00e-1,2.00e-1,
    4.00e-1,1.00e0,2.00e0,4.00e0,1.00e1,2.00e1,4.00e1,1.00e2,2.00e2};
  std::vector<double> z1 (arrz1, arrz1 + sizeof(arrz1) / sizeof(arrz1[0]) );

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x1.size(); ++i)
  {
    component.setMetaValue("native_id","ser-L.ser-L_1.Light");
    component.setMetaValue("peak_apex_int",y1[i]);
    IS_component.setMetaValue("native_id","ser-L.ser-L_1.Heavy");
    IS_component.setMetaValue("peak_apex_int",x1[i]);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = z1[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }
  return component_concentrations;
}

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_amp_standards()
{
  // TEST 2: amp
  static const double arrx2[] = {2.15e5,2.32e5,2.69e5,2.53e5,2.50e5,
  2.75e5,2.67e5,3.31e5,3.15e5,3.04e5,3.45e5,3.91e5,4.62e5,3.18e5};
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = {4.40e2,1.15e3,1.53e3,2.01e3,4.47e3,
  7.36e3,2.18e4,4.46e4,8.50e4,2.33e5,5.04e5,1.09e6,2.54e6,3.64e6};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );
  static const double arrz2[] = {2.00e-3,4.00e-3,8.00e-3,2.00e-2,4.00e-2,
  8.00e-2,2.00e-1,4.00e-1,8.00e-1,2.00e0,4.00e0,8.00e0,2.00e1,4.00e1};
  std::vector<double> z2 (arrz2, arrz2 + sizeof(arrz2) / sizeof(arrz2[0]) );

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x2.size(); ++i)
  {
    component.setMetaValue("native_id","amp.amp_1.Light");
    component.setMetaValue("peak_apex_int",y2[i]);
    IS_component.setMetaValue("native_id","amp.amp_1.Heavy");
    IS_component.setMetaValue("peak_apex_int",x2[i]);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = z2[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }
  return component_concentrations;
}

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_atp_standards()
{
  // TEST 3: atp
  static const double arrx3[] = {8.28e2,1.32e3,1.57e3,1.63e3,1.48e3,
  2.43e3,4.44e3,1.03e4,1.75e4,6.92e4,1.97e5,2.69e5,3.20e5,3.22e5};
  std::vector<double> x3 (arrx3, arrx3 + sizeof(arrx3) / sizeof(arrx3[0]) );
  static const double arry3[] = {2.21e2,4.41e2,3.31e2,2.21e2,3.09e2,
  5.96e2,1.26e3,2.49e3,1.12e4,8.79e4,4.68e5,1.38e6,3.46e6,4.19e6};
  std::vector<double> y3 (arry3, arry3 + sizeof(arry3) / sizeof(arry3[0]) );
  static const double arrz3[] = {2.00e-3,4.00e-3,8.00e-3,2.00e-2,4.00e-2,
  8.00e-2,2.00e-1,4.00e-1,8.00e-1,2.00e0,4.00e0,8.00e0,2.00e1,4.00e1};
  std::vector<double> z3 (arrz3, arrz3 + sizeof(arrz3) / sizeof(arrz3[0]) );

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x3.size(); ++i)
  {
    component.setMetaValue("native_id","atp.atp_1.Light");
    component.setMetaValue("peak_apex_int",y3[i]);
    IS_component.setMetaValue("native_id","atp.atp_1.Heavy");
    IS_component.setMetaValue("peak_apex_int",x3[i]);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = z3[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  return component_concentrations;
}

/////////////////////////////////////////////////////////////

START_TEST(AbsoluteQuantitation, "$Id$")

/////////////////////////////////////////////////////////////

class AbsoluteQuantitation_test : public AbsoluteQuantitation
{
  public :

    int jackknifeOutlierCandidate_(const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params)
    {
      return AbsoluteQuantitation::jackknifeOutlierCandidate_(component_concentrations,
        feature_name,
        transformation_model,
        transformation_model_params);
    }

    int residualOutlierCandidate_(const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params)
    {
      return AbsoluteQuantitation::residualOutlierCandidate_(component_concentrations,
        feature_name,
        transformation_model,
        transformation_model_params);
    }

    std::vector<AbsoluteQuantitationStandards::featureConcentration> extractComponents_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      std::vector<size_t> component_concentrations_indices)
    {
      return AbsoluteQuantitation::extractComponents_(
        component_concentrations,
        component_concentrations_indices);
    }

};

/////////////////////////////////////////////////////////////

AbsoluteQuantitation* ptr = nullptr;
AbsoluteQuantitation* nullPointer = nullptr;


START_SECTION((AbsoluteQuantitation()))
  ptr = new AbsoluteQuantitation();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~AbsoluteQuantitation()))
  delete ptr;
END_SECTION

START_SECTION((double calculateRatio(const Feature & component_1, const Feature & component_2, const String feature_name)))
  AbsoluteQuantitation absquant;
  String feature_name = "peak_apex_int";
  double inf = std::numeric_limits<double>::infinity();
  // dummy features
  OpenMS::Feature component_1, component_2;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id","component1");
  component_2.setMetaValue(feature_name, 5.0);
  component_2.setMetaValue("native_id","component2");
  // tests
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_1,component_2,feature_name),1.0);
  component_2.setMetaValue(feature_name, 0.0);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_1,component_2,feature_name),inf);
  // dummy features
  OpenMS::Feature component_3, component_4;
  component_3.setMetaValue("peak_area", 5.0);
  component_3.setMetaValue("native_id","component3");
  component_4.setMetaValue("peak_area", 5.0);
  component_4.setMetaValue("native_id","component4");
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_1,component_4,feature_name),5.0);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_3,component_4,feature_name),0.0);
  // feature_name == "intensity"
  Feature component_5, component_6, component_7, component_8;
  feature_name = "intensity";
  component_5.setMetaValue("native_id", "component5");
  component_6.setMetaValue("native_id", "component6");
  component_5.setIntensity(3.0);
  component_6.setIntensity(4.0);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_5, component_6, feature_name), 0.75);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_6, component_5, feature_name), 1.33333333333333);
  component_7.setMetaValue("native_id", "component7");
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_5, component_7, feature_name), inf);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_5, component_8, feature_name), 3.0);
END_SECTION

START_SECTION((double calculateBias(const double & actual_concentration, const double & calculated_concentration)))
  AbsoluteQuantitation absquant;
  double actual_concentration = 5.0;
  double calculated_concentration = 5.0;
  TEST_REAL_SIMILAR(absquant.calculateBias(actual_concentration,calculated_concentration),0.0);
  calculated_concentration = 4.0;
  TEST_REAL_SIMILAR(absquant.calculateBias(actual_concentration,calculated_concentration),20.0);
END_SECTION

START_SECTION((double applyCalibration(const Feature & component,
  const Feature & IS_component,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params)))

  AbsoluteQuantitation absquant;

  // set-up the features
  Feature component, IS_component;
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",2.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  String feature_name = "peak_apex_int";

  // set-up the model and params
  // y = m*x + b
  // x = (y - b)/m
  String transformation_model;
  Param param;
  // transformation_model = "TransformationModelLinear";
  transformation_model = "linear";
  param.setValue("slope",2.0);
  param.setValue("intercept",1.0);

  TEST_REAL_SIMILAR(absquant.applyCalibration(component,
    IS_component,
    feature_name,
    transformation_model,
    param),0.5);
END_SECTION

START_SECTION((void quantifyComponents(std::vector<FeatureMap>& unknowns)))

  AbsoluteQuantitation absquant;

  // set-up the unknown FeatureMap
  FeatureMap unknown_feature_map;
  // set-up the features and sub-features
  std::vector<Feature> unknown_feature_subordinates;
  Feature unknown_feature, component, IS_component;
  String feature_name = "peak_apex_int";
  // component 1
  unknown_feature.setMetaValue("PeptideRef","component_group1");
  component.setMetaValue("native_id","component1");
  component.setMetaValue(feature_name,2.0);
  IS_component.setMetaValue("native_id","IS1");
  IS_component.setMetaValue(feature_name,2.0);
  unknown_feature_subordinates.push_back(IS_component);
  unknown_feature_subordinates.push_back(component);
  unknown_feature.setSubordinates(unknown_feature_subordinates);
  unknown_feature_map.push_back(unknown_feature);
  unknown_feature_subordinates.clear();
  // component 2
  unknown_feature.setMetaValue("PeptideRef","component_group2");
  component.setMetaValue("native_id","component2");
  component.setMetaValue(feature_name,4.0);
  IS_component.setMetaValue("native_id","IS2");
  IS_component.setMetaValue(feature_name,4.0);
  unknown_feature_subordinates.push_back(IS_component);
  unknown_feature_subordinates.push_back(component);
  unknown_feature.setSubordinates(unknown_feature_subordinates);
  unknown_feature_map.push_back(unknown_feature);
  unknown_feature_subordinates.clear();
  // component 3
  unknown_feature.setMetaValue("PeptideRef","component_group3");
  component.setMetaValue("native_id","component3");
  component.setMetaValue(feature_name,6.0);
  IS_component.setMetaValue("native_id","IS3");
  IS_component.setMetaValue(feature_name,6.0);
  unknown_feature_subordinates.push_back(component);  // test order change
  unknown_feature_subordinates.push_back(IS_component);
  unknown_feature.setSubordinates(unknown_feature_subordinates);
  unknown_feature_map.push_back(unknown_feature);
  unknown_feature_subordinates.clear();
  // // set-up the unknowns
  // std::vector<FeatureMap> unknowns;
  // unknowns.push_back(unknown_feature_map);

  // set-up the model and params
  AbsoluteQuantitationMethod aqm;
  String transformation_model;
  Param param;
  // transformation_model = "TransformationModelLinear";
  transformation_model = "linear";
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  aqm.setTransformationModel(transformation_model);
  aqm.setTransformationModelParams(param);
  // set-up the quant_method map
  std::vector<AbsoluteQuantitationMethod> quant_methods;
  // component_1
  aqm.setComponentName("component1");
  aqm.setISName("IS1");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_2
  aqm.setComponentName("component2");
  aqm.setISName("IS1");
  aqm.setFeatureName(feature_name); // test IS outside component_group
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_3
  aqm.setComponentName("component3");
  aqm.setISName("IS3");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);

  absquant.setQuantMethods(quant_methods);
  absquant.quantifyComponents(unknown_feature_map);

  // DEBUGGING:
  // for (size_t i = 0; i < unknowns.size(); ++i)
  // {
  //   for (size_t j = 0; j < unknowns[i].size(); ++j)
  //   {
  //     for (size_t k = 0; k < unknowns[i][j].getSubordinates().size(); ++k)
  //     {
  //       std::cout << "component = " << unknowns[i][j].getSubordinates()[k].getMetaValue("native_id") << std::endl;
  //       std::cout << "calculated_concentration = " << unknowns[i][j].getSubordinates()[k].getMetaValue("calculated_concentration") << std::endl;
  //     }
  //   }
  // }

  TEST_EQUAL(unknown_feature_map[0].getSubordinates()[0].getMetaValue("calculated_concentration"),"");
  TEST_STRING_EQUAL(unknown_feature_map[0].getSubordinates()[0].getMetaValue("concentration_units"),"");
  TEST_REAL_SIMILAR(unknown_feature_map[0].getSubordinates()[1].getMetaValue("calculated_concentration"),1.0);
  TEST_STRING_EQUAL(unknown_feature_map[0].getSubordinates()[1].getMetaValue("concentration_units"),"uM");
  TEST_REAL_SIMILAR(unknown_feature_map[1].getSubordinates()[1].getMetaValue("calculated_concentration"),2.0);
  TEST_STRING_EQUAL(unknown_feature_map[1].getSubordinates()[1].getMetaValue("concentration_units"),"uM");
  TEST_REAL_SIMILAR(unknown_feature_map[2].getSubordinates()[0].getMetaValue("calculated_concentration"),1.0);
  TEST_STRING_EQUAL(unknown_feature_map[2].getSubordinates()[0].getMetaValue("concentration_units"),"uM");
END_SECTION

START_SECTION((void (
  const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params,
  std::vector<double> & biases,
  double & correlation_coefficient)))

  AbsoluteQuantitation absquant;

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  // point #1
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",1.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  component_concentration.feature = component;
  component_concentration.IS_feature = IS_component;
  component_concentration.actual_concentration = 2.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentration.dilution_factor = 2.0;
  component_concentrations.push_back(component_concentration);
  // point #2
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",2.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  component_concentration.feature = component;
  component_concentration.IS_feature = IS_component;
  component_concentration.actual_concentration = 4.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentration.dilution_factor = 2.0;
  component_concentrations.push_back(component_concentration);
  // point #3
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",3.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  component_concentration.feature = component;
  component_concentration.IS_feature = IS_component;
  component_concentration.actual_concentration = 6.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentration.dilution_factor = 2.0;
  component_concentrations.push_back(component_concentration);

  String feature_name = "peak_apex_int";

  // set-up the model and params
  // y = m*x + b
  // x = (y - b)/m
  String transformation_model;
  Param param;
  // transformation_model = "TransformationModelLinear";
  transformation_model = "linear";
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  std::vector<double> biases;
  double correlation_coefficient;

  absquant.calculateBiasAndR(
    component_concentrations,
    feature_name,
    transformation_model,
    param,
    biases,
    correlation_coefficient);

  TEST_REAL_SIMILAR(biases[0],0.0);
  TEST_REAL_SIMILAR(correlation_coefficient,1.0);

END_SECTION

START_SECTION((Param AbsoluteQuantitation::fitCalibration(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params)))

  AbsoluteQuantitation absquant;

  // TEST 1:
  static const double arrx1[] = {-1, -2, -3, 1, 2, 3};
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = {1, 1, 1, 1, 1, 1};
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );
  static const double arrz1[] = {-4, -8, -12, 4, 8, 12};
  std::vector<double> z1 (arrz1, arrz1 + sizeof(arrz1) / sizeof(arrz1[0]) );

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x1.size(); ++i)
  {
    component.setMetaValue("native_id","ser-L.ser-L_1.Light");
    component.setMetaValue("peak_apex_int",x1[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",y1[i]);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = z1[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 2.0;
    component_concentrations.push_back(component_concentration);
  }

  String feature_name = "peak_apex_int";
  Param transformation_model_params;
  transformation_model_params.setValue("x_datum_min", -1e12);
  transformation_model_params.setValue("x_datum_max", 1e12);
  transformation_model_params.setValue("y_datum_min", -1e12);
  transformation_model_params.setValue("y_datum_max", 1e12);
  // String transformation_model = "TransformationModelLinear";
  String transformation_model = "linear";

  Param param = absquant.fitCalibration(component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);

  TEST_REAL_SIMILAR(param.getValue("slope"),0.5);
  TEST_REAL_SIMILAR(param.getValue("intercept"),0.0);

  // TEST 2:
  static const double arrx2[] = {0.25,0.5,1,2,3,4,5,6};
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = {1,1,1,1,1,1,1,1};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );
  static const double arrz2[] = {0.5,1,2,4,6,8,10,12};
  std::vector<double> z2 (arrz2, arrz2 + sizeof(arrz2) / sizeof(arrz2[0]) );

  // set-up the features
  component_concentrations.clear();
  for (size_t i = 0; i < x2.size(); ++i)
  {
    component.setMetaValue("native_id","ser-L.ser-L_1.Light");
    component.setMetaValue("peak_apex_int",x2[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",y2[i]);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = z2[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  transformation_model_params.setValue("x_weight", "ln(x)");
  transformation_model_params.setValue("y_weight", "ln(y)");

  param = absquant.fitCalibration(component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);

  TEST_REAL_SIMILAR(param.getValue("slope"), 1.0);
  TEST_REAL_SIMILAR(param.getValue("intercept"), -0.69314718);

END_SECTION

START_SECTION((bool optimizeCalibrationCurveIterative(
  std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params,
  Param & optimized_params)))

  AbsoluteQuantitation absquant;

  // TEST 1: ser-L
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations = make_serL_standards();

  // set-up the class parameters
  Param absquant_params;
  absquant_params.setValue("min_points", 4);
  absquant_params.setValue("max_bias", 30.0);
  absquant_params.setValue("min_correlation_coefficient", 0.9);
  absquant_params.setValue("max_iters", 100);
  absquant_params.setValue("outlier_detection_method", "iter_jackknife");
  absquant_params.setValue("use_chauvenet", "false");
  absquant.setParameters(absquant_params);

  // set-up the function parameters
  const String feature_name = "peak_apex_int";
  const String transformation_model = "linear";
  Param transformation_model_params;
  transformation_model_params.setValue("x_weight", "ln(x)");
  transformation_model_params.setValue("y_weight", "ln(y)");
  transformation_model_params.setValue("slope", "1.0");
  transformation_model_params.setValue("intercept", "0.0");
  transformation_model_params.setValue("x_datum_min", -1e12);
  transformation_model_params.setValue("x_datum_max", 1e12);
  transformation_model_params.setValue("y_datum_min", -1e12);
  transformation_model_params.setValue("y_datum_max", 1e12);
  Param optimized_params;

  bool optimal_fit_found = absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.04);
  TEST_REAL_SIMILAR(component_concentrations[8].actual_concentration, 40.0);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), 0.9011392589);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), 1.870185076);
  TEST_EQUAL(optimal_fit_found, true);

  // TEST 2: amp
  component_concentrations = make_amp_standards();

  optimal_fit_found = absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(component_concentrations[8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), 0.95799683);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), -1.047543387);
  TEST_EQUAL(optimal_fit_found, true);

  // TEST 3: atp
  component_concentrations = make_atp_standards();

  optimal_fit_found = absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(component_concentrations[3].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), 0.623040824);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), 0.36130172586);
  TEST_EQUAL(optimal_fit_found, true);

  // TEST 4: atp with too stringent of a criteria (should fail)
  component_concentrations = make_atp_standards();
  
  absquant_params.setValue("min_points", 12);
  absquant_params.setValue("max_bias", 5.0);
  absquant_params.setValue("min_correlation_coefficient", 0.99);
  absquant_params.setValue("max_iters", 100);
  absquant_params.setValue("outlier_detection_method", "iter_jackknife");
  absquant_params.setValue("use_chauvenet", "false");
  absquant.setParameters(absquant_params);

  optimal_fit_found = absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), 0.482400123016735);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), 0.305125368008987);
  TEST_EQUAL(optimal_fit_found, false);

END_SECTION

START_SECTION((void optimizeCalibrationCurves(AbsoluteQuantitationStandards::components_to_concentrations & components_concentrations)))

  AbsoluteQuantitation absquant;

  // set-up the class parameters
  Param absquant_params;
  absquant_params.setValue("min_points", 4);
  absquant_params.setValue("max_bias", 30.0);
  absquant_params.setValue("min_correlation_coefficient", 0.9);
  absquant_params.setValue("max_iters", 100);
  absquant_params.setValue("outlier_detection_method", "iter_jackknife");
  absquant_params.setValue("use_chauvenet", "false");
  absquant.setParameters(absquant_params);

  // set up the quantitation method
  AbsoluteQuantitationMethod aqm;
  String feature_name = "peak_apex_int";
  String transformation_model;
  Param param;
  transformation_model = "linear";
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  param.setValue("x_weight", "ln(x)");
  param.setValue("y_weight", "ln(y)");
  param.setValue("x_datum_min", -1e12);
  param.setValue("x_datum_max", 1e12);
  param.setValue("y_datum_min", -1e12);
  param.setValue("y_datum_max", 1e12);
  aqm.setTransformationModel(transformation_model);
  aqm.setTransformationModelParams(param);
  // set-up the quant_method map
  std::vector<AbsoluteQuantitationMethod> quant_methods;
  // component_1
  aqm.setComponentName("ser-L.ser-L_1.Light");
  aqm.setISName("ser-L.ser-L_1.Heavy");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_2
  aqm.setComponentName("amp.amp_1.Light");
  aqm.setISName("amp.amp_1.Heavy");
  aqm.setFeatureName(feature_name); // test IS outside component_group
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_3
  aqm.setComponentName("atp.atp_1.Light");
  aqm.setISName("atp.atp_1.Heavy");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);

  absquant.setQuantMethods(quant_methods);

  // set up the standards
  std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> components_concentrations;
  components_concentrations["ser-L.ser-L_1.Light"] = make_serL_standards();
  components_concentrations["amp.amp_1.Light"] = make_amp_standards();
  components_concentrations["atp.atp_1.Light"] = make_atp_standards();

  absquant.optimizeCalibrationCurves(components_concentrations);
  std::map<String, AbsoluteQuantitationMethod> quant_methods_map = absquant.getQuantMethodsAsMap();

  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][0].actual_concentration, 0.04);
  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][8].actual_concentration, 40.0);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("slope"), 0.9011392589);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("intercept"), 1.87018507);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getCorrelationCoefficient(), 0.999320072);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getLLOQ(), 0.04);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getULOQ(), 200);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getNPoints(), 11);

  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("slope"), 0.95799683);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("intercept"), -1.047543387);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getCorrelationCoefficient(), 0.99916926);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getLLOQ(), 0.02);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getULOQ(), 40.0);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getNPoints(), 11);

  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][3].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("slope"), 0.623040824);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("intercept"), 0.36130172586);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getCorrelationCoefficient(), 0.998208402);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getLLOQ(), 0.02);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getULOQ(), 40.0);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getNPoints(), 6);

  components_concentrations.clear();
  components_concentrations["ser-L.ser-L_1.Light"] = make_serL_standards();
  absquant_params.setValue("min_points", 20); // so that an optimal calibration curve is not found, and the sorting is not saved
  absquant.setParameters(absquant_params);
  absquant.optimizeCalibrationCurves(components_concentrations);
  quant_methods_map = absquant.getQuantMethodsAsMap();

  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][0].actual_concentration, 0.01);
  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][8].actual_concentration, 4.0);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("slope"), 0.9011392589);  // previous supplied
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("intercept"), 1.87018507);  // previous supplied
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getCorrelationCoefficient(), 0.0);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getLLOQ(), 0.0);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getULOQ(), 0.0);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getNPoints(), 0.0);
END_SECTION

START_SECTION(void optimizeSingleCalibrationCurve(
  const String& component_name,
  std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations
))
  // set up the quantitation method
  Param param;
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  param.setValue("x_weight", "ln(x)");
  param.setValue("y_weight", "ln(y)");
  param.setValue("x_datum_min", -1e12);
  param.setValue("x_datum_max", 1e12);
  param.setValue("y_datum_min", -1e12);
  param.setValue("y_datum_max", 1e12);
  AbsoluteQuantitationMethod aqm;
  aqm.setTransformationModel("linear");
  aqm.setTransformationModelParams(param);
  // set-up the quant_method map
  std::vector<AbsoluteQuantitationMethod> quant_methods;
  const String feature_name = "peak_apex_int";
  // component_1
  aqm.setComponentName("ser-L.ser-L_1.Light");
  aqm.setISName("ser-L.ser-L_1.Heavy");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_2
  aqm.setComponentName("amp.amp_1.Light");
  aqm.setISName("amp.amp_1.Heavy");
  aqm.setFeatureName(feature_name); // test IS outside component_group
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);
  // component_3
  aqm.setComponentName("atp.atp_1.Light");
  aqm.setISName("atp.atp_1.Heavy");
  aqm.setFeatureName(feature_name);
  aqm.setConcentrationUnits("uM");
  quant_methods.push_back(aqm);

  AbsoluteQuantitation absquant;
  // set-up the class parameters
  Param absquant_params;
  absquant_params.setValue("min_points", 4);
  absquant_params.setValue("max_bias", 30.0);
  absquant_params.setValue("min_correlation_coefficient", 0.9);
  absquant_params.setValue("max_iters", 100);
  absquant_params.setValue("outlier_detection_method", "iter_jackknife");
  absquant_params.setValue("use_chauvenet", "false");
  absquant.setParameters(absquant_params);
  absquant.setQuantMethods(quant_methods);

  // set up the standards
  std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> components_concentrations;
  components_concentrations["ser-L.ser-L_1.Light"] = make_serL_standards();
  components_concentrations["amp.amp_1.Light"] = make_amp_standards();
  components_concentrations["atp.atp_1.Light"] = make_atp_standards();

  absquant.optimizeSingleCalibrationCurve("ser-L.ser-L_1.Light", components_concentrations.at("ser-L.ser-L_1.Light"));
  absquant.optimizeSingleCalibrationCurve("amp.amp_1.Light", components_concentrations.at("amp.amp_1.Light"));
  absquant.optimizeSingleCalibrationCurve("atp.atp_1.Light", components_concentrations.at("atp.atp_1.Light"));
  std::map<String, AbsoluteQuantitationMethod> quant_methods_map = absquant.getQuantMethodsAsMap();

  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][0].actual_concentration, 0.04);
  TEST_REAL_SIMILAR(components_concentrations["ser-L.ser-L_1.Light"][8].actual_concentration, 40.0);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("slope"), 0.9011392589);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("intercept"), 1.87018507);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getCorrelationCoefficient(), 0.999320072);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getLLOQ(), 0.04);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getULOQ(), 200);
  TEST_EQUAL(quant_methods_map["ser-L.ser-L_1.Light"].getNPoints(), 11);

  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("slope"), 0.95799683);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("intercept"), -1.047543387);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getCorrelationCoefficient(), 0.99916926);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getLLOQ(), 0.02);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getULOQ(), 40.0);
  TEST_EQUAL(quant_methods_map["amp.amp_1.Light"].getNPoints(), 11);

  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][3].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][3].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("slope"), 0.623040824);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("intercept"), 0.36130172586);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getCorrelationCoefficient(), 0.998208402);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getLLOQ(), 0.02);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getULOQ(), 40.0);
  TEST_EQUAL(quant_methods_map["atp.atp_1.Light"].getNPoints(), 6);
END_SECTION

/////////////////////////////
/* Protected Members **/
/////////////////////////////

START_SECTION((std::vector<AbsoluteQuantitationStandards::featureConcentration> extractComponents_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      const std::vector<size_t>& component_concentrations_indices)))

  AbsoluteQuantitation_test absquant;
  // make the components_concentrations
  static const double arrx1[] = { 1.1, 2.0, 3.3, 3.9, 4.9, 6.2  };
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = { 0.9, 1.9, 3.0, 3.7, 5.2, 6.1  };
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );
  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x1.size(); ++i)
  {
    component.setMetaValue("native_id","component" + std::to_string(i));
    component.setMetaValue("peak_apex_int",y1[i]);
    IS_component.setMetaValue("native_id","IS" + std::to_string(i));
    IS_component.setMetaValue("peak_apex_int",1.0);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = x1[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  // make the indices to extract
  static const size_t arrx2[] = { 0, 1, 3  };
  std::vector<size_t> component_concentrations_indices(arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );

  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sub = absquant.extractComponents_(
    component_concentrations, component_concentrations_indices);

  TEST_EQUAL(component_concentrations_sub[0].feature.getMetaValue("native_id"), "component0");
  TEST_REAL_SIMILAR(component_concentrations_sub[0].actual_concentration, 1.1);

  TEST_EQUAL(component_concentrations_sub[1].feature.getMetaValue("native_id"), "component1");
  TEST_REAL_SIMILAR(component_concentrations_sub[1].actual_concentration, 2.0);

  TEST_EQUAL(component_concentrations_sub[2].feature.getMetaValue("native_id"), "component3");
  TEST_REAL_SIMILAR(component_concentrations_sub[2].actual_concentration, 3.9);

END_SECTION

START_SECTION((int jackknifeOutlierCandidate_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params)))

  AbsoluteQuantitation_test absquant;

  static const double arrx1[] = { 1.1, 2.0,3.3,3.9,4.9,6.2  };
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = { 0.9, 1.9,3.0,3.7,5.2,6.1  };
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );
  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x1.size(); ++i)
  {
    component.setMetaValue("native_id","component");
    component.setMetaValue("peak_apex_int",y1[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",1.0);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = x1[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  String feature_name = "peak_apex_int";

  // set-up the model and params
  // y = m*x + b
  // x = (y - b)/m
  Param transformation_model_params;
 //   String transformation_model = "TransformationModelLinear";
  String transformation_model = "linear";

  int c1 = absquant.jackknifeOutlierCandidate_(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);
  TEST_EQUAL(c1,4);

  static const double arrx2[] = { 1,2,3,4,5,6  };
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = { 1,2,3,4,5,6};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );
  component_concentrations.clear();
  for (size_t i = 0; i < x2.size(); ++i)
  {
    component.setMetaValue("native_id","component");
    component.setMetaValue("peak_apex_int",y2[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",1.0);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = x2[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  int c2 = absquant.jackknifeOutlierCandidate_(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);
  TEST_EQUAL(c2,0);

END_SECTION

START_SECTION((int residualOutlierCandidate_(
  const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params)))

  AbsoluteQuantitation_test absquant;

  static const double arrx1[] = { 1.1, 2.0,3.3,3.9,4.9,6.2  };
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = { 0.9, 1.9,3.0,3.7,5.2,6.1  };
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );
  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x1.size(); ++i)
  {
    component.setMetaValue("native_id","component");
    component.setMetaValue("peak_apex_int",y1[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",1.0);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = x1[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  String feature_name = "peak_apex_int";

  // set-up the model and params
  // y = m*x + b
  // x = (y - b)/m
  Param transformation_model_params;
  // String transformation_model = "TransformationModelLinear";
  String transformation_model = "linear";

  int c1 = absquant.residualOutlierCandidate_(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);
  TEST_EQUAL(c1,4);

  static const double arrx2[] = { 1,2,3,4,5,6  };
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = { 1,2,3,4,5,6};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );
  component_concentrations.clear();
  for (size_t i = 0; i < x2.size(); ++i)
  {
    component.setMetaValue("native_id","component");
    component.setMetaValue("peak_apex_int",y2[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",1.0);
    component_concentration.feature = component;
    component_concentration.IS_feature = IS_component;
    component_concentration.actual_concentration = x2[i];
    component_concentration.IS_actual_concentration = 1.0;
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration);
  }

  int c2 = absquant.residualOutlierCandidate_(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params);
  TEST_EQUAL(c2,0);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
