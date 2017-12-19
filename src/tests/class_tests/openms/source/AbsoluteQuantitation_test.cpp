// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AbsoluteQuantitation, "$Id$")

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
  transformation_model = "TransformationModelLinear";  
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
  transformation_model = "TransformationModelLinear";  
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

START_SECTION((void calculateBiasAndR2(
  const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params,
  std::vector<double> & biases,
  double & r2)))
  
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
  component_concentration.actual_concentration = 1.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentrations.push_back(component_concentration);  
  // point #2
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",2.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  component_concentration.feature = component;
  component_concentration.IS_feature = IS_component;
  component_concentration.actual_concentration = 2.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentrations.push_back(component_concentration);  
  // point #3
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",3.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",1.0);
  component_concentration.feature = component;
  component_concentration.IS_feature = IS_component;
  component_concentration.actual_concentration = 3.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentrations.push_back(component_concentration);  

  String feature_name = "peak_apex_int";

  // set-up the model and params
  // y = m*x + b
  // x = (y - b)/m
  String transformation_model;
  Param param;
  transformation_model = "TransformationModelLinear"; 
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  std::vector<double> biases;
  double r2;

  absquant.calculateBiasAndR2(
    component_concentrations,
    feature_name, 
    transformation_model, 
    param,
    biases, 
    r2);
  
  TEST_REAL_SIMILAR(biases[0],0.0);  
  TEST_REAL_SIMILAR(r2,1.0);  
  
END_SECTION

START_SECTION((void optimizeCalibrationCurveBruteForce(
  const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
  const String & feature_name,
  const String & transformation_model,
  const Param & transformation_model_params)))
  
  AbsoluteQuantitation absquant;

  //TODO
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
