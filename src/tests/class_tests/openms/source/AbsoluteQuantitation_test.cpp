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

#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_serL_standards()
{
  // TEST 1: ser-L
  static const double arrx1[] = {2.32E+04,2.45E+04,1.78E+04,2.11E+04,1.91E+04,
    2.06E+04,1.85E+04,1.53E+04,1.40E+04,1.03E+04,1.07E+04,6.68E+03,5.27E+03,2.83E+03};
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = {4.94E+03,6.55E+03,7.37E+03,1.54E+04,2.87E+04,
    5.41E+04,1.16E+05,1.85E+05,3.41E+05,7.54E+05,9.76E+05,1.42E+06,1.93E+06,2.23E+06};
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) ); 
  static const double arrz1[] = {1.00E-02,2.00E-02,4.00E-02,1.00E-01,2.00E-01,
    4.00E-01,1.00E+00,2.00E+00,4.00E+00,1.00E+01,2.00E+01,4.00E+01,1.00E+02,2.00E+02};
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
    component_concentration.dilution_factor = 1.0;
    component_concentrations.push_back(component_concentration); 
  } 
  return component_concentrations;
}

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_amp_standards()
{
  // TEST 2: amp
  static const double arrx2[] = {2.15E+05,2.32E+05,2.69E+05,2.53E+05,2.50E+05,
  2.75E+05,2.67E+05,3.31E+05,3.15E+05,3.04E+05,3.45E+05,3.91E+05,4.62E+05,3.18E+05};
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = {4.40E+02,1.15E+03,1.53E+03,2.01E+03,4.47E+03,
  7.36E+03,2.18E+04,4.46E+04,8.50E+04,2.33E+05,5.04E+05,1.09E+06,2.54E+06,3.64E+06};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) ); 
  static const double arrz2[] = {2.00E-03,4.00E-03,8.00E-03,2.00E-02,4.00E-02,
  8.00E-02,2.00E-01,4.00E-01,8.00E-01,2.00E+00,4.00E+00,8.00E+00,2.00E+01,4.00E+01};
  std::vector<double> z2 (arrz2, arrz2 + sizeof(arrz2) / sizeof(arrz2[0]) );

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x2.size(); ++i)
  {
    component.setMetaValue("native_id","amp.amp_1.Light");
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
  return component_concentrations;
}

std::vector<AbsoluteQuantitationStandards::featureConcentration> make_atp_standards()
{
  // TEST 3: atp
  static const double arrx3[] = {8.28E+02,1.32E+03,1.57E+03,1.63E+03,1.48E+03,
  2.43E+03,4.44E+03,1.03E+04,1.75E+04,6.92E+04,1.97E+05,2.69E+05,3.20E+05,3.22E+05};
  std::vector<double> x3 (arrx3, arrx3 + sizeof(arrx3) / sizeof(arrx3[0]) );
  static const double arry3[] = {2.21E+02,4.41E+02,3.31E+02,2.21E+02,3.09E+02,
  5.96E+02,1.26E+03,2.49E+03,1.12E+04,8.79E+04,4.68E+05,1.38E+06,3.46E+06,4.19E+06};
  std::vector<double> y3 (arry3, arry3 + sizeof(arry3) / sizeof(arry3[0]) ); 
  static const double arrz3[] = {2.00E-03,4.00E-03,8.00E-03,2.00E-02,4.00E-02,
  8.00E-02,2.00E-01,4.00E-01,8.00E-01,2.00E+00,4.00E+00,8.00E+00,2.00E+01,4.00E+01};
  std::vector<double> z3 (arrz3, arrz3 + sizeof(arrz3) / sizeof(arrz3[0]) ); 

  // set-up the features
  std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations;
  AbsoluteQuantitationStandards::featureConcentration component_concentration;
  Feature component, IS_component;
  for (size_t i = 0; i < x3.size(); ++i)
  {
    component.setMetaValue("native_id","atp.atp_1.Light");
    component.setMetaValue("peak_apex_int",x3[i]);
    IS_component.setMetaValue("native_id","IS");
    IS_component.setMetaValue("peak_apex_int",y3[i]);
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

AbsoluteQuantitation* ptr = 0;
AbsoluteQuantitation* nullPointer = 0;

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
  component_concentration.actual_concentration = 1.0;
  component_concentration.IS_actual_concentration = 1.0;
  component_concentration.dilution_factor = 1.0;
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
  component_concentration.dilution_factor = 1.0;
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
  component_concentration.dilution_factor = 1.0;
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
  static const double arrz1[] = {-2, -4, -6, 2, 4, 6};
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
    component_concentration.dilution_factor = 1.0;
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

START_SECTION((void optimizeCalibrationCurveIterative(
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
  transformation_model_params.setValue("x_datum_min", -1e12);
  transformation_model_params.setValue("x_datum_max", 1e12);
  transformation_model_params.setValue("y_datum_min", -1e12);
  transformation_model_params.setValue("y_datum_max", 1e12);
  Param optimized_params;

  absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.04);
  TEST_REAL_SIMILAR(component_concentrations[8].actual_concentration, 40.0);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), -0.9011392589);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), -1.870185076);

  // TEST 2: amp
  component_concentrations = make_amp_standards();

  absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(component_concentrations[8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), -0.95799683);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), 1.047543387);

  // TEST 3: atp
  component_concentrations = make_atp_standards();

  absquant.optimizeCalibrationCurveIterative(
    component_concentrations,
    feature_name,
    transformation_model,
    transformation_model_params,
    optimized_params);

  TEST_REAL_SIMILAR(component_concentrations[0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(component_concentrations[8].actual_concentration, 0.8);
  TEST_REAL_SIMILAR(optimized_params.getValue("slope"), -0.623040824);
  TEST_REAL_SIMILAR(optimized_params.getValue("intercept"), -0.36130172586);  

END_SECTION

START_SECTION((void optimizeCalibrationCurves(AbsoluteQuantitationStandards::components_to_concentrations & components_concentrations)))
  
  AbsoluteQuantitation absquant;

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
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("slope"), -0.9011392589);
  TEST_REAL_SIMILAR(quant_methods_map["ser-L.ser-L_1.Light"].getTransformationModelParams().getValue("intercept"), -1.87018507); 

  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["amp.amp_1.Light"][8].actual_concentration, 8.0);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("slope"), -0.95799683);
  TEST_REAL_SIMILAR(quant_methods_map["amp.amp_1.Light"].getTransformationModelParams().getValue("intercept"), 1.047543387); 

  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][0].actual_concentration, 0.02);
  TEST_REAL_SIMILAR(components_concentrations["atp.atp_1.Light"][8].actual_concentration, 0.8);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("slope"), -0.623040824);
  TEST_REAL_SIMILAR(quant_methods_map["atp.atp_1.Light"].getTransformationModelParams().getValue("intercept"), -0.36130172586); 

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
