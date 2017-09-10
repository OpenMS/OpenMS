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

AbsoluteQuantitation* ptr = 0;
AbsoluteQuantitation* nullPointer = 0;

START_SECTION((AbsoluteQuantitation()))
	ptr = new AbsoluteQuantitation();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~AbsoluteQuantitation()))
	delete ptr;
END_SECTION

START_SECTION((double calculateRatio(Feature & component_1, Feature & component_2, std::string feature_name)))
  AbsoluteQuantitation absquant;
  std::string feature_name = "peak_apex_int";
  double inf = 1.0/0.0;
  // dummy features
  OpenMS::Feature component_1, component_2;
  component_1.setMetaValue(feature_name, 5.0);
  component_2.setMetaValue(feature_name, 5.0);
  // tests
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_1,component_2,feature_name),1.0);
  component_2.setMetaValue(feature_name, 0.0);
  TEST_REAL_SIMILAR(absquant.calculateRatio(component_1,component_2,feature_name),inf);
END_SECTION

START_SECTION((double calculateBias(double & actual_concentration, double & calculated_concentration)))
  AbsoluteQuantitation absquant;
  double actual_concentration = 5.0;
  double calculated_concentration = 5.0;
  TEST_REAL_SIMILAR(absquant.calculateBias(actual_concentration,calculated_concentration),0.0);
  calculated_concentration = 4.0;
  TEST_REAL_SIMILAR(absquant.calculateBias(actual_concentration,calculated_concentration),20.0);
END_SECTION

START_SECTION((double applyCalibration(Feature & component,
  Feature & IS_component,
  std::string & feature_name,
  std::string & transformation_model,
  Param & transformation_model_params)))

  AbsoluteQuantitation absquant;

  // set-up the features
  Feature component, IS_component;
  component.setMetaValue("native_id","component");
  component.setMetaValue("peak_apex_int",2.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue("peak_apex_int",2.0);
  std::string feature_name = "peak_apex_int";

  // set-up the model and params
  std::string transformation_model;
  Param param;
  transformation_model = "TransformationModelLinear";  
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);

  TEST_REAL_SIMILAR(absquant.applyCalibration(component,
    IS_component,
    feature_name,
    transformation_model,
    param),1.0);
END_SECTION

START_SECTION((void quantifyComponents(std::vector<FeatureMap>& unknowns)))

  AbsoluteQuantitation absquant;

  // set-up the unknowns
  std::vector<FeatureMap> unknowns;
  // set-up the unknown FeatureMap
  FeatureMap unknown_feature_map;
  // set-up the features and sub-features
  std::vector<Feature> unknown_feature_subordinates;
  Feature unknown_feature, component, IS_component;
  std::string feature_name = "peak_apex_int";
  unknown_feature.setMetaValue("PeptideRef","component_group");
  component.setMetaValue("native_id","component");
  component.setMetaValue(feature_name,2.0);
  IS_component.setMetaValue("native_id","IS");
  IS_component.setMetaValue(feature_name,2.0);
  unknown_feature_subordinates.push_back(IS_component);
  unknown_feature_subordinates.push_back(component);
  unknown_feature.setSubordinates(unknown_feature_subordinates);
  unknown_feature_map.push_back(unknown_feature);
  unknowns.push_back(unknown_feature_map);

  // set-up the quant_method map
  std::vector<AbsoluteQuantitationMethod> quant_methods;
  AbsoluteQuantitationMethod aqm;
  aqm.setComponentISFeatureNames("component","IS",feature_name);
  aqm.setConcentrationUnits("uM");
  // set-up the model and params
  std::string transformation_model;
  Param param;
  transformation_model = "TransformationModelLinear";  
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  aqm.setTransformationModel(transformation_model, param);
  quant_methods.push_back(aqm);

  absquant.setQuantMethods(quant_methods);
  absquant.quantifyComponents(unknowns);

  TEST_EQUAL(unknowns[0][0].getSubordinates()[0].getMetaValue("calculated_concentration"),"");
  TEST_STRING_EQUAL(unknowns[0][0].getSubordinates()[0].getMetaValue("concentration_units"),"");
  TEST_REAL_SIMILAR(unknowns[0][0].getSubordinates()[1].getMetaValue("calculated_concentration"),1.0);
  TEST_STRING_EQUAL(unknowns[0][0].getSubordinates()[1].getMetaValue("concentration_units"),"uM");
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
