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

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AbsoluteQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethod* ptr = 0;
AbsoluteQuantitationMethod* nullPointer = 0;
START_SECTION((AbsoluteQuantitationMethod()))
	ptr = new AbsoluteQuantitationMethod();
	TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION((~AbsoluteQuantitationMethod()))
	delete ptr;
END_SECTION

START_SECTION((bool checkLOD(const double & value)))

  AbsoluteQuantitationMethod aqm;
  double value = 2.0;

  // tests
  aqm.setLLOD(0.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value),true);
  aqm.setLLOD(0.0);
  aqm.setULOD(1.0);
  TEST_EQUAL(aqm.checkLOD(value),false);
  aqm.setLLOD(3.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value),false);
END_SECTION

START_SECTION((bool checkLOQ(const double & value)))

  AbsoluteQuantitationMethod aqm;
  double value = 2.0;

  // tests
  aqm.setLLOQ(0.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value),true);
  aqm.setLLOQ(0.0);
  aqm.setULOQ(1.0);
  TEST_EQUAL(aqm.checkLOQ(value),false);
  aqm.setLLOQ(3.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value),false);
END_SECTION

START_SECTION((Param fitTransformationModel(const String & transformation_model,
  const TransformationModel::DataPoints& data,
  const Param& transformation_model_params)))
  
  TransformationModel::DataPoints data;
  data.push_back(make_pair(0.0, 0.0));
  data.push_back(make_pair(1.0, 1.0));
  data.push_back(make_pair(2.0, 2.0));
  data.push_back(make_pair(3.0, 3.0));
  data.push_back(make_pair(4.0, 4.0));

  AbsoluteQuantitationMethod aqm;
  String transformation_model;
  Param param, test;

  transformation_model = "TransformationModelLinear";  
  TransformationModelLinear tmlinear(data, param);
  test = aqm.fitTransformationModel(transformation_model,
    data,param);
  TEST_REAL_SIMILAR(test.getValue("slope"), 1.0);
  TEST_REAL_SIMILAR(test.getValue("intercept"), 0.0);
  test.clear();
  param.clear();
  
  transformation_model = "TransformationModelBSpline";
  TransformationModelBSpline tmbspline(data, param);
  test = aqm.fitTransformationModel(transformation_model,
    data,param);
  TEST_EQUAL(test.getValue("extrapolate"), "linear");
  TEST_REAL_SIMILAR(test.getValue("wavelength"), 0.0);
  TEST_REAL_SIMILAR(test.getValue("num_nodes"), 5);
  TEST_REAL_SIMILAR(test.getValue("boundary_condition"), 2);
  test.clear();
  param.clear();
  
  transformation_model = "TransformationModelInterpolated";
  TransformationModelInterpolated tminterpolated(data, param);
  test = aqm.fitTransformationModel(transformation_model,
    data,param);
  TEST_EQUAL(test.getValue("interpolation_type"), "cspline");
  TEST_EQUAL(test.getValue("extrapolation_type"), "two-point-linear");
  test.clear();
  param.clear();
  
  transformation_model = "TransformationModelLowess";
  TransformationModelLowess tmlowess(data, param);
  test = aqm.fitTransformationModel(transformation_model,
    data,param);
  TEST_EQUAL(test.getValue("interpolation_type"), "cspline");
  TEST_REAL_SIMILAR(test.getValue("num_iterations"), 3.0);
  TEST_REAL_SIMILAR(test.getValue("span"), 2/3.0);
  test.clear();
  param.clear();
  
  transformation_model = "";
  TransformationModel tm(data, param);
  test = aqm.fitTransformationModel(transformation_model,
    data,param);
  TEST_EQUAL(test.empty(), true);
END_SECTION

START_SECTION((double evaluateTransformationModel(const String & transformation_model,
  const double& datum,
  const Param& transformation_model_params)))
  
  TransformationModel::DataPoints data;
  double datum = 2.0;
  AbsoluteQuantitationMethod aqm;
  String transformation_model;
  Param param;

  transformation_model = "TransformationModelLinear";  
  param.setValue("slope",1.0);
  param.setValue("intercept",0.0);
  TransformationModelLinear tmlinear(data, param);
  TEST_REAL_SIMILAR(aqm.evaluateTransformationModel(transformation_model,
    datum,param), 2.0);
  param.clear();
  
  // TODO:  No support yet for the following TransformationModels
  // transformation_model = "TransformationModelBSpline";
  // //TODO: update param
  // TransformationModelBSpline tmbspline(data, param);
  // //TODO: update test
  // TEST_REAL_SIMILAR(aqm.evaluateTransformationModel(transformation_model,
  //   datum,param), 2.0);
  // param.clear();
  
  // transformation_model = "TransformationModelInterpolated";
  // //TODO: update param
  // TransformationModelInterpolated tminterpolated(data, param);
  // //TODO: update test
  // TEST_REAL_SIMILAR(aqm.evaluateTransformationModel(transformation_model,
  //   datum,param), 2.0);
  // param.clear();
  
  // transformation_model = "TransformationModelLowess";
  // //TODO: update param
  // TransformationModelLowess tmlowess(data, param);
  // //TODO: update test
  // TEST_REAL_SIMILAR(aqm.evaluateTransformationModel(transformation_model,
  //   datum,param), 2.0);
  // param.clear();
  
  transformation_model = "";
  TransformationModel tm(data, param);
  TEST_REAL_SIMILAR(aqm.evaluateTransformationModel(transformation_model,
    datum,param), 2.0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST