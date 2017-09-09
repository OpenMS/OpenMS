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

START_SECTION((fitTransformationModel(const std::string & transformation_model,
  TransformationModel::DataPoints& data,
  Param& transformation_model_params)))
  
  TransformationModel::DataPoints data;
  data.push_back(make_pair(0.0, 1.0));
  data.push_back(make_pair(1.0, 2.0));
  data.push_back(make_pair(1.0, 4.0));

  AbsoluteQuantitationMethod aqm;
  std::string transformation_model;
  Param param, test;

  transformation_model = "TransformationModelLinear";  
  TransformationModelLinear tmlinear(data, param);
  //TODO: update test
  TEST_EQUAL(aqm.fitTransformationModel(transformation_model,
    data,param), test);
  
  transformation_model = "TransformationModelBSpline";
  TransformationModelBSpline tmbspline(data, param);
  TEST_EQUAL(aqm.fitTransformationModel(transformation_model,
    data,param), test);
  
  transformation_model = "TransformationModelInterpolated";
  TransformationModelInterpolated tminterpolated(data, param);
  TEST_EQUAL(aqm.fitTransformationModel(transformation_model,
    data,param), test);
  
  transformation_model = "TransformationModelLowess";
  TransformationModelLowess tmlowess(data, param);
  TEST_EQUAL(aqm.fitTransformationModel(transformation_model,
    data,param), test);
  
  transformation_model = "";
  TransformationModel tm(data, param);
  TEST_EQUAL(aqm.fitTransformationModel(transformation_model,
    data,param), test);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST