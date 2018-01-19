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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

using namespace OpenMS;
using namespace std;

class AbsoluteQuantitationMethodFile_facade : AbsoluteQuantitationMethodFile
{
  public:

    void parseHeader_(StringList & line, std::map<String, int> & headers,
    std::map<String, int> & params_headers)
    {
      AbsoluteQuantitationMethodFile::parseHeader_(line, headers, params_headers);
    }

    void parseLine_(StringList & line, std::map<String,int> & headers, 
    std::map<String,int> & params_headers,
    AbsoluteQuantitationMethod & aqm)
    {
      AbsoluteQuantitationMethodFile::parseLine_(line, headers, params_headers, aqm);
    }
};

///////////////////////////

START_TEST(AbsoluteQuantitationMethodFile, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethodFile* ptr = nullptr;
AbsoluteQuantitationMethodFile* nullPointer = nullptr;

START_SECTION((AbsoluteQuantitationMethodFile()))
	ptr = new AbsoluteQuantitationMethodFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~AbsoluteQuantitationMethodFile()))
	delete ptr;
END_SECTION

START_SECTION((void parseHeader_(StringList & line, std::map<String,int> & headers,
  std::map<String,int> & params_headers)))
  
  AbsoluteQuantitationMethodFile_facade aqmf;
  
  std::map<String,int> headers;
  std::map<String,int> params_headers;

  // header test 1
  StringList header1; 
  header1.push_back("IS_name");
  header1.push_back("component_name");
  header1.push_back("feature_name");
  header1.push_back("concentration_units");
  header1.push_back("llod");
  header1.push_back("ulod");
  header1.push_back("lloq");
  header1.push_back("uloq");
  header1.push_back("correlation_coefficient");
  header1.push_back("n_points");
  header1.push_back("transformation_model");
  header1.push_back("transformation_model_param_slope");
  header1.push_back("transformation_model_param_intercept");

  aqmf.parseHeader_(header1, headers, params_headers);

  TEST_EQUAL(headers["IS_name"], 0);
  TEST_EQUAL(headers["transformation_model"], 10);
  TEST_EQUAL(params_headers["slope"], 11);
  TEST_EQUAL(params_headers["intercept"], 12);

  headers.clear();
  params_headers.clear();
  
  // header test 2
  StringList header2; 
  header2.push_back("IS_name");
  header2.push_back("component_name");
  header2.push_back("feature_name");
  header2.push_back("concentration_units");
  // header2.push_back("llod"); //test missing value
  header2.push_back("ulod");
  header2.push_back("lloq");
  header2.push_back("uloq");
  header2.push_back("correlation_coefficient");
  header2.push_back("n_points");
  header2.push_back("transformation_model");
  header2.push_back("transformation_model_param_slope");
  header2.push_back("transformation_model_param_intercept");

  aqmf.parseHeader_(header2, headers, params_headers);

  TEST_EQUAL(headers["IS_name"], 0);
  TEST_EQUAL(headers["llod"], -1);
  TEST_EQUAL(headers["transformation_model"], 9);
  TEST_EQUAL(params_headers["slope"], 10);
  TEST_EQUAL(params_headers["intercept"], 11);
  
END_SECTION

START_SECTION((void parseLine_(StringList & line, std::map<String,int> & headers, 
  std::map<String,int> & params_headers, AbsoluteQuantitationMethod & aqm)))
  
  AbsoluteQuantitationMethodFile_facade aqmf;
  AbsoluteQuantitationMethod aqm;
  
  // headers
  std::map<String,int> headers;
  std::map<String,int> params_headers;  
  headers["IS_name"] = 0;
  headers["component_name"] = 1;
  headers["feature_name"] = 2;
  headers["concentration_units"] = 3;
  headers["llod"] = 4;
  headers["ulod"] = 5;
  headers["lloq"] = 6;
  headers["uloq"] = 7;
  headers["correlation_coefficient"] = 8;
  headers["n_points"] = 9;
  headers["transformation_model"] = 10;
  params_headers["slope"] = 11;
  params_headers["intercept"] = 12;

  // line test 1
  StringList line1; 
  line1.push_back("IS1");
  line1.push_back("component1");
  line1.push_back("feature1");
  line1.push_back("uM");
  line1.push_back("0.0");
  line1.push_back(""); //test for empty string
  line1.push_back(" 2.0  "); //test for leading and trailing white spaces
  line1.push_back("8.0");
  line1.push_back("0.99");
  line1.push_back("5");
  line1.push_back("TransformationModelLinear");
  line1.push_back("2.0");
  line1.push_back("1.0");

  aqmf.parseLine_(line1, headers, params_headers, aqm);

  String component_name = aqm.getComponentName();
  String IS_name = aqm.getISName();
  String feature_name = aqm.getFeatureName();
  TEST_EQUAL(component_name, "component1");
  TEST_EQUAL(IS_name, "IS1");
  TEST_EQUAL(feature_name, "feature1");
  double llod = aqm.getLLOD();
  double ulod = aqm.getULOD();
  TEST_REAL_SIMILAR(llod, 0.0);
  TEST_REAL_SIMILAR(ulod, 0.0);
  double lloq = aqm.getLLOQ();
  double uloq = aqm.getULOQ();
  TEST_REAL_SIMILAR(lloq, 2.0);
  TEST_REAL_SIMILAR(uloq, 8.0);
  String concentration_units = aqm.getConcentrationUnits();
  TEST_EQUAL(concentration_units, "uM");  
  int n_points = aqm.getNPoints();
  double correlation_coefficient = aqm.getCorrelationCoefficient();
  TEST_EQUAL(n_points, 5);
  TEST_REAL_SIMILAR(correlation_coefficient, 0.99);
  String transformation_model = aqm.getTransformationModel();
  Param transformation_model_params = aqm.getTransformationModelParams();
  TEST_EQUAL(transformation_model, "TransformationModelLinear");
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"),2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"),1.0);

  headers.clear();
  params_headers.clear();
  
END_SECTION

START_SECTION((void load(const String & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)))
  AbsoluteQuantitationMethodFile aqmf;
  std::vector<AbsoluteQuantitationMethod> aqm_list;

  aqmf.load(OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_1.csv"), aqm_list);
  String component_name, IS_name, feature_name;
  component_name = aqm_list[0].getComponentName();
  IS_name = aqm_list[0].getISName();
  feature_name = aqm_list[0].getFeatureName();
  TEST_EQUAL(component_name, "component1");
  TEST_EQUAL(IS_name, "IS1");
  TEST_EQUAL(feature_name, "feature1");
  String transformation_model;
  Param transformation_model_params;
  transformation_model = aqm_list[0].getTransformationModel();
  transformation_model_params = aqm_list[0].getTransformationModelParams();
  TEST_EQUAL(transformation_model, "TransformationModelLinear");
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"),2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"),1.0);
  transformation_model_params.clear();

  component_name = aqm_list[1].getComponentName();
  IS_name = aqm_list[1].getISName();
  feature_name = aqm_list[1].getFeatureName();
  TEST_EQUAL(component_name, "component2");
  TEST_EQUAL(IS_name, "IS2");
  TEST_EQUAL(feature_name, "feature2");
  transformation_model = aqm_list[1].getTransformationModel();
  transformation_model_params = aqm_list[1].getTransformationModelParams();
  TEST_EQUAL(transformation_model, "TransformationModelLinear");
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"),2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"),2.0);
  transformation_model_params.clear();
  
  component_name = aqm_list[2].getComponentName();
  IS_name = aqm_list[2].getISName();
  feature_name = aqm_list[2].getFeatureName();
  TEST_EQUAL(component_name, "component3");
  TEST_EQUAL(IS_name, "IS3");
  TEST_EQUAL(feature_name, "feature3");
  transformation_model = aqm_list[2].getTransformationModel();
  transformation_model_params = aqm_list[2].getTransformationModelParams();
  TEST_EQUAL(transformation_model, "TransformationModelLinear");
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"),1.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"),2.0);
  transformation_model_params.clear();

END_SECTION

START_SECTION(void store(const String & filename, const std::vector<AbsoluteQuantitationMethod> & aqm_list) const)
  AbsoluteQuantitationMethodFile aqmf;
  vector<AbsoluteQuantitationMethod> aqm_list1, aqm_list2;
  aqmf.load(OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_1.csv"), aqm_list1);
  aqmf.store(OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_2.csv"), aqm_list1);
  aqmf.load(OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_2.csv"), aqm_list2);
  TEST_EQUAL(aqm_list1.size(), aqm_list2.size())
  for (Size i = 0; i < aqm_list1.size(); ++i)
  {
    TEST_EQUAL(aqm_list1[i].getISName(), aqm_list2[i].getISName())
    TEST_EQUAL(aqm_list1[i].getComponentName(), aqm_list2[i].getComponentName())
    TEST_EQUAL(aqm_list1[i].getFeatureName(), aqm_list2[i].getFeatureName())
    TEST_EQUAL(aqm_list1[i].getConcentrationUnits(), aqm_list2[i].getConcentrationUnits())
    TEST_REAL_SIMILAR(aqm_list1[i].getLLOD(), aqm_list2[i].getLLOD())
    TEST_REAL_SIMILAR(aqm_list1[i].getULOD(), aqm_list2[i].getULOD())
    TEST_REAL_SIMILAR(aqm_list1[i].getLLOQ(), aqm_list2[i].getLLOQ())
    TEST_REAL_SIMILAR(aqm_list1[i].getULOQ(), aqm_list2[i].getULOQ())
    TEST_REAL_SIMILAR(aqm_list1[i].getCorrelationCoefficient(), aqm_list2[i].getCorrelationCoefficient())
    TEST_EQUAL(aqm_list1[i].getNPoints(), aqm_list2[i].getNPoints())
    TEST_EQUAL(aqm_list1[i].getTransformationModel(), aqm_list2[i].getTransformationModel())
    const Param& tm_params1 = aqm_list1[i].getTransformationModelParams();
    const Param& tm_params2 = aqm_list2[i].getTransformationModelParams();
    TEST_EQUAL(tm_params1.size(), tm_params2.size())
    TEST_REAL_SIMILAR(tm_params1.getValue("slope"), tm_params2.getValue("slope"))
    TEST_REAL_SIMILAR(tm_params1.getValue("intercept"), tm_params2.getValue("intercept"))
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
