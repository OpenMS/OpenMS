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
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

#include <OpenMS/FORMAT/MRMFeatureQCFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureFilter* ptr = 0;
MRMFeatureFilter* nullPointer = 0;

START_SECTION(MRMFeatureFilter())
{
  ptr = new MRMFeatureFilter();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureFilter())
{
  delete ptr;
}
END_SECTION

START_SECTION((void parseHeader(StringList & line, std::map<std::string,int> & headers,
    std::map<std::string,int> & params_headers)))
    //TODO

    MRMFeatureQCFile mrmfqcfile;
    
    std::map<std::string,int> headers;
    std::map<std::string,int> params_headers;
  
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
    header1.push_back("actual_concentration");
    header1.push_back("n_points");
    header1.push_back("transformation_model");
    header1.push_back("transformation_model_param_slope");
    header1.push_back("transformation_model_param_intercept");
  
    // mrmfqcfile.parseHeader(header1, headers, params_headers);
  
    // TEST_EQUAL(headers["IS_name"], 0);
    // TEST_EQUAL(headers["transformation_model"], 11);
    // TEST_EQUAL(params_headers["slope"], 12);
    // TEST_EQUAL(params_headers["intercept"], 13);
  
    // headers.clear();
    // params_headers.clear();
    
    // // header test 2
    // StringList header2; 
    // header2.push_back("IS_name");
    // header2.push_back("component_name");
    // header2.push_back("feature_name");
    // header2.push_back("concentration_units");
    // // header2.push_back("llod"); //test missing value
    // header2.push_back("ulod");
    // header2.push_back("lloq");
    // header2.push_back("uloq");
    // header2.push_back("correlation_coefficient");
    // header2.push_back("actual_concentration");
    // header2.push_back("n_points");
    // header2.push_back("transformation_model");
    // header2.push_back("transformation_model_param_slope");
    // header2.push_back("transformation_model_param_intercept");
  
    // mrmfqcfile.parseHeader(header2, headers, params_headers);
  
    // TEST_EQUAL(headers["IS_name"], 0);
    // TEST_EQUAL(headers["llod"], -1);
    // TEST_EQUAL(headers["transformation_model"], 10);
    // TEST_EQUAL(params_headers["slope"], 11);
    // TEST_EQUAL(params_headers["intercept"], 12);
    
  END_SECTION
  
  START_SECTION((void parseLine(StringList & line, std::map<std::string,int> & headers, 
    std::map<std::string,int> & params_headers, MRMFeatureQC & mrmfqc)))
    
    MRMFeatureQCFile mrmfqcfile;
    MRMFeatureQC mrmfqc;
    
    // headers
    std::map<std::string,int> headers;
    std::map<std::string,int> params_headers;  
    headers["IS_name"] = 0;
    headers["component_name"] = 1;
    headers["feature_name"] = 2;
    headers["concentration_units"] = 3;
    headers["llod"] = 4;
    headers["ulod"] = 5;
    headers["lloq"] = 6;
    headers["uloq"] = 7;
    headers["correlation_coefficient"] = 8;
    headers["actual_concentration"] = 9;
    headers["n_points"] = 10;
    headers["transformation_model"] = 11;
    params_headers["slope"] = 12;
    params_headers["intercept"] = 13;
  
    // line test 1
    StringList line1; 
    line1.push_back("IS1");
    line1.push_back("component1");
    line1.push_back("feature1");
    line1.push_back("uM");
    line1.push_back(0.0);
    line1.push_back(10.0);
    line1.push_back(2.0);
    line1.push_back(8.0);
    line1.push_back(0.99);
    line1.push_back(1.0);
    line1.push_back(5);
    line1.push_back("TransformationModelLinear");
    line1.push_back(2.0);
    line1.push_back(1.0);
  
    mrmfqcfile.parseLine(line1, headers, params_headers, mrmfqc);
  
    // std::string component_name, IS_name, feature_name;
    // mrmfqc.getComponentISFeatureNames(component_name, IS_name, feature_name);
    // TEST_EQUAL(component_name, "component1");
    // TEST_EQUAL(IS_name, "IS1");
    // TEST_EQUAL(feature_name, "feature1");
  
    headers.clear();
    
  END_SECTION
  
  START_SECTION((void load(const String & filename, MRMFeatureQC & mrmfqc)))
    MRMFeatureQCFile mrmfqcfile;
    MRMFeatureQC mrmfqc;
  
    // mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_1.csv"), mrmfqc);
    // first line
    // TEST_EQUAL(component_name, "component1");
    // TEST_EQUAL(IS_name, "IS1");
    // TEST_EQUAL(feature_name, "feature1");
    
    // second line
    // TEST_EQUAL(component_name, "component2");
    // TEST_EQUAL(IS_name, "IS2");
    // TEST_EQUAL(feature_name, "feature2");

    // third line
    // TEST_EQUAL(component_name, "component3");
    // TEST_EQUAL(IS_name, "IS3");
    // TEST_EQUAL(feature_name, "feature3");
  
  END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

