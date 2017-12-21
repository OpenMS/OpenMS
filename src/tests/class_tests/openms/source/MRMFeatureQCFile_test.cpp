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

class MRMFeatureQCFile_facade : MRMFeatureQCFile
{
  public:

    void parseHeader_(StringList & line, std::map<String, int> & headers,
    std::map<String, int> & params_headers)
    {
      MRMFeatureQCFile::parseHeader_(line, headers, params_headers);
    }

    void parseLine_(StringList & line, std::map<String,int> & headers, 
    std::map<String,int> & params_headers,
    MRMFeatureQC & mrmfqc)
    {
      MRMFeatureQCFile::parseLine_(line, headers, params_headers, mrmfqc);
    }
};

START_TEST(MRMFeatureQCFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureQCFile* ptr = nullptr;
MRMFeatureQCFile* nullPointer = nullptr;

START_SECTION(MRMFeatureQCFile())
{
  ptr = new MRMFeatureQCFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureQCFile())
{
  delete ptr;
}
END_SECTION

START_SECTION((void parseHeader_(StringList & line, std::map<String,int> & headers,
    std::map<String,int> & params_headers)))
    //TODO

    MRMFeatureQCFile_facade mrmfqcfile;
    
    std::map<String,int> headers;
    std::map<String,int> params_headers;
  
    // header test 1
    StringList header1; 
    header1.push_back("component_name");
    header1.push_back("component_group_name");
    header1.push_back("n_heavy_l"); 
    header1.push_back("n_heavy_u");
    header1.push_back("n_light_l");
    header1.push_back("n_light_u");
    header1.push_back("n_detecting_l");
    header1.push_back("n_detecting_u");
    header1.push_back("n_quantifying_l");
    header1.push_back("n_quantifying_u");
    header1.push_back("n_identifying_l");
    header1.push_back("n_identifying_u");
    header1.push_back("n_transitions_l");
    header1.push_back("n_transitions_u");
    header1.push_back("ion_ratio_pair_name_1");
    header1.push_back("ion_ratio_pair_name_2");
    header1.push_back("ion_ratio_l");
    header1.push_back("ion_ratio_u");
    header1.push_back("ion_ratio_feature_name");
    header1.push_back("retention_time_l");
    header1.push_back("retention_time_u");
    header1.push_back("intensity_l");
    header1.push_back("intensity_u");
    header1.push_back("overall_quality_l");
    header1.push_back("overall_quality_u");
    header1.push_back("metaValue_peak_apex_int_l");
    header1.push_back("metaValue_peak_apex_int_u");
    header1.push_back("metaValue_sn_score_l");
    header1.push_back("metaValue_sn_score_u");
  
    mrmfqcfile.parseHeader_(header1, headers, params_headers);
  
    TEST_EQUAL(headers["component_name"], 0);
    TEST_EQUAL(headers["n_detecting_u"], 7);
    TEST_EQUAL(headers["overall_quality_u"], 24);
    TEST_EQUAL(params_headers["peak_apex_int_l"], 25);
    TEST_EQUAL(params_headers["sn_score_u"], 28);
  
    headers.clear();
    params_headers.clear();
    
    // header test 2
    StringList header2; 
    header2.push_back("component_name");
    header2.push_back("component_group_name");
    header2.push_back("n_heavy_l"); 
    header2.push_back("n_heavy_u");
    header2.push_back("n_light_l");
    header2.push_back("n_light_u");
    header2.push_back("n_detecting_l");
    // header2.push_back("n_detecting_u");
    header2.push_back("n_quantifying_l");
    header2.push_back("n_quantifying_u");
    header2.push_back("n_identifying_l");
    header2.push_back("n_identifying_u");
    header2.push_back("n_transitions_l");
    header2.push_back("n_transitions_u");
    header2.push_back("ion_ratio_pair_name_1");
    header2.push_back("ion_ratio_pair_name_2");
    header2.push_back("ion_ratio_l");
    header2.push_back("ion_ratio_u");
    // header2.push_back("ion_ratio_feature_name");
    header2.push_back("retention_time_l");
    header2.push_back("retention_time_u");
    header2.push_back("intensity_l");
    header2.push_back("intensity_u");
    header2.push_back("overall_quality_l");
    header2.push_back("overall_quality_u");
    header2.push_back("metaValue_peak_apex_int_l");
    header2.push_back("metaValue_peak_apex_int_u");
    header2.push_back("metaValue_sn_score_l");
    header2.push_back("metaValue_sn_score_u");
  
    mrmfqcfile.parseHeader_(header2, headers, params_headers);
    
    TEST_EQUAL(headers["component_name"], 0);
    TEST_EQUAL(headers["n_detecting_u"], -1);
    TEST_EQUAL(headers["overall_quality_u"], 22);
    TEST_EQUAL(params_headers["peak_apex_int_l"], 23);
    TEST_EQUAL(params_headers["sn_score_u"], 26);
    
  END_SECTION
  
  START_SECTION((void parseLine_(StringList & line, std::map<String,int> & headers, 
    std::map<String,int> & params_headers, MRMFeatureQC & mrmfqc)))
    
    MRMFeatureQCFile_facade mrmfqcfile;
    MRMFeatureQC mrmfqc;
    
    // headers
    std::map<String,int> headers;
    std::map<String,int> params_headers;
    headers["component_name"] = 0;
    headers["component_group_name"] = 1;
    headers["n_heavy_l"] = 2;
    headers["n_heavy_u"] = 3;
    headers["n_light_l"] = 4;
    headers["n_light_u"] = 5;
    headers["n_detecting_l"] = 6;
    headers["n_detecting_u"] = 7;
    headers["n_quantifying_l"] = 8;
    headers["n_quantifying_u"] = 9;
    headers["n_identifying_l"] = 10;
    headers["n_identifying_u"] = 11;
    headers["n_transitions_l"] = 12;
    headers["n_transitions_u"] = 13;
    headers["ion_ratio_pair_name_1"] = 14;
    headers["ion_ratio_pair_name_2"] = 15;
    headers["ion_ratio_l"] = 16;
    headers["ion_ratio_u"] = 17;
    headers["ion_ratio_feature_name"] = 18;
    headers["retention_time_l"] = 19;
    headers["retention_time_u"] = 20;
    headers["intensity_l"] = 21;
    headers["intensity_u"] = 22;
    headers["overall_quality_l"] = 23;
    headers["overall_quality_u"] = 24;
    params_headers["peak_apex_int_l"] = 25;
    params_headers["peak_apex_int_u"] = 26;
    params_headers["sn_score_l"] = 27;
    params_headers["sn_score_u"] = 28;
  
    // line test 1
    StringList line1; 
    line1.push_back("component1");
    line1.push_back("component_group1");
    line1.push_back(1);
    line1.push_back(1);
    line1.push_back(2);
    line1.push_back(2);
    line1.push_back(0);
    line1.push_back(0);
    line1.push_back(1);
    line1.push_back(1);
    line1.push_back(2);
    line1.push_back(2);
    line1.push_back(3);
    line1.push_back(3);
    line1.push_back("component1");
    line1.push_back("component2");
    line1.push_back(0.5);
    line1.push_back(0.75);
    line1.push_back("peak_apex_int");
    line1.push_back(1.0);
    line1.push_back(2.0);
    line1.push_back(1.0e3);
    line1.push_back(1.0e5);
    line1.push_back(2.0);
    line1.push_back(5.0);  
    line1.push_back(1.1e3);
    line1.push_back(1.1e5);
    line1.push_back(2.0);
    line1.push_back(10.0);

    mrmfqcfile.parseLine_(line1, headers, params_headers, mrmfqc);
  
    TEST_EQUAL(mrmfqc.component_group_qcs[0].component_group_name, "component_group1");
    TEST_EQUAL(mrmfqc.component_group_qcs[0].n_quantifying_u, 1);

    TEST_EQUAL(mrmfqc.component_qcs[0].component_name, "component1");
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].retention_time_l, 1.0);
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].overall_quality_u, 5.0);
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].meta_value_qc["peak_apex_int"].first, 1.1e3);
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].meta_value_qc["sn_score"].second, 10.0);
  
    headers.clear();
    
  END_SECTION
  
  START_SECTION((void load(const String & filename, MRMFeatureQC & mrmfqc)))
    MRMFeatureQCFile mrmfqcfile;
    MRMFeatureQC mrmfqc;
  
    mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_1.csv"), mrmfqc);
    // first line
    TEST_EQUAL(mrmfqc.component_group_qcs[0].component_group_name, "componentGroup1");
    TEST_EQUAL(mrmfqc.component_qcs[0].component_name, "component1");
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].meta_value_qc["sn_score"].second, 10.0);
    
    // second line
    TEST_EQUAL(mrmfqc.component_group_qcs[1].component_group_name, "componentGroup2");
    TEST_EQUAL(mrmfqc.component_qcs[1].component_name, "component2");
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[1].meta_value_qc["sn_score"].second, 20.0);

    // third line
    TEST_EQUAL(mrmfqc.component_group_qcs[2].component_group_name, "componentGroup3");
    TEST_EQUAL(mrmfqc.component_qcs[2].component_name, "component3");
    TEST_REAL_SIMILAR(mrmfqc.component_qcs[2].meta_value_qc["sn_score"].second, 50.0);
  
  END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

