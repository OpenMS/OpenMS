// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/MRMFeatureQCFile.h>

using namespace OpenMS;
using namespace std;

class MRMFeatureQCFile_facade : MRMFeatureQCFile
{
  public:
    void pushValuesFromLine_(
      const StringList& line,
      const std::map<String, Size>& headers,
      std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
    ) const
    {
      MRMFeatureQCFile::pushValuesFromLine_(line, headers, c_qcs);
    }

    void pushValuesFromLine_(
      const StringList& line,
      const std::map<String, Size>& headers,
      std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
    ) const
    {
      MRMFeatureQCFile::pushValuesFromLine_(line, headers, cg_qcs);
    }

    void setPairValue_(
      const String& key,
      const String& value,
      const String& boundary,
      std::map<String, std::pair<double,double>>& meta_values_qc
    ) const
    {
      MRMFeatureQCFile::setPairValue_(key, value, boundary, meta_values_qc);
    }

  Int getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const Int default_value
  ) const
  {
    return MRMFeatureQCFile::getCastValue_(headers, line, header, default_value);
  }

  double getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const double default_value
  ) const
  {
    return MRMFeatureQCFile::getCastValue_(headers, line, header, default_value);
  }

  String getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const String& default_value
  ) const
  {
    return MRMFeatureQCFile::getCastValue_(headers, line, header, default_value);
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

START_SECTION(void load(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) const)
{
  MRMFeatureQCFile mrmfqcfile;
  MRMFeatureQC mrmfqc;
  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_1.csv"), mrmfqc, false); // components file
  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_2.csv"), mrmfqc, true); // component groups file
  const std::vector<MRMFeatureQC::ComponentQCs>& c_qcs = mrmfqc.component_qcs;
  const std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs = mrmfqc.component_group_qcs;
  TEST_EQUAL(c_qcs[0].component_name, "component1");
  TEST_EQUAL(c_qcs[1].component_name, "component2");
  TEST_EQUAL(c_qcs[2].component_name, "component3");
  TEST_EQUAL(c_qcs[3].component_name, "component4"); // note that the previous line within the file is skipped because component_name is empty
  TEST_REAL_SIMILAR(c_qcs[0].retention_time_l, 1.0);
  TEST_REAL_SIMILAR(c_qcs[1].retention_time_l, 3.0);
  TEST_REAL_SIMILAR(c_qcs[2].retention_time_l, 5.0);
  TEST_REAL_SIMILAR(c_qcs[3].retention_time_l, 0.0);   // default value
  TEST_REAL_SIMILAR(c_qcs[0].retention_time_u, 2.0);
  TEST_REAL_SIMILAR(c_qcs[1].retention_time_u, 4.0);
  TEST_REAL_SIMILAR(c_qcs[2].retention_time_u, 6.0);
  TEST_REAL_SIMILAR(c_qcs[3].retention_time_u, 1e12); // default value
  TEST_REAL_SIMILAR(c_qcs[0].intensity_l, 1000.0);
  TEST_REAL_SIMILAR(c_qcs[1].intensity_l, 2000.0);
  TEST_REAL_SIMILAR(c_qcs[2].intensity_l, 3000.0);
  TEST_REAL_SIMILAR(c_qcs[3].intensity_l, 0.0);        // default value
  TEST_REAL_SIMILAR(c_qcs[0].intensity_u, 1000000.0);
  TEST_REAL_SIMILAR(c_qcs[1].intensity_u, 2000000.0);
  TEST_REAL_SIMILAR(c_qcs[2].intensity_u, 3000000.0);
  TEST_REAL_SIMILAR(c_qcs[3].intensity_u, 1e12);       // default value
  TEST_REAL_SIMILAR(c_qcs[0].overall_quality_l, 2.0);
  TEST_REAL_SIMILAR(c_qcs[1].overall_quality_l, 3.0);
  TEST_REAL_SIMILAR(c_qcs[2].overall_quality_l, 4.0);
  TEST_REAL_SIMILAR(c_qcs[3].overall_quality_l, 0.0);  // default value
  TEST_REAL_SIMILAR(c_qcs[0].overall_quality_u, 5.0);
  TEST_REAL_SIMILAR(c_qcs[1].overall_quality_u, 6.0);
  TEST_REAL_SIMILAR(c_qcs[2].overall_quality_u, 7.0);
  TEST_REAL_SIMILAR(c_qcs[3].overall_quality_u, 1e12); // default value
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc.at("peak_apex_int").first, 1000.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc.at("peak_apex_int").first, 2000.0);
  TEST_REAL_SIMILAR(c_qcs[2].meta_value_qc.at("peak_apex_int").first, 3000.0);
  TEST_REAL_SIMILAR(c_qcs[3].meta_value_qc.at("peak_apex_int").first, 0.0);        // default value
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc.at("peak_apex_int").second, 1000000.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc.at("peak_apex_int").second, 2000000.0);
  TEST_REAL_SIMILAR(c_qcs[2].meta_value_qc.at("peak_apex_int").second, 3000000.0);
  TEST_REAL_SIMILAR(c_qcs[3].meta_value_qc.at("peak_apex_int").second, 1e12);      // default value
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc.at("sn_score").first, 2.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc.at("sn_score").first, 5.0);
  TEST_REAL_SIMILAR(c_qcs[2].meta_value_qc.at("sn_score").first, 10.0);
  TEST_REAL_SIMILAR(c_qcs[3].meta_value_qc.at("sn_score").first, 0.0);             // default value
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc.at("sn_score").second, 10.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc.at("sn_score").second, 20.0);
  TEST_REAL_SIMILAR(c_qcs[2].meta_value_qc.at("sn_score").second, 50.0);
  TEST_REAL_SIMILAR(c_qcs[3].meta_value_qc.at("sn_score").second, 1e12);           // default value
  TEST_EQUAL(cg_qcs[0].component_group_name, "componentGroup1");
  TEST_EQUAL(cg_qcs[1].component_group_name, "componentGroup2");
  TEST_EQUAL(cg_qcs[2].component_group_name, "componentGroup3");
  TEST_EQUAL(cg_qcs[3].component_group_name, "componentGroup5");
  TEST_EQUAL(cg_qcs[0].n_heavy_l, 1);
  TEST_EQUAL(cg_qcs[2].n_heavy_l, 3);
  TEST_EQUAL(cg_qcs[0].n_heavy_u, 2);
  TEST_EQUAL(cg_qcs[2].n_heavy_u, 4);
  TEST_EQUAL(cg_qcs[0].n_light_l, 3);
  TEST_EQUAL(cg_qcs[2].n_light_l, 5);
  TEST_EQUAL(cg_qcs[0].n_light_u, 4);
  TEST_EQUAL(cg_qcs[2].n_light_u, 6);
  TEST_EQUAL(cg_qcs[0].n_detecting_l, 5);
  TEST_EQUAL(cg_qcs[2].n_detecting_l, 7);
  TEST_EQUAL(cg_qcs[0].n_detecting_u, 6);
  TEST_EQUAL(cg_qcs[2].n_detecting_u, 8);
  TEST_EQUAL(cg_qcs[0].n_quantifying_l, 7);
  TEST_EQUAL(cg_qcs[2].n_quantifying_l, 9);
  TEST_EQUAL(cg_qcs[0].n_quantifying_u, 8);
  TEST_EQUAL(cg_qcs[2].n_quantifying_u, 10);
  TEST_EQUAL(cg_qcs[0].n_identifying_l, 9);
  TEST_EQUAL(cg_qcs[2].n_identifying_l, 11);
  TEST_EQUAL(cg_qcs[0].n_identifying_u, 10);
  TEST_EQUAL(cg_qcs[2].n_identifying_u, 12);
  TEST_EQUAL(cg_qcs[0].n_transitions_l, 11);
  TEST_EQUAL(cg_qcs[2].n_transitions_l, 13);
  TEST_EQUAL(cg_qcs[0].n_transitions_u, 12);
  TEST_EQUAL(cg_qcs[2].n_transitions_u, 14);
  TEST_EQUAL(cg_qcs[0].ion_ratio_pair_name_1, "component1");
  TEST_EQUAL(cg_qcs[2].ion_ratio_pair_name_1, "component5");
  TEST_EQUAL(cg_qcs[0].ion_ratio_pair_name_2, "component2");
  TEST_EQUAL(cg_qcs[2].ion_ratio_pair_name_2, "component6");
  TEST_REAL_SIMILAR(cg_qcs[0].ion_ratio_l, 0.5);
  TEST_REAL_SIMILAR(cg_qcs[2].ion_ratio_l, 2.5);
  TEST_REAL_SIMILAR(cg_qcs[0].ion_ratio_u, 0.6);
  TEST_REAL_SIMILAR(cg_qcs[2].ion_ratio_u, 2.6);
  TEST_EQUAL(cg_qcs[0].ion_ratio_feature_name, "feature1");
  TEST_EQUAL(cg_qcs[2].ion_ratio_feature_name, "feature3");
  TEST_REAL_SIMILAR(cg_qcs[0].retention_time_l, 1.0);
  TEST_REAL_SIMILAR(cg_qcs[1].retention_time_l, 2.0);
  TEST_REAL_SIMILAR(cg_qcs[2].retention_time_l, 3.0);
  TEST_REAL_SIMILAR(cg_qcs[0].retention_time_u, 2.0);
  TEST_REAL_SIMILAR(cg_qcs[1].retention_time_u, 3.0);
  TEST_REAL_SIMILAR(cg_qcs[2].retention_time_u, 4.0);
  TEST_REAL_SIMILAR(cg_qcs[0].intensity_l, 1000.0);
  TEST_REAL_SIMILAR(cg_qcs[1].intensity_l, 1001.0);
  TEST_REAL_SIMILAR(cg_qcs[2].intensity_l, 1002.0);
  TEST_REAL_SIMILAR(cg_qcs[0].intensity_u, 1000000.0);
  TEST_REAL_SIMILAR(cg_qcs[1].intensity_u, 1000001.0);
  TEST_REAL_SIMILAR(cg_qcs[2].intensity_u, 1000002.0);
  TEST_REAL_SIMILAR(cg_qcs[0].overall_quality_l, 2.0);
  TEST_REAL_SIMILAR(cg_qcs[1].overall_quality_l, 3.0);
  TEST_REAL_SIMILAR(cg_qcs[2].overall_quality_l, 4.0);
  TEST_REAL_SIMILAR(cg_qcs[0].overall_quality_u, 5.0);
  TEST_REAL_SIMILAR(cg_qcs[1].overall_quality_u, 6.0);
  TEST_REAL_SIMILAR(cg_qcs[2].overall_quality_u, 7.0);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc.at("peak_apex_int").first, 1000.0);
  TEST_REAL_SIMILAR(cg_qcs[2].meta_value_qc.at("peak_apex_int").first, 1002.0);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc.at("peak_apex_int").second, 1000000.0);
  TEST_REAL_SIMILAR(cg_qcs[2].meta_value_qc.at("peak_apex_int").second, 1000002.0);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc.at("sn_score").first, 2.0);
  TEST_REAL_SIMILAR(cg_qcs[2].meta_value_qc.at("sn_score").first, 10.0);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc.at("sn_score").second, 10.0);
  TEST_REAL_SIMILAR(cg_qcs[2].meta_value_qc.at("sn_score").second, 50.0);
  TEST_EQUAL(cg_qcs[3].component_group_name, "componentGroup5");
  TEST_EQUAL(cg_qcs[3].n_heavy_l, 0);
  TEST_EQUAL(cg_qcs[3].n_heavy_u, 100);
  TEST_EQUAL(cg_qcs[3].n_light_l, 0);
  TEST_EQUAL(cg_qcs[3].n_light_u, 100);
  TEST_EQUAL(cg_qcs[3].n_detecting_l, 0);
  TEST_EQUAL(cg_qcs[3].n_detecting_u, 100);
  TEST_EQUAL(cg_qcs[3].n_quantifying_l, 0);
  TEST_EQUAL(cg_qcs[3].n_quantifying_u, 100);
  TEST_EQUAL(cg_qcs[3].n_identifying_l, 0);
  TEST_EQUAL(cg_qcs[3].n_identifying_u, 100);
  TEST_EQUAL(cg_qcs[3].n_transitions_l, 0);
  TEST_EQUAL(cg_qcs[3].n_transitions_u, 100);
  TEST_EQUAL(cg_qcs[3].ion_ratio_pair_name_1, "");
  TEST_EQUAL(cg_qcs[3].ion_ratio_pair_name_2, "");
  TEST_REAL_SIMILAR(cg_qcs[3].ion_ratio_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].ion_ratio_u, 1e12);
  TEST_EQUAL(cg_qcs[3].ion_ratio_feature_name, "");
  TEST_REAL_SIMILAR(cg_qcs[3].retention_time_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].retention_time_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[3].intensity_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].intensity_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[3].overall_quality_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].overall_quality_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[3].meta_value_qc.at("peak_apex_int").first, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].meta_value_qc.at("peak_apex_int").second, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[3].meta_value_qc.at("sn_score").first, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[3].meta_value_qc.at("sn_score").second, 1e12);
}
END_SECTION

START_SECTION(void store(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group))
{
  MRMFeatureQCFile mrmfqcfile;
  MRMFeatureQC mrmfqc, mrmfqc_test;
  String file_comp = File::getTemporaryFile();
  String file_comp_group = File::getTemporaryFile();

  mrmfqcfile.store(file_comp, mrmfqc, false); // empty components file
  mrmfqcfile.store(file_comp_group, mrmfqc, true); // empty component groups file

  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_1.csv"), mrmfqc, false); // components file
  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_2.csv"), mrmfqc, true); // component groups file
  mrmfqcfile.store(file_comp, mrmfqc, false); // components file
  mrmfqcfile.store(file_comp_group, mrmfqc, true); // component groups file
  mrmfqcfile.load(file_comp, mrmfqc_test, false); // components file
  mrmfqcfile.load(file_comp_group, mrmfqc_test, true); // component groups file
  TEST_EQUAL(mrmfqc.component_qcs.size(), mrmfqc_test.component_qcs.size());
  for (size_t i = 0; i < mrmfqc.component_qcs.size(); ++i) {
    TEST_EQUAL(mrmfqc.component_qcs.at(i) == mrmfqc_test.component_qcs.at(i), true);
  }
  TEST_EQUAL(mrmfqc.component_group_qcs.size(), mrmfqc_test.component_group_qcs.size());
  for (size_t i = 0; i < mrmfqc.component_group_qcs.size(); ++i) {
    TEST_EQUAL(mrmfqc.component_group_qcs.at(i) == mrmfqc_test.component_group_qcs.at(i), true);
  }
}
END_SECTION

START_SECTION(void pushValuesFromLine_(
  const StringList& line,
  const std::map<String, Size>& headers,
  std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
) const)
{
  const std::map<String, Size> headers {
    {"component_name", 0},
    {"retention_time_l", 1},
    {"retention_time_u", 2},
    {"intensity_l", 3},
    {"intensity_u", 4},
    {"overall_quality_l", 5},
    {"overall_quality_u", 6},
    {"metaValue_peak_apex_int_l", 7},
    {"metaValue_peak_apex_int_u", 8},
    {"metaValue_sn_score_l", 9},
    {"metaValue_sn_score_u", 10}
  };

  const std::vector<String> sl1 { // all info are present
    "component1",
    "0.1", "0.2",
    "0.3", "0.4",
    "0.5", "0.6",
    "0.7", "0.8",
    "0.9", "1.0"
  };

  const std::vector<String> sl2 { // component_name is empty
    "",
    "0.1", "0.2",
    "0.3", "0.4",
    "0.5", "0.6",
    "0.7", "0.8",
    "0.9", "1.0"
  };

  const std::vector<String> sl3 { // testing defaults
    "component3",
    "", "",
    "", "",
    "", "",
    "", "",
    "", ""
  };

  MRMFeatureQCFile_facade mrmfqcfile_f;
  std::vector<MRMFeatureQC::ComponentQCs> c_qcs;

  mrmfqcfile_f.pushValuesFromLine_(sl1, headers, c_qcs);
  TEST_EQUAL(c_qcs.size(), 1);
  TEST_EQUAL(c_qcs[0].component_name, "component1");
  TEST_REAL_SIMILAR(c_qcs[0].retention_time_l, 0.1);
  TEST_REAL_SIMILAR(c_qcs[0].retention_time_u, 0.2);
  TEST_REAL_SIMILAR(c_qcs[0].intensity_l, 0.3);
  TEST_REAL_SIMILAR(c_qcs[0].intensity_u, 0.4);
  TEST_REAL_SIMILAR(c_qcs[0].overall_quality_l, 0.5);
  TEST_REAL_SIMILAR(c_qcs[0].overall_quality_u, 0.6);
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc["peak_apex_int"].first, 0.7);
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc["peak_apex_int"].second, 0.8);
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc["sn_score"].first, 0.9);
  TEST_REAL_SIMILAR(c_qcs[0].meta_value_qc["sn_score"].second, 1.0);
  mrmfqcfile_f.pushValuesFromLine_(sl2, headers, c_qcs);
  TEST_EQUAL(c_qcs.size(), 1);
  mrmfqcfile_f.pushValuesFromLine_(sl3, headers, c_qcs);
  TEST_EQUAL(c_qcs.size(), 2);
  TEST_EQUAL(c_qcs[1].component_name, "component3");
  TEST_REAL_SIMILAR(c_qcs[1].retention_time_l, 0.0);
  TEST_REAL_SIMILAR(c_qcs[1].retention_time_u, 1e12);
  TEST_REAL_SIMILAR(c_qcs[1].intensity_l, 0.0);
  TEST_REAL_SIMILAR(c_qcs[1].intensity_u, 1e12);
  TEST_REAL_SIMILAR(c_qcs[1].overall_quality_l, 0.0);
  TEST_REAL_SIMILAR(c_qcs[1].overall_quality_u, 1e12);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc["peak_apex_int"].first, 0.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc["peak_apex_int"].second, 1e12);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc["sn_score"].first, 0.0);
  TEST_REAL_SIMILAR(c_qcs[1].meta_value_qc["sn_score"].second, 1e12);
}
END_SECTION

START_SECTION(void pushValuesFromLine_(
  const StringList& line,
  const std::map<String, Size>& headers,
  std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
) const)
{
  const std::map<String, Size> headers {
    {"component_group_name", 0},
    {"n_heavy_l", 1},
    {"n_heavy_u", 2},
    {"n_light_l", 3},
    {"n_light_u", 4},
    {"n_detecting_l", 5},
    {"n_detecting_u", 6},
    {"n_quantifying_l", 7},
    {"n_quantifying_u", 8},
    {"n_identifying_l", 9},
    {"n_identifying_u", 10},
    {"n_transitions_l", 11},
    {"n_transitions_u", 12},
    {"ion_ratio_pair_name_1", 13},
    {"ion_ratio_pair_name_2", 14},
    {"ion_ratio_l", 15},
    {"ion_ratio_u", 16},
    {"ion_ratio_feature_name", 17},
    {"retention_time_l", 18},
    {"retention_time_u", 19},
    {"intensity_l", 20},
    {"intensity_u", 21},
    {"overall_quality_l", 22},
    {"overall_quality_u", 23},
    {"metaValue_peak_apex_int_l", 24},
    {"metaValue_peak_apex_int_u", 25},
    {"metaValue_sn_score_l", 26},
    {"metaValue_sn_score_u", 27}
  };

  const std::vector<String> sl1 { // all info are present
    "component_group_1",
    "1", "2",
    "3", "4",
    "5", "6",
    "7", "8",
    "9", "10",
    "11", "12",
    "ionRatioPairName1", "ionRatioPairName2",
    "1.1", "1.2",
    "ionRatioFeatureName",
    "0.1", "0.2",
    "0.3", "0.4",
    "0.5", "0.6",
    "0.7", "0.8",
    "0.9", "1.0"
  };

  const std::vector<String> sl2 { // component_name is empty
    "",
    "1", "2",
    "3", "4",
    "5", "6",
    "7", "8",
    "9", "10",
    "11", "12",
    "ionRatioPairName1", "ionRatioPairName2",
    "1.1", "1.2",
    "ionRatioFeatureName",
    "0.1", "0.2",
    "0.3", "0.4",
    "0.5", "0.6",
    "0.7", "0.8",
    "0.9", "1.0"
  };

  const std::vector<String> sl3 { // testing defaults
    "component_group_3",
    "", "",
    "", "",
    "", "",
    "", "",
    "", "",
    "", "",
    "", "",
    "", "",
    "",
    "", "",
    "", "",
    "", "",
    "", "",
    "", ""
  };

  MRMFeatureQCFile_facade mrmfqcfile_f;
  std::vector<MRMFeatureQC::ComponentGroupQCs> cg_qcs;

  mrmfqcfile_f.pushValuesFromLine_(sl1, headers, cg_qcs);
  TEST_EQUAL(cg_qcs.size(), 1);
  TEST_EQUAL(cg_qcs[0].component_group_name, "component_group_1");
  TEST_EQUAL(cg_qcs[0].n_heavy_l, 1);
  TEST_EQUAL(cg_qcs[0].n_heavy_u, 2);
  TEST_EQUAL(cg_qcs[0].n_light_l, 3);
  TEST_EQUAL(cg_qcs[0].n_light_u, 4);
  TEST_EQUAL(cg_qcs[0].n_detecting_l, 5);
  TEST_EQUAL(cg_qcs[0].n_detecting_u, 6);
  TEST_EQUAL(cg_qcs[0].n_quantifying_l, 7);
  TEST_EQUAL(cg_qcs[0].n_quantifying_u, 8);
  TEST_EQUAL(cg_qcs[0].n_identifying_l, 9);
  TEST_EQUAL(cg_qcs[0].n_identifying_u, 10);
  TEST_EQUAL(cg_qcs[0].n_transitions_l, 11);
  TEST_EQUAL(cg_qcs[0].n_transitions_u, 12);
  TEST_EQUAL(cg_qcs[0].ion_ratio_pair_name_1, "ionRatioPairName1");
  TEST_EQUAL(cg_qcs[0].ion_ratio_pair_name_2, "ionRatioPairName2");
  TEST_REAL_SIMILAR(cg_qcs[0].ion_ratio_l, 1.1);
  TEST_REAL_SIMILAR(cg_qcs[0].ion_ratio_u, 1.2);
  TEST_EQUAL(cg_qcs[0].ion_ratio_feature_name, "ionRatioFeatureName");
  TEST_REAL_SIMILAR(cg_qcs[0].retention_time_l, 0.1);
  TEST_REAL_SIMILAR(cg_qcs[0].retention_time_u, 0.2);
  TEST_REAL_SIMILAR(cg_qcs[0].intensity_l, 0.3);
  TEST_REAL_SIMILAR(cg_qcs[0].intensity_u, 0.4);
  TEST_REAL_SIMILAR(cg_qcs[0].overall_quality_l, 0.5);
  TEST_REAL_SIMILAR(cg_qcs[0].overall_quality_u, 0.6);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc["peak_apex_int"].first, 0.7);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc["peak_apex_int"].second, 0.8);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc["sn_score"].first, 0.9);
  TEST_REAL_SIMILAR(cg_qcs[0].meta_value_qc["sn_score"].second, 1.0);
  mrmfqcfile_f.pushValuesFromLine_(sl2, headers, cg_qcs);
  TEST_EQUAL(cg_qcs.size(), 1);
  mrmfqcfile_f.pushValuesFromLine_(sl3, headers, cg_qcs);
  TEST_EQUAL(cg_qcs.size(), 2);
  TEST_EQUAL(cg_qcs[1].component_group_name, "component_group_3");
  TEST_EQUAL(cg_qcs[1].n_heavy_l, 0);
  TEST_EQUAL(cg_qcs[1].n_heavy_u, 100);
  TEST_EQUAL(cg_qcs[1].n_light_l, 0);
  TEST_EQUAL(cg_qcs[1].n_light_u, 100);
  TEST_EQUAL(cg_qcs[1].n_detecting_l, 0);
  TEST_EQUAL(cg_qcs[1].n_detecting_u, 100);
  TEST_EQUAL(cg_qcs[1].n_quantifying_l, 0);
  TEST_EQUAL(cg_qcs[1].n_quantifying_u, 100);
  TEST_EQUAL(cg_qcs[1].n_identifying_l, 0);
  TEST_EQUAL(cg_qcs[1].n_identifying_u, 100);
  TEST_EQUAL(cg_qcs[1].n_transitions_l, 0);
  TEST_EQUAL(cg_qcs[1].n_transitions_u, 100);
  TEST_EQUAL(cg_qcs[1].ion_ratio_pair_name_1, "");
  TEST_EQUAL(cg_qcs[1].ion_ratio_pair_name_2, "");
  TEST_REAL_SIMILAR(cg_qcs[1].ion_ratio_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].ion_ratio_u, 1e12);
  TEST_EQUAL(cg_qcs[1].ion_ratio_feature_name, "");
  TEST_REAL_SIMILAR(cg_qcs[1].retention_time_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].retention_time_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[1].intensity_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].intensity_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[1].overall_quality_l, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].overall_quality_u, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[1].meta_value_qc["peak_apex_int"].first, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].meta_value_qc["peak_apex_int"].second, 1e12);
  TEST_REAL_SIMILAR(cg_qcs[1].meta_value_qc["sn_score"].first, 0.0);
  TEST_REAL_SIMILAR(cg_qcs[1].meta_value_qc["sn_score"].second, 1e12);
}
END_SECTION

START_SECTION(void setPairValue_(
  const String& key,
  const String& value,
  const String& boundary,
  std::map<String, std::pair<double,double>>& meta_values_qc
) const)
{
  std::map<String, std::pair<double,double>> metavalues;
  MRMFeatureQCFile_facade mrmfqcfile_f;
  mrmfqcfile_f.setPairValue_("meta1", "0.123", "u", metavalues); // first pair (initializing the upper bound)
  TEST_EQUAL(metavalues.size(), 1);
  TEST_REAL_SIMILAR(metavalues["meta1"].first, 0.0);             // default lower bound value
  TEST_REAL_SIMILAR(metavalues["meta1"].second, 0.123);
  mrmfqcfile_f.setPairValue_("meta1", "0.456", "l", metavalues); // overwrite the lower bound value
  TEST_EQUAL(metavalues.size(), 1);                              // the size of the map doesn't change
  TEST_REAL_SIMILAR(metavalues["meta1"].first, 0.456);
  TEST_REAL_SIMILAR(metavalues["meta1"].second, 0.123);
  mrmfqcfile_f.setPairValue_("meta1", "0.789", "u", metavalues); // overwrite the upper bound value
  TEST_EQUAL(metavalues.size(), 1);
  TEST_REAL_SIMILAR(metavalues["meta1"].first, 0.456);
  TEST_REAL_SIMILAR(metavalues["meta1"].second, 0.789);
  mrmfqcfile_f.setPairValue_("meta2", "0.111", "l", metavalues); // create another pair (initializing the lower bound)
  TEST_EQUAL(metavalues.size(), 2);                              // the size of the map changes
  TEST_REAL_SIMILAR(metavalues["meta2"].first, 0.111);
  TEST_REAL_SIMILAR(metavalues["meta2"].second, 1e12);           // default upper bound value
  mrmfqcfile_f.setPairValue_("meta3", "0.222", "u", metavalues); // just another pair
  TEST_EQUAL(metavalues.size(), 3);
  TEST_REAL_SIMILAR(metavalues["meta3"].first, 0.0);
  TEST_REAL_SIMILAR(metavalues["meta3"].second, 0.222);
}
END_SECTION

START_SECTION(Int getCastValue_(
  const std::map<String, Size>& headers,
  const StringList& line,
  const String& header,
  const Int default_value
) const)
{
  MRMFeatureQCFile_facade mrmfqcfile_f;
  const std::map<String, Size> headers {
    {"component_group_name", 0},
    {"n_heavy_l", 1},
    {"n_heavy_u", 2}
  };
  const StringList line1 {"componentgroup1", "2", "3"}; // all info are present
  const StringList line2 {"componentgroup2", "", "3"}; // some info is missing
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "n_heavy_l", 3), 2) // info is found, converted and returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "n_light_l", 4), 4) // the requested column is not present in the headers, default value is returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line2, "n_heavy_l", 5), 5) // the requested column is present in the headers, but the value is empty. Default value is returned
}
END_SECTION

START_SECTION(double getCastValue_(
  const std::map<String, Size>& headers,
  const StringList& line,
  const String& header,
  const double default_value
) const)
{
  MRMFeatureQCFile_facade mrmfqcfile_f;
  const std::map<String, Size> headers {
    {"component_name", 0},
    {"retention_time_l", 1},
    {"retention_time_u", 2}
  };
  const StringList line1 {"component1", "1.2", "1.3"}; // all info are present
  const StringList line2 {"component2", "", "1.3"}; // some info is missing
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "retention_time_l", 3.1), 1.2) // info is found, converted and returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "intensity_l", 4.1), 4.1) // the requested column is not present in the headers, default value is returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line2, "retention_time_l", 5.1), 5.1) // the requested column is present in the headers, but the value is empty. Default value is returned
}
END_SECTION

START_SECTION(String getCastValue_(
  const std::map<String, Size>& headers,
  const StringList& line,
  const String& header,
  const String& default_value
) const)
{
  MRMFeatureQCFile_facade mrmfqcfile_f;
  const std::map<String, Size> headers {
    {"component_name", 0},
    {"ion_ratio_feature_name", 1}
  };
  const StringList line1 {"component1", "name1"}; // all info are present
  const StringList line2 {"component2", ""}; // some info is missing
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "ion_ratio_feature_name", "name30"), "name1") // info is found, converted and returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line1, "intensity_l", "name30"), "name30") // the requested column is not present in the headers, default value is returned
  TEST_EQUAL(mrmfqcfile_f.getCastValue_(headers, line2, "ion_ratio_feature_name", "name30"), "name30") // the requested column is present in the headers, but the value is empty. Default value is returned
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
