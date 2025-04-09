// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
///////////////////////////

#define TRANSITIONTSVREADER_TESTING 1 

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureSelector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const String features_path = OPENMS_GET_TEST_DATA_PATH("MRMFeatureSelector_150601_0_BloodProject01_PLT_QC_Broth-1_1_reduced.featureXML");
const String features_path_small = OPENMS_GET_TEST_DATA_PATH("MRMFeatureSelector_100ug.featureXML");
MRMFeatureSelectorScore* ptr = nullptr;
MRMFeatureSelectorScore* null_ptr = nullptr;

START_SECTION(MRMFeatureSelectorScore())
{
  ptr = new MRMFeatureSelectorScore();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MRMFeatureSelectorScore())
{
  delete ptr;
}
END_SECTION

START_SECTION(MRMFeatureSelectorScore::selectMRMFeature())
{
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);
  TEST_EQUAL(feature_map.size(), 210);

  MRMFeatureSelector::SelectorParameters parameters;

  parameters.select_transition_group = true;
  parameters.segment_window_length = -1;
  parameters.segment_step_length = -1;
  parameters.variable_type = MRMFeatureSelector::VariableType::INTEGER;
  parameters.optimal_threshold = 0.5;
  parameters.score_weights = {
    {"sn_ratio", MRMFeatureSelector::LambdaScore::LOG},
    {"peak_apices_sum", MRMFeatureSelector::LambdaScore::LOG}
  };

  MRMFeatureSelectorScore selectorScore;
  FeatureMap output_selected;
  selectorScore.selectMRMFeature(feature_map, output_selected, parameters);

  TEST_EQUAL(output_selected.size(), 37);

  const Feature& f1 = output_selected[0].getSubordinates()[0];
  TEST_REAL_SIMILAR(f1.getMetaValue("peak_apex_int"), 286.0);
  TEST_STRING_EQUAL(f1.getMetaValue("native_id").toString(), "23dpg.23dpg_1.Heavy");
  TEST_REAL_SIMILAR(f1.getRT(), 16.7592102584839);

  const Feature& f2 = output_selected[9].getSubordinates()[0];
  TEST_REAL_SIMILAR(f2.getMetaValue("peak_apex_int"), 8671.5);
  TEST_STRING_EQUAL(f2.getMetaValue("native_id").toString(), "Pool_2pg_3pg.Pool_2pg_3pg_1.Heavy");
  TEST_REAL_SIMILAR(f2.getRT(), 16.1933587515513);
}
END_SECTION

START_SECTION(removeSpaces_())
{
  MRMFeatureSelector_test selector;
  TEST_STRING_EQUAL(selector.removeSpaces_("h e ll o"), "hello");
  TEST_STRING_EQUAL(selector.removeSpaces_("hello"), "hello");
  TEST_STRING_EQUAL(selector.removeSpaces_(""), "");
  TEST_STRING_EQUAL(selector.removeSpaces_("A    B"), "AB");
}
END_SECTION

START_SECTION(constructTargTransList_())
{
  MRMFeatureSelector_test selector;
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);

  vector<pair<double, String>> time_to_name;
  map<String, vector<Feature>> feature_name_map;

  const bool select_transition_group = true;

  selector.constructTargTransList_(feature_map, time_to_name, feature_name_map, select_transition_group);

  TEST_EQUAL(time_to_name.size(), 37)
  TEST_EQUAL(feature_name_map.size(), 37)

  sort(time_to_name.begin(), time_to_name.end());

  pair<double, String> *p = nullptr;

  p = &time_to_name.front();
  TEST_REAL_SIMILAR(p->first, 0)
  TEST_STRING_EQUAL(p->second, "arg-L")

  p = &time_to_name[1];
  TEST_REAL_SIMILAR(p->first, 1.314232178)
  TEST_STRING_EQUAL(p->second, "asn-L")

  p = &time_to_name[5];
  TEST_REAL_SIMILAR(p->first, 1.421901525)
  TEST_STRING_EQUAL(p->second, "citr-L")

  p = &time_to_name[7];
  TEST_REAL_SIMILAR(p->first, 2.667749413)
  TEST_STRING_EQUAL(p->second, "AICAr")

  p = &time_to_name.back();
  TEST_REAL_SIMILAR(p->first, 99.98770892)
  TEST_STRING_EQUAL(p->second, "accoa")
}
END_SECTION

START_SECTION(weightScore_())
{
  MRMFeatureSelector_test selector;
  double score;

  score = selector.weightScore_(3413.0, MRMFeatureSelector::LambdaScore::LINEAR);
  TEST_REAL_SIMILAR(score, 3413.0)
  score = selector.weightScore_(341.0, MRMFeatureSelector::LambdaScore::INVERSE);
  TEST_REAL_SIMILAR(score, 0.002932551)
  score = selector.weightScore_(341.0, MRMFeatureSelector::LambdaScore::LOG);
  TEST_REAL_SIMILAR(score, 5.831882477)
  score = selector.weightScore_(96640.0, MRMFeatureSelector::LambdaScore::INVERSE_LOG);
  TEST_REAL_SIMILAR(score, 0.087117)
  score = selector.weightScore_(341.0, MRMFeatureSelector::LambdaScore::INVERSE_LOG10);
  TEST_REAL_SIMILAR(score, 0.394827074)
}
END_SECTION

START_SECTION(computeScore_())
{
  MRMFeatureSelector_test selector;
  double score;
  Feature feature;
  feature.setMetaValue("sn_ratio", 6.84619503982874);
  feature.setMetaValue("peak_apices_sum", 96640.0);

  score = selector.computeScore_(feature, {{"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE_LOG}});
  TEST_REAL_SIMILAR(score, 0.5198334582314795)

  score = selector.computeScore_(feature, {{"peak_apices_sum", MRMFeatureSelector::LambdaScore::INVERSE_LOG10}});
  TEST_REAL_SIMILAR(score, 0.20059549093267626)

  score = selector.computeScore_(feature, {
    {"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE_LOG},
    {"peak_apices_sum", MRMFeatureSelector::LambdaScore::INVERSE_LOG10}
  });
  TEST_REAL_SIMILAR(score, 0.10427624775717449)

  // Checks for bad input
  feature.setMetaValue("sn_ratio", 0.0);
  feature.setMetaValue("peak_apices_sum", 0.0);
  feature.setMetaValue("var_xcorr_coelution", -1.0);

  score = selector.computeScore_(feature, { {"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE} });
  TEST_REAL_SIMILAR(score, 1.0)

  score = selector.computeScore_(feature, { {"peak_apices_sum", MRMFeatureSelector::LambdaScore::LOG} });
  TEST_REAL_SIMILAR(score, 1.0)

  score = selector.computeScore_(feature, { {"var_xcorr_coelution", MRMFeatureSelector::LambdaScore::LOG} });
  TEST_REAL_SIMILAR(score, 1.0)
}
END_SECTION

START_SECTION(batchMRMFeaturesQMIP() integer) // integer variable type
{
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path_small, feature_map);

  MRMFeatureSelector::SelectorParameters params1;
  params1.nn_threshold = 4;
  params1.locality_weight = "false";
  params1.select_transition_group = "true";
  params1.segment_window_length = 8;
  params1.segment_step_length = 4;
  params1.variable_type = MRMFeatureSelector::VariableType::INTEGER;
  params1.optimal_threshold = 0.5;
  params1.score_weights = {
    // {"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE_LOG},
    // {"peak_apices_sum", MRMFeatureSelector::LambdaScore::INVERSE_LOG10}
    {"sn_ratio", MRMFeatureSelector::LambdaScore::LINEAR},
    {"peak_apices_sum", MRMFeatureSelector::LambdaScore::LINEAR}
  };

  MRMFeatureSelector::SelectorParameters params2 = params1;
  params2.segment_window_length = -1;
  params2.segment_step_length = -1;

  FeatureMap output_selected;
  MRMBatchFeatureSelector::batchMRMFeaturesQMIP(feature_map, output_selected, {params1, params2});

  sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){
    return a.getMetaValue("PeptideRef").toString() < b.getMetaValue("PeptideRef").toString(); });

  TEST_EQUAL(output_selected.size(), 8);
  TEST_STRING_EQUAL(output_selected[0].getMetaValue("PeptideRef"), "5-HTP");
  TEST_REAL_SIMILAR(output_selected[0].getRT(), 2.03215546258545);
  TEST_STRING_EQUAL(output_selected[1].getMetaValue("PeptideRef"), "Acetylserotonin");
  TEST_REAL_SIMILAR(output_selected[1].getRT(), 5.07082551965332);
  TEST_STRING_EQUAL(output_selected[2].getMetaValue("PeptideRef"), "Acetyltryptamine");
  TEST_REAL_SIMILAR(output_selected[2].getRT(), 6.96528256036377);
  TEST_STRING_EQUAL(output_selected[3].getMetaValue("PeptideRef"), "Melatonin");
  TEST_REAL_SIMILAR(output_selected[3].getRT(), 6.96528256036377);
  TEST_STRING_EQUAL(output_selected[4].getMetaValue("PeptideRef"), "Riboflavin");
  TEST_REAL_SIMILAR(output_selected[4].getRT(), 5.07082551965332);
  TEST_STRING_EQUAL(output_selected[5].getMetaValue("PeptideRef"), "Serotonin");
  TEST_REAL_SIMILAR(output_selected[5].getRT(), 1.78708603594971);
  TEST_STRING_EQUAL(output_selected[6].getMetaValue("PeptideRef"), "Tryptamine");
  TEST_REAL_SIMILAR(output_selected[6].getRT(), 3.43251273956299);
  TEST_STRING_EQUAL(output_selected[7].getMetaValue("PeptideRef"), "Tryptophan");
  TEST_REAL_SIMILAR(output_selected[7].getRT(), 3.43251273956299);

  // DEBUG
  // sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){ return a.getRT() < b.getRT(); });
  cout << "\n\nSTART DEBUG INFO" << endl;
  for (size_t i = 0; i < output_selected.size(); ++i)
  {
    const Feature& f = output_selected[i];
    cout << "[" << i << "]\t\t" << f.getMetaValue("PeptideRef") << "\t\t" << f.getRT() << endl;
    for (size_t j = 0; j < f.getSubordinates().size(); ++j)
    {
      cout << "[" << i << "][" << j << "]\t\t" << f.getSubordinates()[j].getMetaValue("native_id") << "\t\t" << f.getSubordinates()[j].getMetaValue("peak_apex_int") << endl;
    }
  }
  cout << "END   DEBUG INFO\n" << endl;
}
END_SECTION

START_SECTION(batchMRMFeaturesQMIP() continuous) // continuous variable type
{
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);

  MRMFeatureSelector::SelectorParameters params1;
  params1.nn_threshold = 4;
  params1.locality_weight = false;
  params1.select_transition_group = true;
  params1.segment_window_length = 8;
  params1.segment_step_length = 4;
  params1.variable_type = MRMFeatureSelector::VariableType::CONTINUOUS;
  params1.optimal_threshold = 0.5;
  params1.score_weights = {
    // {"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE_LOG},
    // {"peak_apices_sum", MRMFeatureSelector::LambdaScore::INVERSE_LOG10}
    {"sn_ratio", MRMFeatureSelector::LambdaScore::LINEAR},
    {"peak_apices_sum", MRMFeatureSelector::LambdaScore::LINEAR}
  };

  MRMFeatureSelector::SelectorParameters params2 = params1;
  params2.segment_window_length = -1;
  params2.segment_step_length = -1;

  FeatureMap output_selected;
  MRMBatchFeatureSelector::batchMRMFeaturesQMIP(feature_map, output_selected, {params1, params2});

  sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){
    return a.getMetaValue("PeptideRef").toString() < b.getMetaValue("PeptideRef").toString(); });

  /// @todo WARNING: This test is flaky on clang vs gcc. Deactivated.
  //TEST_EQUAL(output_selected.size(), 11);

  const Feature* f = &output_selected[0].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 0.0);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "AICAr.AICAr_1.Heavy");
  TEST_REAL_SIMILAR(f->getRT(), 1.19311977717082);

  f = &output_selected[1].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 0.0);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Hexose_Pool_fru_glc-D.Hexose_Pool_fru_glc-D_1.Heavy");
  TEST_REAL_SIMILAR(f->getRT(), 1.52517738800049);

  f = &output_selected[2].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 318.5);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Lcystin.Lcystin_1.Heavy");
  TEST_REAL_SIMILAR(f->getRT(), 0.796409679158529);

  /// @todo WARNING: This test is flaky on clang vs gcc. Deactivated.
  //f = &output_selected[10].getSubordinates()[0];
  //TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 0.0);
  //TEST_STRING_EQUAL(f->getMetaValue("native_id"), "cytd.cytd_1.Heavy");
  //TEST_REAL_SIMILAR(f->getRT(), 1.4385963780721);

  // // DEBUG
  // // sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){ return a.getRT() < b.getRT(); });
  // for (const Feature& f : output_selected)
  // {
  //   cout << f.getMetaValue("PeptideRef") << "\t" << f << endl;
  // }
}
END_SECTION

START_SECTION(batchMRMFeaturesQMIP() continuous but smaller experiment) // continuous variable type
{
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path_small, feature_map);

  MRMFeatureSelector::SelectorParameters params1;
  params1.nn_threshold = 4;
  params1.locality_weight = false;
  params1.select_transition_group = true;
  params1.segment_window_length = 8;
  params1.segment_step_length = 4;
  params1.variable_type = MRMFeatureSelector::VariableType::CONTINUOUS;
  params1.optimal_threshold = 0.5;
  params1.score_weights = {
    // {"sn_ratio", MRMFeatureSelector::LambdaScore::INVERSE_LOG},
    // {"peak_apices_sum", MRMFeatureSelector::LambdaScore::INVERSE_LOG10}
    {"sn_ratio", MRMFeatureSelector::LambdaScore::LINEAR},
    {"peak_apices_sum", MRMFeatureSelector::LambdaScore::LINEAR}
  };

  MRMFeatureSelector::SelectorParameters params2 = params1;
  params2.segment_window_length = -1;
  params2.segment_step_length = -1;

  FeatureMap output_selected;
  MRMBatchFeatureSelector::batchMRMFeaturesQMIP(feature_map, output_selected, {params1, params2});

  sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){
    return a.getMetaValue("PeptideRef").toString() < b.getMetaValue("PeptideRef").toString(); });

  TEST_EQUAL(output_selected.size(), 8);

  const Feature* f = &output_selected[0].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 29.5228353632885);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "5-HTP");
  TEST_REAL_SIMILAR(f->getRT(), 2.03215546258545);

  f = &output_selected[1].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 30.7684884945637);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Acetylserotonin");
  TEST_REAL_SIMILAR(f->getRT(), 5.07082551965332);

  f = &output_selected[2].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 28.4325753928028);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Acetyltryptamine");
  TEST_REAL_SIMILAR(f->getRT(), 6.96528256036377);

  f = &output_selected[3].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 28.4325753928028);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Melatonin");
  TEST_REAL_SIMILAR(f->getRT(), 6.96528256036377);

  f = &output_selected[4].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 30.7684884945637);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Riboflavin");
  TEST_REAL_SIMILAR(f->getRT(), 5.07082551965332);

  f = &output_selected[5].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 22.6054459245013);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Serotonin");
  TEST_REAL_SIMILAR(f->getRT(), 1.78708603594971);

  f = &output_selected[6].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 37.9693079695627);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Tryptamine");
  TEST_REAL_SIMILAR(f->getRT(), 3.43251273956299);

  f = &output_selected[7].getSubordinates()[0];
  TEST_REAL_SIMILAR(f->getMetaValue("peak_apex_int"), 37.9693079695627);
  TEST_STRING_EQUAL(f->getMetaValue("native_id"), "Tryptophan");
  TEST_REAL_SIMILAR(f->getRT(), 3.43251273956299);

  // // DEBUG
  // // sort(output_selected.begin(), output_selected.end(), [](const Feature& a, const Feature& b){ return a.getRT() < b.getRT(); });
  cout << "\n\nSTART DEBUG INFO" << endl;
  for (size_t i = 0; i < output_selected.size(); ++i)
  {
    const Feature& f = output_selected[i];
    cout << "[" << i << "]\t\t" << f.getMetaValue("PeptideRef") << "\t\t" << f.getRT() << endl;
    for (size_t j = 0; j < f.getSubordinates().size(); ++j)
    {
      cout << "[" << i << "][" << j << "]\t\t" << f.getSubordinates()[j].getMetaValue("native_id") << "\t\t" << f.getSubordinates()[j].getMetaValue("peak_apex_int") << endl;
    }
  }
  cout << "END   DEBUG INFO\n" << endl;
}
END_SECTION

END_TEST
