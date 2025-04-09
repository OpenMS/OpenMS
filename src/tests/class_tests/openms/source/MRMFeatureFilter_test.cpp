// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureFilter* ptr = nullptr;
MRMFeatureFilter* nullPointer = nullptr;

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

//START_SECTION(template <typename T> bool checkRange(T const& value, T const& value_l, T const& value_u))
//{
//  MRMFeatureFilter mrmff;
//  // tests
//  TEST_EQUAL(mrmff.checkRange(2.0, 1.0, 2.0), true);
//  TEST_EQUAL(mrmff.checkRange(0.0, 1.0, 2.0), false);
//  TEST_EQUAL(mrmff.checkRange(3.0, 1.0, 2.0), false);
//  TEST_EQUAL(mrmff.checkRange(2, 1, 2), true);
//  TEST_EQUAL(mrmff.checkRange(0, 1, 2), false);
//  TEST_EQUAL(mrmff.checkRange(3, 1, 2), false);
//}
//END_SECTION
//
//START_SECTION(template <typename T> void updateRange(T const& value, T & value_l, T & value_u))
//{
//  MRMFeatureFilter mrmff;
//  double value_l = 0;
//  double value_u = 12;
//  // tests
//  mrmff.updateRange(6.0, value_l, value_u);
//  TEST_EQUAL(value_l, 0);
//  TEST_EQUAL(value_u, 12);
//  mrmff.updateRange(13.0, value_l, value_u);
//  TEST_EQUAL(value_l, 0);
//  TEST_EQUAL(value_u, 13);
//  mrmff.updateRange(-1.0, value_l, value_u);
//  TEST_EQUAL(value_l, -1);
//  TEST_EQUAL(value_u, 13);
//}
//END_SECTION
//
//START_SECTION(template <typename T> void setRange(T const& value, T & value_l, T & value_u))
//{
//  MRMFeatureFilter mrmff;
//  double value_l = 2;
//  double value_u = 12;
//  // tests
//  mrmff.setRange(6.0, value_l, value_u);
//  TEST_EQUAL(value_l, 0);
//  TEST_EQUAL(value_u, 6);
//  value_l = 2;
//  value_u = 12;
//  mrmff.setRange(13.0, value_l, value_u);
//  TEST_EQUAL(value_l, 2);
//  TEST_EQUAL(value_u, 13);
//  mrmff.setRange(-1.0, value_l, value_u);
//  TEST_EQUAL(value_l, -1);
//  TEST_EQUAL(value_u, 0);
//}
//END_SECTION
//
//START_SECTION(template <typename T> void initRange(T const& value, T & value_l, T & value_u))
//{
//  MRMFeatureFilter mrmff;
//  double value_l = 2;
//  double value_u = 12;
//  // tests
//  mrmff.initRange(6.0, value_l, value_u);
//  TEST_EQUAL(value_l, 6);
//  TEST_EQUAL(value_u, 6);
//}
//END_SECTION

START_SECTION(double calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name) const)
{
  MRMFeatureFilter mrmff;
  String feature_name = "peak_apex_int";
  double inf = std::numeric_limits<double>::infinity();
  // dummy features
  OpenMS::Feature component_1, component_2;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id","component1");
  component_2.setMetaValue(feature_name, 5.0);
  component_2.setMetaValue("native_id","component2");
  // tests
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_1,component_2,feature_name),1.0);
  component_2.setMetaValue(feature_name, 0.0);
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_1,component_2,feature_name),inf);
  // dummy features
  OpenMS::Feature component_3, component_4;
  component_3.setMetaValue("peak_area", 5.0);
  component_3.setMetaValue("native_id","component3");
  component_4.setMetaValue("peak_area", 5.0);
  component_4.setMetaValue("native_id","component4");
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_1,component_4,feature_name),5.0);
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_3,component_4,feature_name),0.0);
  // feature_name == "intensity"
  // feature_name == "intensity"
  Feature component_5, component_6, component_7, component_8;
  feature_name = "intensity";
  component_5.setMetaValue("native_id", "component5");
  component_6.setMetaValue("native_id", "component6");
  component_5.setIntensity(3.0);
  component_6.setIntensity(4.0);
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_5, component_6, feature_name), 0.75);
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_6, component_5, feature_name), 1.33333333333333);
  component_7.setMetaValue("native_id", "component7");
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_5, component_7, feature_name), inf);
  TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_5, component_8, feature_name), 3.0);
}
END_SECTION

START_SECTION(bool checkMetaValue(
  const Feature & component,
  const String & meta_value_key,
  const double & meta_value_l,
  const double & meta_value_u,
  bool & key_exists
) const)
{
  MRMFeatureFilter mrmff;
  bool metavalue_exists;

  //make test feature
  String feature_name = "peak_apex_int";
  OpenMS::Feature component_1;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id","component1");

  // test parameters
  double meta_value_l(4.0), meta_value_u(6.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists), true); // pass case
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, 6.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u,metavalue_exists), true); // edge pass case
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, 3.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists), false); // fail case
  TEST_EQUAL(metavalue_exists, true);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, "peak_area", meta_value_l, meta_value_u, metavalue_exists), true); // not found case
  TEST_EQUAL(metavalue_exists, false);
}
END_SECTION

START_SECTION(void updateMetaValue(
  const Feature & component,
  const String & meta_value_key,
  double & meta_value_l,
  double & meta_value_u,
  bool & key_exists
) const)
{
  MRMFeatureFilter mrmff;
  bool metavalue_exists;

  //make test feature
  String feature_name = "peak_apex_int";
  OpenMS::Feature component_1;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id", "component1");

  // test parameters
  double meta_value_l(4.0), meta_value_u(6.0);
  mrmff.updateMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 4); // no change case
  TEST_EQUAL(meta_value_u, 6); // no change case
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, 7.0);
  mrmff.updateMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 4); // no change case
  TEST_EQUAL(meta_value_u, 7); // change case
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, 3.0);
  mrmff.updateMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 3); // change case
  TEST_EQUAL(meta_value_u, 7); // no change case
  TEST_EQUAL(metavalue_exists, true);
  mrmff.updateMetaValue(component_1, "peak_area", meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 3); // no change case
  TEST_EQUAL(meta_value_u, 7); // no change case
  TEST_EQUAL(metavalue_exists, false); // not found case
}
END_SECTION

START_SECTION(void setMetaValue(
  const Feature & component,
  const String & meta_value_key,
  double & meta_value_l,
  double & meta_value_u,
  bool & key_exists
) const)
{
  MRMFeatureFilter mrmff;
  bool metavalue_exists;

  //make test feature
  String feature_name = "peak_apex_int";
  OpenMS::Feature component_1;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id", "component1");

  // test parameters
  double meta_value_l(4.0), meta_value_u(6.0);
  mrmff.setMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 0);
  TEST_EQUAL(meta_value_u, 5);
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, 7.0);
  meta_value_l = 4.0; meta_value_u = 6.0;
  mrmff.setMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 0);
  TEST_EQUAL(meta_value_u, 7);
  TEST_EQUAL(metavalue_exists, true);
  component_1.setMetaValue(feature_name, -1.0);
  mrmff.setMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, -1);
  TEST_EQUAL(meta_value_u, 0);
  TEST_EQUAL(metavalue_exists, true);
  mrmff.setMetaValue(component_1, "peak_area", meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, -1); // no change case
  TEST_EQUAL(meta_value_u, 0); // no change case
  TEST_EQUAL(metavalue_exists, false); // not found case
}
END_SECTION

START_SECTION(void initMetaValue(
  const Feature & component,
  const String & meta_value_key,
  double & meta_value_l,
  double & meta_value_u,
  bool & key_exists
) const)
{
  MRMFeatureFilter mrmff;
  bool metavalue_exists;

  //make test feature
  String feature_name = "peak_apex_int";
  OpenMS::Feature component_1;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id", "component1");

  // test parameters
  double meta_value_l(4.0), meta_value_u(6.0);
  mrmff.initMetaValue(component_1, feature_name, meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 5);
  TEST_EQUAL(meta_value_u, 5);
  TEST_EQUAL(metavalue_exists, true);
  meta_value_l = 4.0; meta_value_u = 6.0;
  mrmff.initMetaValue(component_1, "peak_area", meta_value_l, meta_value_u, metavalue_exists);
  TEST_EQUAL(meta_value_l, 4); // no change case
  TEST_EQUAL(meta_value_u, 6); // no change case
  TEST_EQUAL(metavalue_exists, false); // not found case
}
END_SECTION

START_SECTION(countLabelsAndTransitionTypes(const Feature & component_group, const TargetedExperiment & transitions) const)
{
  MRMFeatureFilter mrmff;

  // make the feature
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setMetaValue("LabelType","Heavy");
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setMetaValue("LabelType","Light");
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setMetaValue("LabelType","Light");
  subordinates.push_back(subordinate);
  // component_1.setPeptideRef("component_group1");
  component_1.setSubordinates(subordinates);
  // // transition group 2
  // // transition 1
  // subordinate.setMetaValue("native_id","component2.1.Heavy")
  // subordinate.setMetaValue("LabelType","Heavy");
  // subordinates.push_back(subordinate);
  // // transition 2
  // subordinate.setMetaValue("native_id","component2.1.Light")
  // subordinate.setMetaValue("LabelType","Light");
  // subordinates.push_back(subordinate);

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);
  // transition group 2
  // transition 1
  transition.setNativeID("component2.1.Heavy");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component2.1.Light");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);

  std::map<String,int> test1 = mrmff.countLabelsAndTransitionTypes(component_1, transitions);
  TEST_EQUAL(test1["n_heavy"], 1);
  TEST_EQUAL(test1["n_light"], 2);
  TEST_EQUAL(test1["n_quantifying"], 2);
  TEST_EQUAL(test1["n_identifying"], 0);
  TEST_EQUAL(test1["n_detecting"], 3);
  TEST_EQUAL(test1["n_transitions"], 3);

}
END_SECTION

START_SECTION(void FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria,
  const TargetedExperiment & transitions))
{
  /** FilterFeatureMap Test 1: basic ability to flag or filter transitions or transition groups */

  MRMFeatureFilter mrmff;

  //make the FeatureMap
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",500); //should fail
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  // transition group 2
  // transition 1
  subordinate.setMetaValue("native_id","component2.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component2.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",1000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group2");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);
  // transition group 2
  // transition 1
  transition.setNativeID("component2.1.Heavy");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component2.1.Light");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);

  //make the QCs
  MRMFeatureQC qc_criteria;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  std::pair<double,double> lbub(501, 4e6);
  // transition group 1
  cgqcs.component_group_name =  "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 2;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 2;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 10.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // transition group 2
  cgqcs.component_group_name =  "component_group2";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 2; //should fail
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;  //should fail
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component2.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component2.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 10.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component2.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  // transition 2
  cqcs.component_name = "component2.1.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);

  //test flag mode
  Param params;
  params.setValue("flag_or_filter", "flag");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);

  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList()[0], "peak_apex_int");
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_message").toStringList().size(), 2);
  TEST_STRING_EQUAL(components[1].getMetaValue("QC_transition_group_message").toStringList()[0], "n_light");
  TEST_STRING_EQUAL(components[1].getMetaValue("QC_transition_group_message").toStringList()[1], "n_transitions");
  TEST_EQUAL(components[1].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[1].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 0.75);
  TEST_REAL_SIMILAR(components[1].getMetaValue("QC_transition_group_score"), 0.777777777777778);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);

  //test filter mode
  params.setValue("flag_or_filter", "filter");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);

  TEST_EQUAL(components.size(), 1);
  TEST_EQUAL(components[0].getSubordinates().size(), 2);
}
END_SECTION

START_SECTION(void FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria,
  const TargetedExperiment & transitions))
{
  /** FilterFeatureMap Test 2: tests for individual checks on each transition and transition group */
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);

  //make the QCs
  MRMFeatureQC qc_criteria;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  std::pair<double,double> lbub(500, 4e6);
  // transition group 1
  cgqcs.component_group_name =  "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 2.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  cqcs.meta_value_qc["peak_area"] = lbub;
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);

  //test all possible comparisons
  Param params;
  params.setValue("flag_or_filter", "flag");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);

  // control
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 1.0);
  components.clear();

  // RT
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(6.0); // should fail
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList()[0], "retention_time");

  components.clear();

  // Intensity
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(0.0); // should fail
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList()[0], "intensity");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 0.75);
  components.clear();

  // OverallQuality
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(0.0); //should fail
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList()[0], "overall_quality");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 0.75);
  components.clear();

  // MetaValue
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",500);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",400); //should fail
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), false);
  // TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "metaValue[peak_apex_int]");
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message").toStringList()[0], "peak_apex_int");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 0.75);
  components.clear();

  // n_heavy
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList()[0], "n_heavy");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 0.892857142857143);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 1.0);
  components.clear();

  // n_light
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList()[0], "n_light");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 0.892857142857143);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 1.0);
  components.clear();

  // n_transitions
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList()[0], "n_transitions");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 0.888888888888889);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  components.clear();

  // ion_ratio_pair
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id","component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Heavy");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType","Light");
  subordinate.setMetaValue("peak_apex_int",500);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  mrmff.FilterFeatureMap(components, qc_criteria, transitions);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getMetaValue("QC_transition_group_message").toStringList()[0], "ion_ratio_pair[component1.1.Light/component1.2.Light]");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_score"), 0.964285714285714);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_score"), 1.0);
  components.clear();
}
END_SECTION

START_SECTION(void EstimateDefaultMRMFeatureQCValues(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions))
{
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  std::vector<FeatureMap> samples;
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // sample 1
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(2.5);
  component_1.setIntensity(15000);
  component_1.setOverallQuality(300);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();
  // sample 2
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(3.0);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(200);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(3);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(400);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 2000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(3.1);
  subordinate.setIntensity(2000);
  subordinate.setOverallQuality(300);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 800);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(3);
  component_1.setIntensity(5000);
  component_1.setOverallQuality(200);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();
  // sample 3
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2);
  subordinate.setIntensity(4000);
  subordinate.setOverallQuality(150);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 6000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(300);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 100);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(2);
  component_1.setIntensity(15000);
  component_1.setOverallQuality(500);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);

  // make the expected QCs values
  std::vector<MRMFeatureQC> filter_values;
  MRMFeatureQC qc_criteria1, qc_criteria2;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;

  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 0;
  cgqcs.n_heavy_u = 0;
  cgqcs.n_light_l = 0;
  cgqcs.n_light_u = 0;
  cgqcs.n_detecting_l = 0;
  cgqcs.n_detecting_u = 0;
  cgqcs.n_quantifying_l = 0;
  cgqcs.n_quantifying_u = 0;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 0;
  cgqcs.n_transitions_l = 0;
  cgqcs.n_transitions_u = 0;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 0);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(0, 0); // should not change
  qc_criteria1.component_group_qcs.push_back(cgqcs);
  qc_criteria1.component_qcs.push_back(cqcs);
  qc_criteria2 = qc_criteria1;

  // Test calculateFilterValuesMean without initialization of values
  mrmff.EstimateDefaultMRMFeatureQCValues(samples, qc_criteria1, transitions, false);

  // transition group 1
  TEST_STRING_EQUAL(qc_criteria1.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_heavy_l, 0); // lower limits are not changed
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_heavy_u, 1);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_light_l, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_light_u, 2);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_detecting_u, 3);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_quantifying_u, 2);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_transitions_l, 0);
  TEST_EQUAL(qc_criteria1.component_group_qcs.at(0).n_transitions_u, 3);
  TEST_STRING_EQUAL(qc_criteria1.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(qc_criteria1.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(qc_criteria1.component_group_qcs.at(0).ion_ratio_l, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_group_qcs.at(0).ion_ratio_u, 60);
  TEST_STRING_EQUAL(qc_criteria1.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(qc_criteria1.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).retention_time_l, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).retention_time_u, 3);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).intensity_l, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).intensity_u, 5000);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).overall_quality_l, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).overall_quality_u, 200);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 5000);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(qc_criteria1.component_qcs.at(0).meta_value_qc["peak_area"].second, 0);

  // Test calculateFilterValuesMean
  std::vector<MRMFeatureQC> filter_values_test;
  mrmff.EstimateDefaultMRMFeatureQCValues(samples, qc_criteria2, transitions, true);

  // transition group 1
  TEST_STRING_EQUAL(qc_criteria2.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_heavy_l, 1);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_heavy_u, 1);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_light_l, 2);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_light_u, 2);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_detecting_l, 3);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_detecting_u, 3);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_quantifying_l, 2);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_quantifying_u, 2);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_transitions_l, 3);
  TEST_EQUAL(qc_criteria2.component_group_qcs.at(0).n_transitions_u, 3);
  TEST_STRING_EQUAL(qc_criteria2.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(qc_criteria2.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(qc_criteria2.component_group_qcs.at(0).ion_ratio_l, 2.5);
  TEST_REAL_SIMILAR(qc_criteria2.component_group_qcs.at(0).ion_ratio_u, 60);
  TEST_STRING_EQUAL(qc_criteria2.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(qc_criteria2.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).retention_time_l, 2);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).retention_time_u, 3);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).intensity_l, 1000);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).intensity_u, 5000);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).overall_quality_l, 100);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).overall_quality_u, 200);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 1000);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 5000);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(qc_criteria2.component_qcs.at(0).meta_value_qc["peak_area"].second, 0);
}
END_SECTION

START_SECTION(void EstimatePercRSD(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions))
{
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  std::vector<FeatureMap> samples;
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  // transition group 2
  // transition 1
  subordinate.setMetaValue("native_id", "component2.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component2.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group2");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  // simulate triplicates with idententical values
  // (sufficient to test for differences in the means/vars/rsds)
  samples.push_back(components);
  samples.push_back(components);
  samples.push_back(components);

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);
  // transition group 2
  // transition 1
  transition.setNativeID("component2.1.Heavy");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component2.1.Light");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);

  //make the QCs
  MRMFeatureQC qc_criteria;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  std::pair<double, double> lbub(500, 4e6);
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 2.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  cqcs.meta_value_qc["peak_area"] = lbub;
  qc_criteria.component_qcs.push_back(cqcs);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  cqcs.meta_value_qc["peak_area"] = lbub;
  qc_criteria.component_qcs.push_back(cqcs);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  cqcs.meta_value_qc["peak_area"] = lbub;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);

  MRMFeatureQC filter_zeros = qc_criteria;
  mrmff.EstimatePercRSD(samples, filter_zeros, transitions);
  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 0);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 0);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 0);
  // transition 2
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(1).component_name, "component1.1.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).retention_time_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).intensity_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).overall_quality_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_apex_int"].second, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_area"].second, 0);
  // transition 3
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(2).component_name, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).retention_time_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).intensity_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).overall_quality_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_apex_int"].second, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_area"].second, 0);
}
END_SECTION

START_SECTION(void EstimateBackgroundInterferences(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions))
{
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  std::vector<FeatureMap> samples;
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  // transition group 2
  // transition 1
  subordinate.setMetaValue("native_id", "component2.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component2.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group2");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  // simulate triplicates with idententical values
  // (sufficient to test for differences in the means/vars/rsds)
  samples.push_back(components);
  samples.push_back(components);
  samples.push_back(components);

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);
  // transition group 2
  // transition 1
  transition.setNativeID("component2.1.Heavy");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component2.1.Light");
  transition.setPeptideRef("component_group2");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);

  //make the QCs
  MRMFeatureQC qc_criteria;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  std::pair<double, double> lbub(500, 4e6);
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 2.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  cqcs.meta_value_qc["peak_area"] = lbub;
  qc_criteria.component_qcs.push_back(cqcs);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  qc_criteria.component_qcs.push_back(cqcs);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = lbub;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);

  MRMFeatureQC filter_zeros = qc_criteria;
  mrmff.EstimateBackgroundInterferences(samples, filter_zeros, transitions);
  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 1);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 2);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 3);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 2);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 3);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 10);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 2.5);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 5000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 100);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 5000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 500);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 4e6);
  // transition 2
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(1).component_name, "component1.1.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).retention_time_u, 2.5);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).intensity_u, 5000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).overall_quality_u, 100);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_apex_int"].second, 5000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_area"].first, 500);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(1).meta_value_qc["peak_area"].second, 4e6);
  // transition 3
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(2).component_name, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).retention_time_u, 2.5);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).intensity_u, 5000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).overall_quality_u, 100);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_apex_int"].second, 500);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_area"].first, 500);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(2).meta_value_qc["peak_area"].second, 4e6);
}
END_SECTION

START_SECTION(void accumulateFilterValues(std::vector<MRMFeatureQC>& filter_values, const std::vector<FeatureMap>& samples, const MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const)
{
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  std::vector<FeatureMap> samples;
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // sample 1
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(2.5);
  component_1.setIntensity(15000);
  component_1.setOverallQuality(300);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();
  // sample 2
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(3.0);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(200);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(3);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(400);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 2000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(3.1);
  subordinate.setIntensity(2000);
  subordinate.setOverallQuality(300);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 800);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(3);
  component_1.setIntensity(5000);
  component_1.setOverallQuality(200);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();
  // sample 3
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2);
  subordinate.setIntensity(4000);
  subordinate.setOverallQuality(150);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 6000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2);
  subordinate.setIntensity(1000);
  subordinate.setOverallQuality(300);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 100);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setRT(2);
  component_1.setIntensity(15000);
  component_1.setOverallQuality(500);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  samples.push_back(components);
  components.clear();

  // make the targeted experiment
  TargetedExperiment transitions;
  ReactionMonitoringTransition transition;
  // transition group 1
  // transition 1
  transition.setNativeID("component1.1.Heavy");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 2
  transition.setNativeID("component1.1.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(true);
  transitions.addTransition(transition);
  // transition 3
  transition.setNativeID("component1.2.Light");
  transition.setPeptideRef("component_group1");
  transition.setDetectingTransition(true);
  transition.setIdentifyingTransition(false);
  transition.setQuantifyingTransition(false);
  transitions.addTransition(transition);

  // make the expected QCs values
  std::vector<MRMFeatureQC> filter_values;
  MRMFeatureQC qc_criteria1, qc_criteria2, qc_criteria3;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;

  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 0;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 0;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 0;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 0;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 0;
  cgqcs.n_transitions_l = 0;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 10;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 2.5;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 5000;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 100;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 5000);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 4e6);
  qc_criteria1.component_group_qcs.push_back(cgqcs);
  qc_criteria1.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria1);
  // transition group 1;
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 0;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 0;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 0;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 0;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 0;
  cgqcs.n_transitions_l = 0;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 2.5;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 3;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 1000;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 200;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1000);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 4e6);
  qc_criteria2.component_group_qcs.push_back(cgqcs);
  qc_criteria2.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria2);
  // transition group 1;
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 0;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 0;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 0;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 0;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 0;
  cgqcs.n_transitions_l = 0;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 60;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 2;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 5000;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 100;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 5000);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 4e6);
  qc_criteria3.component_group_qcs.push_back(cgqcs);
  qc_criteria3.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria3);

  // Test calculateFilterValuesMean
  std::vector<MRMFeatureQC> filter_values_test;
  mrmff.accumulateFilterValues(filter_values_test, samples, qc_criteria1, transitions);
  for (size_t i = 0; i < filter_values.size(); ++i) 
  {
    // transition group 1
    TEST_STRING_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).component_group_name, filter_values.at(i).component_group_qcs.at(0).component_group_name);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_heavy_l, filter_values.at(i).component_group_qcs.at(0).n_heavy_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_heavy_u, filter_values.at(i).component_group_qcs.at(0).n_heavy_u);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_light_l, filter_values.at(i).component_group_qcs.at(0).n_light_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_light_u, filter_values.at(i).component_group_qcs.at(0).n_light_u);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_detecting_l, filter_values.at(i).component_group_qcs.at(0).n_detecting_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_detecting_u, filter_values.at(i).component_group_qcs.at(0).n_detecting_u);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_quantifying_l, filter_values.at(i).component_group_qcs.at(0).n_quantifying_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_quantifying_u, filter_values.at(i).component_group_qcs.at(0).n_quantifying_u);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_identifying_l, filter_values.at(i).component_group_qcs.at(0).n_identifying_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_identifying_u, filter_values.at(i).component_group_qcs.at(0).n_identifying_u);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_transitions_l, filter_values.at(i).component_group_qcs.at(0).n_transitions_l);
    TEST_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).n_transitions_u, filter_values.at(i).component_group_qcs.at(0).n_transitions_u);
    TEST_STRING_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).ion_ratio_pair_name_1, filter_values.at(i).component_group_qcs.at(0).ion_ratio_pair_name_1);
    TEST_STRING_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).ion_ratio_pair_name_2, filter_values.at(i).component_group_qcs.at(0).ion_ratio_pair_name_2);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_group_qcs.at(0).ion_ratio_l, filter_values.at(i).component_group_qcs.at(0).ion_ratio_l);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_group_qcs.at(0).ion_ratio_u, filter_values.at(i).component_group_qcs.at(0).ion_ratio_u);
    TEST_STRING_EQUAL(filter_values_test.at(i).component_group_qcs.at(0).ion_ratio_feature_name, filter_values.at(i).component_group_qcs.at(0).ion_ratio_feature_name);
    // transition 1
    TEST_STRING_EQUAL(filter_values_test.at(i).component_qcs.at(0).component_name, filter_values.at(i).component_qcs.at(0).component_name);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).retention_time_l, filter_values.at(i).component_qcs.at(0).retention_time_l);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).retention_time_u, filter_values.at(i).component_qcs.at(0).retention_time_u);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).intensity_l, filter_values.at(i).component_qcs.at(0).intensity_l);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).intensity_u, filter_values.at(i).component_qcs.at(0).intensity_u);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).overall_quality_l, filter_values.at(i).component_qcs.at(0).overall_quality_l);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).overall_quality_u, filter_values.at(i).component_qcs.at(0).overall_quality_u);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).meta_value_qc["peak_apex_int"].first, filter_values.at(i).component_qcs.at(0).meta_value_qc["peak_apex_int"].first);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).meta_value_qc["peak_apex_int"].second, filter_values.at(i).component_qcs.at(0).meta_value_qc["peak_apex_int"].second);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).meta_value_qc["peak_area"].first, filter_values.at(i).component_qcs.at(0).meta_value_qc["peak_area"].first);
    TEST_REAL_SIMILAR(filter_values_test.at(i).component_qcs.at(0).meta_value_qc["peak_area"].second, filter_values.at(i).component_qcs.at(0).meta_value_qc["peak_area"].second);
  }
}
END_SECTION

START_SECTION(void zeroFilterValues(MRMFeatureQC& filter_zeros, const MRMFeatureQC& filter_template) const)
{
  MRMFeatureFilter mrmff;

  //make the QCs
  MRMFeatureQC qc_criteria;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;

  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 2.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2.0;
  cqcs.retention_time_u = 3.0;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4e6;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(500, 4e6);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 4e6);
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);

  //test for all zeros
  MRMFeatureQC filter_zeros;
  mrmff.zeroFilterValues(filter_zeros, qc_criteria);
  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 0);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 0);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 0);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 0);
}
END_SECTION

START_SECTION(void calculateFilterValuesMean(MRMFeatureQC& filter_mean, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_template) const)
{ // Test Mean, Var, and PercRSD
  MRMFeatureFilter mrmff;

  //make the QCs
  std::vector<MRMFeatureQC> filter_values;
  MRMFeatureQC qc_criteria1, qc_criteria2, qc_criteria3;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;

  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 0;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 3;
  cgqcs.n_transitions_u = 3;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.5;
  cgqcs.ion_ratio_u = 2;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 2;
  cqcs.retention_time_u = 3;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 4.00E+06;
  cqcs.overall_quality_l = 100;
  cqcs.overall_quality_u = 500;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(500, 4e6);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 4e6);
  qc_criteria1.component_group_qcs.push_back(cgqcs);
  qc_criteria1.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria1);
  // transition group 1;
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 1;
  cgqcs.n_light_l = 2;
  cgqcs.n_light_u = 3;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 6;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 2;
  cgqcs.n_identifying_l = 0;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 1;
  cgqcs.n_transitions_u = 2;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.25;
  cgqcs.ion_ratio_u = 3;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 1;
  cqcs.retention_time_u = 2;
  cqcs.intensity_l = 400;
  cqcs.intensity_u = 5.00E+05;
  cqcs.overall_quality_l = 50;
  cqcs.overall_quality_u = 700;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(400, 5e5);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(400, 5e5);
  qc_criteria2.component_group_qcs.push_back(cgqcs);
  qc_criteria2.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria2);
  // transition group 1;
  cgqcs.component_group_name = "component_group1";
  cgqcs.n_heavy_l = 1;
  cgqcs.n_heavy_u = 2;
  cgqcs.n_light_l = 1;
  cgqcs.n_light_u = 2;
  cgqcs.n_detecting_l = 2;
  cgqcs.n_detecting_u = 3;
  cgqcs.n_quantifying_l = 2;
  cgqcs.n_quantifying_u = 4;
  cgqcs.n_identifying_l = 1;
  cgqcs.n_identifying_u = 3;
  cgqcs.n_transitions_l = 0;
  cgqcs.n_transitions_u = 4;
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0.4;
  cgqcs.ion_ratio_u = 2;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1;
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 1;
  cqcs.retention_time_u = 4;
  cqcs.intensity_l = 500;
  cqcs.intensity_u = 3.00E+06;
  cqcs.overall_quality_l = 10;
  cqcs.overall_quality_u = 600;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(500, 3e6);
  cqcs.meta_value_qc["peak_area"] = std::make_pair(500, 3e6);
  qc_criteria3.component_group_qcs.push_back(cgqcs);
  qc_criteria3.component_qcs.push_back(cqcs);
  filter_values.push_back(qc_criteria3);

  // Test calculateFilterValuesMean
  MRMFeatureQC filter_zeros;
  mrmff.calculateFilterValuesMean(filter_zeros, filter_values, qc_criteria1);  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 0.666666667);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 1.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 1.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 2.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 2);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 4);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 2);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 2.666666667);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 3);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 1.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 3);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 0.383333333);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 2.333333333);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 1.333333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 3);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 466.6666667);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 2500000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 53.33333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 600);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 466.6666667);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 2500000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 466.6666667);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 2500000);

  // Test calculateFilterValuesVar
  MRMFeatureQC filter_means = filter_zeros;
  mrmff.calculateFilterValuesVar(filter_zeros, filter_values, filter_means, qc_criteria1);
  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 1);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 0.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 0.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 3);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 2);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 2.333333333);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 1);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 0.015833333);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 0.333333333);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 0.333333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 1);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 3333.333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 3.25e12);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 2033.333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 10000);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 3333.333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 3.25e12);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 3333.333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 3.25e12);

  // Test calculateFilterValuesPercRSD
  MRMFeatureQC filter_vars = filter_zeros;
  mrmff.calculateFilterValuesPercRSD(filter_zeros, filter_means, filter_vars);

  // transition group 1
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).component_group_name, "component_group1");
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_heavy_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_light_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_detecting_u, 43.30127019);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_quantifying_u, 70);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_l, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_identifying_u, 0);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_l, 141);
  TEST_EQUAL(filter_zeros.component_group_qcs.at(0).n_transitions_u, 33.33333333);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_1, "component1.1.Light");
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_pair_name_2, "component1.2.Light");
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_l, 32.82536711);
  TEST_REAL_SIMILAR(filter_zeros.component_group_qcs.at(0).ion_ratio_u, 24.74358297);
  TEST_STRING_EQUAL(filter_zeros.component_group_qcs.at(0).ion_ratio_feature_name, "peak_apex_int");
  // transition 1
  TEST_STRING_EQUAL(filter_zeros.component_qcs.at(0).component_name, "component1.1.Heavy");
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_l, 43.30127019);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).retention_time_u, 33.33333333);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_l, 12.37179148);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).intensity_u, 72.11102551);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_l, 84.54843287);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).overall_quality_u, 16.66666667);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].first, 12.37179148);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_apex_int"].second, 72.11102551);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].first, 12.37179148);
  TEST_REAL_SIMILAR(filter_zeros.component_qcs.at(0).meta_value_qc["peak_area"].second, 72.11102551);
}
END_SECTION

START_SECTION(void FilterFeatureMapPercRSD(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))
{
  /** FilterFeatureMap Test 1: basic ability to flag or filter transitions or transition groups */

  MRMFeatureFilter mrmff;

  //make the FeatureMap
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500); //should fail
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  // transition group 2
  // transition 1
  subordinate.setMetaValue("native_id", "component2.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component2.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group2");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  //make the %RSD filter criteria and %RSD calculated values
  MRMFeatureQC qc_criteria, qc_rsd_values;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  // %RSD filter criteria (30% RSD for all values)
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // transition group 2
  cgqcs.component_group_name = "component_group2";
  cgqcs.ion_ratio_pair_name_1 = "component2.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component2.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component2.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 2
  cqcs.component_name = "component2.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // Calculated %RSD values
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0; 
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 100.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  // transition group 2
  cgqcs.component_group_name = "component_group2";
  cgqcs.ion_ratio_pair_name_1 = "component2.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component2.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 40.0; // should fail
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component2.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 2
  cqcs.component_name = "component2.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 10.0; 
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);

  //test flag mode
  Param params;
  params.setValue("flag_or_filter", "flag");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);

  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList()[0], "peak_apex_int");
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_%RSD_pass"), false);
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[1].getMetaValue("QC_transition_group_%RSD_message").toStringList()[0], "ion_ratio_pair[component2.1.Light/component2.2.Light]");
  TEST_EQUAL(components[1].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[1].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 0.75);
  TEST_REAL_SIMILAR(components[1].getMetaValue("QC_transition_group_%RSD_score"), 0.75);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);

  //test filter mode
  params.setValue("flag_or_filter", "filter");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);

  TEST_EQUAL(components.size(), 1);
  TEST_EQUAL(components[0].getSubordinates().size(), 2);
}
END_SECTION

START_SECTION(void FilterFeatureMapPercRSD(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))
{
  /** FilterFeatureMap Test 2: tests for individual checks on each transition and transition group */
  MRMFeatureFilter mrmff;

  //make the FeatureMap
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  //make the %RSD filter criteria and %RSD calculated values
  MRMFeatureQC qc_criteria, qc_rsd_values;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  // %RSD filter criteria (30% RSD for all values)
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // Calculated %RSD values
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);

  //test all possible comparisons
  Param params;
  params.setValue("flag_or_filter", "flag");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);

  // control
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 1.0);

  // RT
  qc_rsd_values.component_group_qcs.clear();
  qc_rsd_values.component_qcs.clear();
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 80.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList()[0], "retention_time");

  // Intensity
  qc_rsd_values.component_group_qcs.clear();
  qc_rsd_values.component_qcs.clear();
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 100.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList()[0], "intensity");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 0.75);

  // OverallQuality
  qc_rsd_values.component_group_qcs.clear();
  qc_rsd_values.component_qcs.clear();
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 300.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList()[0], "overall_quality");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 0.75);

  // MetaValue
  qc_rsd_values.component_group_qcs.clear();
  qc_rsd_values.component_qcs.clear();
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 30.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 100.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_message").toStringList()[0], "peak_apex_int");
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 0.75);

  // ion_ratio_pair
  qc_rsd_values.component_group_qcs.clear();
  qc_rsd_values.component_qcs.clear();
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.ion_ratio_pair_name_1 = "component1.1.Light";
  cgqcs.ion_ratio_pair_name_2 = "component1.2.Light";
  cgqcs.ion_ratio_l = 0;
  cgqcs.ion_ratio_u = 50.0;
  cgqcs.ion_ratio_feature_name = "peak_apex_int";
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 20.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 10.0);
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 10.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 20.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 30.0);
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.retention_time_l = 0;
  cqcs.retention_time_u = 30.0;
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10.0;
  cqcs.overall_quality_l = 0;
  cqcs.overall_quality_u = 30.0;
  cqcs.meta_value_qc["peak_apex_int"] = std::make_pair(0, 1.0);
  qc_rsd_values.component_group_qcs.push_back(cgqcs);
  qc_rsd_values.component_qcs.push_back(cqcs);
  mrmff.FilterFeatureMapPercRSD(components, qc_criteria, qc_rsd_values);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_pass"), false);
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getMetaValue("QC_transition_group_%RSD_message").toStringList()[0], "ion_ratio_pair[component1.1.Light/component1.2.Light]");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%RSD_score"), 0.75);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%RSD_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%RSD_score"), 1.0);
}
END_SECTION

START_SECTION(void FilterFeatureMapBackgroundInterference(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))
{
  /** FilterFeatureMap Test 1: basic ability to flag or filter transitions or transition groups */

  MRMFeatureFilter mrmff;

  //make the FeatureMap
  FeatureMap components;
  OpenMS::Feature component_1, subordinate;
  std::vector<OpenMS::Feature> subordinates;
  // transition group 1
  // transition 1
  subordinate.setMetaValue("native_id", "component1.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component1.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 5000);
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id", "component1.2.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 500); //should fail
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group1");
  component_1.setIntensity(5000);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();
  // transition group 2
  // transition 1
  subordinate.setMetaValue("native_id", "component2.1.Heavy");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Heavy");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id", "component2.1.Light");
  subordinate.setRT(2.5);
  subordinate.setIntensity(5000);
  subordinate.setOverallQuality(100);
  subordinate.setMetaValue("LabelType", "Light");
  subordinate.setMetaValue("peak_apex_int", 1000);
  subordinates.push_back(subordinate);
  component_1.setMetaValue("PeptideRef", "component_group2");
  component_1.setIntensity(5000);
  component_1.setSubordinates(subordinates);
  components.push_back(component_1);
  subordinates.clear();

  //make the %BackgroundInterference filter criteria and %BackgroundInterference calculated values
  MRMFeatureQC qc_criteria, qc_background_values;
  MRMFeatureQC::ComponentGroupQCs cgqcs;
  MRMFeatureQC::ComponentQCs cqcs;
  // %BackgroundInterference filter criteria (30% BackgroundInterference for all values)
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.intensity_l = 0;
  cgqcs.intensity_u = 30;
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // transition group 2
  cgqcs.component_group_name = "component_group2";
  cgqcs.intensity_l = 0;
  cgqcs.intensity_u = 30;
  // transition 1
  cqcs.component_name = "component2.1.Heavy";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  // transition 2
  cqcs.component_name = "component2.1.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 30.0;
  qc_criteria.component_group_qcs.push_back(cgqcs);
  qc_criteria.component_qcs.push_back(cqcs);
  // Calculated %BackgroundInterference values
  // transition group 1
  cgqcs.component_group_name = "component_group1";
  cgqcs.intensity_l = 0;
  cgqcs.intensity_u = 0;
  // transition 1
  cqcs.component_name = "component1.1.Heavy";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 0.0;
  // transition 2
  cqcs.component_name = "component1.1.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 0.0;
  // transition 3
  cqcs.component_name = "component1.2.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 10000.0;
  qc_background_values.component_group_qcs.push_back(cgqcs);
  qc_background_values.component_qcs.push_back(cqcs);
  // transition group 2
  cgqcs.component_group_name = "component_group2";
  cgqcs.intensity_l = 0;
  cgqcs.intensity_u = 10000;
  // transition 1
  cqcs.component_name = "component2.1.Heavy";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 0;
  // transition 2
  cqcs.component_name = "component2.1.Light";
  cqcs.intensity_l = 0;
  cqcs.intensity_u = 0;
  qc_background_values.component_group_qcs.push_back(cgqcs);
  qc_background_values.component_qcs.push_back(cqcs);

  //test flag mode
  Param params;
  params.setValue("flag_or_filter", "flag");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMapBackgroundInterference(components, qc_criteria, qc_background_values);

  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_%BackgroundInterference_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_%BackgroundInterference_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_%BackgroundInterference_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%BackgroundInterference_pass"), false);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%BackgroundInterference_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_%BackgroundInterference_message").toStringList()[0], "intensity");
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_%BackgroundInterference_pass"), false);
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_%BackgroundInterference_message").toStringList().size(), 1);
  TEST_STRING_EQUAL(components[1].getMetaValue("QC_transition_group_%BackgroundInterference_message").toStringList()[0], "intensity");
  TEST_EQUAL(components[1].getSubordinates()[0].getMetaValue("QC_transition_%BackgroundInterference_pass"), true);
  TEST_EQUAL(components[1].getSubordinates()[1].getMetaValue("QC_transition_%BackgroundInterference_pass"), true);
  TEST_REAL_SIMILAR(components[0].getMetaValue("QC_transition_group_%BackgroundInterference_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[0].getMetaValue("QC_transition_%BackgroundInterference_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[1].getMetaValue("QC_transition_%BackgroundInterference_score"), 1.0);
  TEST_REAL_SIMILAR(components[0].getSubordinates()[2].getMetaValue("QC_transition_%BackgroundInterference_score"), 0.0);
  TEST_REAL_SIMILAR(components[1].getMetaValue("QC_transition_group_%BackgroundInterference_score"), 0.0);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[0].getMetaValue("QC_transition_%BackgroundInterference_score"), 1.0);
  TEST_REAL_SIMILAR(components[1].getSubordinates()[1].getMetaValue("QC_transition_%BackgroundInterference_score"), 1.0);

  //test filter mode
  params.setValue("flag_or_filter", "filter");
  mrmff.setParameters(params);
  mrmff.FilterFeatureMapBackgroundInterference(components, qc_criteria, qc_background_values);

  TEST_EQUAL(components.size(), 1);
  TEST_EQUAL(components[0].getSubordinates().size(), 2);
}
END_SECTION

START_SECTION(double calculateRTDifference(Feature & component_1, Feature & component_2))
{
  MRMFeatureFilter mrmff;
  //TODO

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

