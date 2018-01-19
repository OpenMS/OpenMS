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

// NOTE: private (tested when public)
// START_SECTION(template <typename T> bool checkRange(T const& value, T const& value_l, T const& value_u))
// {
//   MRMFeatureFilter mrmff;
//   // tests
//   TEST_EQUAL(mrmff.checkRange(2.0, 1.0, 2.0), true);
//   TEST_EQUAL(mrmff.checkRange(0.0, 1.0, 2.0), false);
//   TEST_EQUAL(mrmff.checkRange(3.0, 1.0, 2.0), false);
//   TEST_EQUAL(mrmff.checkRange(2, 1, 2), true);
//   TEST_EQUAL(mrmff.checkRange(0, 1, 2), false);
//   TEST_EQUAL(mrmff.checkRange(3, 1, 2), false);

// }
// END_SECTION

START_SECTION(String uniqueJoin(std::vector<String>& str_vec, const String& delim))
{
  MRMFeatureFilter mrmff;
  const String str_vec_c[] = {"hello", "hello", "bye", "bye"};
  std::vector<String> str_vec(str_vec_c, str_vec_c + sizeof(str_vec_c) / sizeof(str_vec_c[0]));

  // tests
  String delim = ";";
  TEST_EQUAL(mrmff.uniqueJoin(str_vec, delim), "bye;hello");
  delim = "||";
  TEST_EQUAL(mrmff.uniqueJoin(str_vec, delim), "bye||hello");

}
END_SECTION

START_SECTION(double calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name))
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

}
END_SECTION

START_SECTION(bool checkMetaValue(const Feature & component, const String & meta_value_key, const double & meta_value_l, const double & meta_value_u))
{
  MRMFeatureFilter mrmff;

  //make test feature
  String feature_name = "peak_apex_int";
  OpenMS::Feature component_1;
  component_1.setMetaValue(feature_name, 5.0);
  component_1.setMetaValue("native_id","component1");

  // test parameters
  double meta_value_l(4.0), meta_value_u(6.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), true); // pass case
  component_1.setMetaValue(feature_name, 6.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), true); // edge pass case
  component_1.setMetaValue(feature_name, 3.0);
  TEST_EQUAL(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), false); // fail case
  TEST_EQUAL(mrmff.checkMetaValue(component_1, "peak_area", meta_value_l, meta_value_u), true);  // not found case

}
END_SECTION

START_SECTION((std::map<String,int> countLabelsAndTransitionTypes(const Feature & component_group, const TargetedExperiment & transitions)))
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
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "metaValue[peak_apex_int]");
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_pass"), false);
  TEST_EQUAL(components[1].getMetaValue("QC_transition_group_message"), "n_light;n_transitions");
  TEST_EQUAL(components[1].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[1].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  
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
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "retention_time");

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
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "intensity");
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
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "overall_quality");
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
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_message"), "metaValue[peak_apex_int]");
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
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message"), "n_heavy");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
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
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message"), "n_light");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
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
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message"), "n_transitions");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
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
  TEST_EQUAL(components[0].getMetaValue("QC_transition_group_message"), "ion_ratio_pair[component1.1.Light/component1.2.Light]");
  TEST_EQUAL(components[0].getSubordinates()[0].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[1].getMetaValue("QC_transition_pass"), true);
  TEST_EQUAL(components[0].getSubordinates()[2].getMetaValue("QC_transition_pass"), true);
  components.clear();

}
END_SECTION

START_SECTION(void FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment))
{
  MRMFeatureFilter mrmff;
  //TODO

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

