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

START_SECTION(void FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria,
  const TargetedExperiment & transitions))
{
  MRMFeatureFilter mrmff;
  //TODO

}
END_SECTION

START_SECTION(void FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment))
{
  MRMFeatureFilter mrmff;
  //TODO

}
END_SECTION

START_SECTION(double calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name))
{
  MRMFeatureFilter mrmff;
  String feature_name = "peak_apex_int";
  double inf = 1.0/0.0;
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

START_SECTION(double calculateRTDifference(Feature & component_1, Feature & component_2))
{
  MRMFeatureFilter mrmff;
  //TODO

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
  TEST_REAL_SIMILAR(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), true); // pass case
  component_1.setMetaValue(feature_name, 6.0);
  TEST_REAL_SIMILAR(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), true); // edge pass case
  component_1.setMetaValue(feature_name, 3.0);
  TEST_REAL_SIMILAR(mrmff.checkMetaValue(component_1, feature_name, meta_value_l, meta_value_u), false); // fail case
  TEST_REAL_SIMILAR(mrmff.checkMetaValue(component_1, "peak_area", meta_value_l, meta_value_u), true);  // not found case

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
  subordinate.setMetaValue("native_id","component1.1.Heavy")
  subordinate.setMetaValue("LabelType","Heavy");
  subordinates.push_back(subordinate);
  // transition 2
  subordinate.setMetaValue("native_id","component1.1.Light")
  subordinate.setMetaValue("LabelType","Light");
  subordinates.push_back(subordinate);
  // transition 3
  subordinate.setMetaValue("native_id","component1.2.Light")
  subordinate.setMetaValue("LabelType","Light");
  subordinates.push_back(subordinate);
  component_1.setPeptideRef("component_group1");  
  // // transition group 2
  // // transition 1
  // subordinate.setMetaValue("native_id","component2.1.Heavy")
  // subordinate.setMetaValue("LabelType","Heavy");
  // subordinates.push_back(subordinate);
  // // transition 2
  // subordinate.setMetaValue("native_id","component2.1.Light")
  // subordinate.setMetaValue("LabelType","Light");
  // subordinates.push_back(subordinate);
  // // transition 3
  // subordinate.setMetaValue("native_id","component2.2.Light")
  // subordinate.setMetaValue("LabelType","Light");
  // subordinates.push_back(subordinate);
  // component_1.setPeptideRef("component_group2");
  
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
  // // transition group 2
  // // transition 1
  // transition.setNativeID("component2.1.Heavy");
  // transition.setPeptideRef("component_group2");
  // transition.setDetectingTransition(true);
  // transition.setIdentifyingTransition(false);
  // transition.setQuantifyingTransition(true);
  // transitions.addTransition(transition);
  // // transition 2
  // transition.setNativeID("component2.1.Light");
  // transition.setPeptideRef("component_group2");
  // transition.setDetectingTransition(true);
  // transition.setIdentifyingTransition(false);
  // transition.setQuantifyingTransition(true);
  // transitions.addTransition(transition);
  // // transition 3
  // transition.setNativeID("component2.2.Light");
  // transition.setPeptideRef("component_group2");
  // transition.setDetectingTransition(true);
  // transition.setIdentifyingTransition(false);
  // transition.setQuantifyingTransition(false);
  // transitions.addTransition(transition);

  std::map<String,int> test1 = countLabelsAndTransitionTypes(component_1, transitions);
  TEST_EQUAL(test1["n_heavy"], 1);
  TEST_EQUAL(test1["n_light"], 2);
  TEST_EQUAL(test1["n_quant"], 2);
  TEST_EQUAL(test1["n_ident"], 0);
  TEST_EQUAL(test1["n_detect"], 3);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

