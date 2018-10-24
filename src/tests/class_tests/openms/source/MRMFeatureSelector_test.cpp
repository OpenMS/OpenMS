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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureScheduler.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
///////////////////////////

#define TRANSITIONTSVREADER_TESTING 1 

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureSelector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const String features_path = OPENMS_GET_TEST_DATA_PATH("MRMFeatureSelector_150601_0_BloodProject01_PLT_QC_Broth-1_1.featureXML");
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

START_SECTION(getParameters().getValue("nn_threshold"))
{
  MRMFeatureSelectorScore selectorScore;
  TEST_REAL_SIMILAR(selectorScore.getParameters().getValue("nn_threshold"), 4.0)
}
END_SECTION

START_SECTION(setNNThreshold())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_REAL_SIMILAR(selectorScore.getNNThreshold(), 4.0)
  selectorScore.setNNThreshold(5.0);
  TEST_REAL_SIMILAR(selectorScore.getNNThreshold(), 5.0)
}
END_SECTION

START_SECTION(getLocalityWeight())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_EQUAL(selectorScore.getLocalityWeight(), false)
  selectorScore.setLocalityWeight(true);
  TEST_EQUAL(selectorScore.getLocalityWeight(), true)
}
END_SECTION

START_SECTION(getSelectTransitionGroup())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_EQUAL(selectorScore.getSelectTransitionGroup(), true)
  selectorScore.setSelectTransitionGroup(false);
  TEST_EQUAL(selectorScore.getSelectTransitionGroup(), false)
}
END_SECTION

START_SECTION(getSegmentWindowLength())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_REAL_SIMILAR(selectorScore.getSegmentWindowLength(), 8.0)
  selectorScore.setSegmentWindowLength(7.0);
  TEST_REAL_SIMILAR(selectorScore.getSegmentWindowLength(), 7.0)
}
END_SECTION

START_SECTION(getSegmentStepLength())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_REAL_SIMILAR(selectorScore.getSegmentStepLength(), 4.0)
  selectorScore.setSegmentStepLength(3.0);
  TEST_REAL_SIMILAR(selectorScore.getSegmentStepLength(), 3.0)
}
END_SECTION

START_SECTION(getSelectHighestCount())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_EQUAL(selectorScore.getSelectHighestCount(), false)
  selectorScore.setSelectHighestCount(true);
  TEST_EQUAL(selectorScore.getSelectHighestCount(), true)
}
END_SECTION

START_SECTION(getVariableType())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_STRING_EQUAL(selectorScore.getVariableType(), "continuous")
  selectorScore.setVariableType("integer");
  TEST_STRING_EQUAL(selectorScore.getVariableType(), "integer")
}
END_SECTION

START_SECTION(getOptimalThreshold())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_REAL_SIMILAR(selectorScore.getOptimalThreshold(), 0.5)
  selectorScore.setOptimalThreshold(0.6);
  TEST_REAL_SIMILAR(selectorScore.getOptimalThreshold(), 0.6)
}
END_SECTION

START_SECTION(select_MRMFeature())
{
  const char* s_integer = MRMFeatureSelector::s_integer;
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);
  TEST_EQUAL(feature_map.size(), 703);

  MRMFeatureSelectorScore selectoreScore;

  Param param;
  param.setValue("nn_threshold", 4);
  param.setValue("locality_weight", "true");
  param.setValue("select_transition_group", "true");
  param.setValue("segment_window_length", -1);
  param.setValue("segment_step_length", -1);
  param.setValue("select_highest_count", "false");
  param.setValue("variable_type", s_integer);
  param.setValue("optimal_threshold", 1.0);
  selectoreScore.setParameters(param);

  FeatureMap output_selected;
  selectoreScore.select_MRMFeature(feature_map, output_selected);
  TEST_EQUAL(output_selected.size(), 117);
  TEST_REAL_SIMILAR(output_selected[0].getSubordinates()[0].getMetaValue("peak_apex_int"), 286.0);                        // NOTE: same result as python, but assert failing
  TEST_STRING_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("native_id").toString(), "23dpg.23dpg_1.Heavy");
  TEST_REAL_SIMILAR(output_selected[0].getSubordinates()[0].getRT(), 16.7592102584839);                                   // NOTE: same result as python, but assert failing
  TEST_REAL_SIMILAR(output_selected[50].getSubordinates()[0].getMetaValue("peak_apex_int"), 391.5);                       // NOTE: same result as python, but assert failing
  TEST_STRING_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("native_id").toString(), "f1p.f1p_1.Heavy");
  TEST_REAL_SIMILAR(output_selected[50].getSubordinates()[0].getRT(), 8.53021852213542);                                  // NOTE: same result as python, but assert failing
}
END_SECTION


START_SECTION(remove_spaces())
{
  MRMFeatureSelectorScore selectorScore;
  TEST_STRING_EQUAL(selectorScore.remove_spaces("h e ll o"), "hello");
  TEST_STRING_EQUAL(selectorScore.remove_spaces("hello"), "hello");
  TEST_STRING_EQUAL(selectorScore.remove_spaces(""), "");
  TEST_STRING_EQUAL(selectorScore.remove_spaces("A    B"), "AB");
}
END_SECTION

START_SECTION(schedule_MRMFeaturesQMIP())
{
  const char* s_continuous = MRMFeatureSelector::s_continuous;
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);

  MRMFeatureScheduler* ptrQMIP = new MRMFeatureScheduler();

  std::vector<double> nn_thresholds {4, 4};
  std::vector<String>   locality_weights {"false", "false", "false", "true"};
  std::vector<String>   select_transition_groups {"true", "true", "true", "true"};
  std::vector<double> segment_window_lengths {8, -1};
  std::vector<double> segment_step_lengths {4, -1};
  std::vector<String>   select_highest_counts {"false", "false", "false", "false"};
  std::vector<String> variable_types {s_continuous, s_continuous, s_continuous, s_continuous};
  std::vector<double> optimal_thresholds {0.5, 0.5, 0.5, 0.5};

  ptrQMIP->setNNThresholds(nn_thresholds);
  ptrQMIP->setLocalityWeights(locality_weights);
  ptrQMIP->setSelectTransitionGroups(select_transition_groups);
  ptrQMIP->setSegmentWindowLengths(segment_window_lengths);
  ptrQMIP->setSegmentStepLengths(segment_step_lengths);
  ptrQMIP->setSelectHighestCounts(select_highest_counts);
  ptrQMIP->setVariableTypes(variable_types);
  ptrQMIP->setOptimalThresholds(optimal_thresholds);

  FeatureMap output_selected;
  ptrQMIP->schedule_MRMFeaturesQMIP(feature_map, output_selected);

  TEST_EQUAL(output_selected.size(), 49);                                                                      // TODO: fails
  TEST_REAL_SIMILAR(output_selected[0].getSubordinates()[0].getMetaValue("peak_apex_int"), 262623.5);          // TODO: fails
  TEST_STRING_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("native_id"), "23dpg.23dpg_1.Heavy"); // TODO: fails
  TEST_REAL_SIMILAR(output_selected[0].getSubordinates()[0].getRT(), 15.8944563381195);                        // TODO: fails
  // TEST_REAL_SIMILAR(output_selected[50].getSubordinates()[0].getMetaValue("peak_apex_int"), 1080.0);
  // TEST_STRING_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("native_id"), "oxa.oxa_1.Heavy");
  // TEST_REAL_SIMILAR(output_selected[50].getSubordinates()[0].getRT(), 13.4963475631714);
}
END_SECTION

END_TEST
