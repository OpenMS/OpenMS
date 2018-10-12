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

START_SECTION(MRMFeatureSelectorScore())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  MRMFeatureSelectorScore* null_ptr = 0;
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MRMFeatureSelectorScore())
{
  MRMFeatureSelectorScore* ptr = 0;
  delete ptr;
}
END_SECTION


START_SECTION(getParameters().getValue("nn_threshold"))
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getParameters().getValue("nn_threshold"), 4.0)
}
END_SECTION

START_SECTION(setNNThreshold())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getNNThreshold(), 4.0)
  ptr->setNNThreshold(5.0);
  TEST_EQUAL(ptr->getNNThreshold(), 5.0)
}
END_SECTION

START_SECTION(getLocalityWeight())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getLocalityWeight(), false)
  ptr->setLocalityWeight(true);
  TEST_EQUAL(ptr->getLocalityWeight(), true)
}
END_SECTION

START_SECTION(getSelectTransitionGroup())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getSelectTransitionGroup(), true)
  ptr->setSelectTransitionGroup(false);
  TEST_EQUAL(ptr->getSelectTransitionGroup(), false)
}
END_SECTION

START_SECTION(getSegmentWindowLength())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getSegmentWindowLength(), 8.0)
  ptr->setSegmentWindowLength(7.0);
  TEST_EQUAL(ptr->getSegmentWindowLength(), 7.0)
}
END_SECTION

START_SECTION(getSegmentStepLength())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getSegmentStepLength(), 4.0)
  ptr->setSegmentStepLength(3.0);
  TEST_EQUAL(ptr->getSegmentStepLength(), 3.0)
}
END_SECTION

START_SECTION(getSelectHighestCount())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getSelectHighestCount(), false)
  ptr->setSelectHighestCount(true);
  TEST_EQUAL(ptr->getSelectHighestCount(), true)
}
END_SECTION

START_SECTION(getVariableType())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getVariableType(), "continuous")
  ptr->setVariableType("integer");
  TEST_EQUAL(ptr->getVariableType(), "integer")
}
END_SECTION

START_SECTION(getOptimalThreshold())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->getOptimalThreshold(), 0.5)
  ptr->setOptimalThreshold(0.6);
  TEST_EQUAL(ptr->getOptimalThreshold(), 0.6)
}
END_SECTION

START_SECTION(select_MRMFeature())
{
  FeatureMap feature_map;
  FeatureXMLFile feature_file;
  feature_file.load(features_path, feature_map);
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();

  Param param;
  param.setValue("nn_threshold", 4);
  param.setValue("locality_weight", "true");
  param.setValue("select_transition_group", "true");
  param.setValue("segment_window_length", -1);
  param.setValue("segment_step_length", -1);
  param.setValue("select_highest_count", "false");
  param.setValue("variable_type", s_integer);
  param.setValue("optimal_threshold", 1.0);
  ptr->setParameters(param);

  FeatureMap output_selected;
  ptr->select_MRMFeature(feature_map, output_selected);
  std::cout << output_selected.size() << std::endl;
  TEST_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("peak_apex_int"), 0.0);
  TEST_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("native_id").toString(), "23dpg.23dpg_1.Heavy");
  TEST_EQUAL(output_selected[0].getSubordinates()[0].getRT(), 17.2147079447428);
  TEST_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("peak_apex_int"), 0.0);
  TEST_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("native_id").toString(), "f1p.f1p_1.Heavy");
  TEST_EQUAL(output_selected[50].getSubordinates()[0].getRT(), 13.4859151489258);
}
END_SECTION


START_SECTION(remove_spaces())
{
  MRMFeatureSelectorScore* ptr = new MRMFeatureSelectorScore();
  TEST_EQUAL(ptr->remove_spaces("h e ll o"), "hello");
  TEST_EQUAL(ptr->remove_spaces("hello"), "hello");
  TEST_EQUAL(ptr->remove_spaces(""), "");
  TEST_EQUAL(ptr->remove_spaces("A    B"), "AB");
}
END_SECTION

START_SECTION(schedule_MRMFeaturesQMIP())
{
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

  TEST_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("peak_apex_int"), 262623.5);
  TEST_EQUAL(output_selected[0].getSubordinates()[0].getMetaValue("native_id"), "23dpg.23dpg_1.Heavy");
  TEST_EQUAL(output_selected[0].getSubordinates()[0].getRT(), 15.8944563381195);
  // TEST_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("peak_apex_int"), 1080.0);
  // TEST_EQUAL(output_selected[50].getSubordinates()[0].getMetaValue("native_id"), "oxa.oxa_1.Heavy");
  // TEST_EQUAL(output_selected[50].getSubordinates()[0].getRT(), 13.4963475631714);
}
END_SECTION

END_TEST
