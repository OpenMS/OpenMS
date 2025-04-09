// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <boost/assign/std/vector.hpp>

///////////////////////////

START_TEST(OpenSwathHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

OpenSwathHelper* ptr = nullptr;
OpenSwathHelper* nullPointer = nullptr;

START_SECTION(OpenSwathHelper())
{
  ptr = new OpenSwathHelper();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathHelper())
{
  delete ptr;
}
END_SECTION

START_SECTION(static String computePrecursorId(const String& transition_group_id, int isotope))
{
  TEST_EQUAL(OpenSwathHelper::computePrecursorId("tr_gr2", 0), "tr_gr2_Precursor_i0")
  TEST_EQUAL(OpenSwathHelper::computePrecursorId("tr_gr2__test", 0), "tr_gr2__test_Precursor_i0")
}
END_SECTION

START_SECTION(static String computeTransitionGroupId(const String& precursor_id))
{
  TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId("tr_gr2_Precursor_i0"), "tr_gr2")
  TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId("tr_gr2__test_Precursor_i0"), "tr_gr2__test")
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenMS::TargetedExperiment &targeted_exp, OpenMS::TargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  TargetedExperiment exp1;
  TargetedExperiment exp2;

  ReactionMonitoringTransition tr1;
  ReactionMonitoringTransition tr2;
  ReactionMonitoringTransition tr3;

  tr1.setPrecursorMZ(100.0);
  tr2.setPrecursorMZ(200.0);
  tr3.setPrecursorMZ(300.0);

  std::vector<ReactionMonitoringTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.setTransitions(transitions);

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenSwath::LightTargetedExperiment &targeted_exp, OpenSwath::LightTargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  LightTargetedExperiment exp1;
  LightTargetedExperiment exp2;

  LightTransition tr1;
  LightTransition tr2;
  LightTransition tr3;

  tr1.precursor_mz = 100.0;
  tr2.precursor_mz = 200.0;
  tr3.precursor_mz = 300.0;

  std::vector<LightTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.transitions = transitions;

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION( (template < class TargetedExperimentT > static bool checkSwathMapAndSelectTransitions(const OpenMS::PeakMap &exp, const TargetedExperimentT &targeted_exp, TargetedExperimentT &transition_exp_used, double min_upper_edge_dist)))
{
  // tested above already
  NOT_TESTABLE
}
END_SECTION

START_SECTION(static void checkSwathMap(const OpenMS::PeakMap &swath_map, double &lower, double &upper))
{
  OpenMS::PeakMap swath_map;
  OpenMS::MSSpectrum spectrum;
  OpenMS::Precursor prec;
  std::vector<Precursor> precursors;
  prec.setMZ(250);
  prec.setIsolationWindowLowerOffset(50);
  prec.setIsolationWindowUpperOffset(50);
  precursors.push_back(prec);
  spectrum.setPrecursors(precursors);
  swath_map.addSpectrum(spectrum);

  double lower, upper, center;
  OpenSwathHelper::checkSwathMap(swath_map, lower, upper, center);

  TEST_REAL_SIMILAR(lower, 200);
  TEST_REAL_SIMILAR(upper, 300);
  TEST_REAL_SIMILAR(center, 250);
}
END_SECTION

START_SECTION((static std::pair<double,double> estimateRTRange(OpenSwath::LightTargetedExperiment & exp)))
{
  LightTargetedExperiment exp;

  LightCompound pep1;
  LightCompound pep2;
  LightCompound pep3;

  pep1.rt = -100.0;
  pep2.rt = 900.0;
  pep3.rt = 300.0;

  std::vector<LightCompound> peptides;
  peptides.push_back(pep1);
  peptides.push_back(pep2);
  peptides.push_back(pep3);

  exp.compounds = peptides;

  std::pair<double, double> range = OpenSwathHelper::estimateRTRange(exp);
  TEST_REAL_SIMILAR(range.first, -100)
  TEST_REAL_SIMILAR(range.second, 900)
}
END_SECTION

START_SECTION((static std::map<std::string, double> simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
        bool useQualCutoff = false, double qualCutoff = 0.0)))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


