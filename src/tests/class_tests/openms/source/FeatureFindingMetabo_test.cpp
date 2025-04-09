// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFindingMetabo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFindingMetabo* ptr = nullptr;
FeatureFindingMetabo* null_ptr = nullptr;
START_SECTION(FeatureFindingMetabo())
{
  ptr = new FeatureFindingMetabo();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFindingMetabo())
{
  delete ptr;
}
END_SECTION

// load a mzML file for testing the algorithm
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_input1.mzML"), input);

FeatureMap test_fm;

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

std::vector<std::vector< OpenMS::MSChromatogram > > chromatograms;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
test_epd.detectPeaks(output_mt, splitted_mt);

FuzzyStringComparator fsc;
fsc.setAcceptableRelative(1.001);
fsc.setAcceptableAbsolute(1);
StringList sl;
sl.push_back("xml-stylesheet");
sl.push_back("<featureMap");
sl.push_back("<feature id");
fsc.setWhitelist(sl);

//std::cout << "\n\n" << fsc.compareStrings("529090", "529091") << "\n\n\n";

START_SECTION((void run(std::vector< MassTrace > &, FeatureMap &, chromatograms &)))
{
  FeatureFindingMetabo test_ffm;
  // run with non-default setting (C13 isotope distance)
  Param p = test_ffm.getParameters();
  p.setValue("mz_scoring_13C", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 84);

  // run with default settings (from paper using charge+isotope# dependent distances)
  p.setValue("report_convex_hulls", "true");
  p.setValue("mz_scoring_13C", "false");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 81);
  // --> this gives less features, i.e. more isotope clusters (but the input data is simulated and highly weird -- should be replaced at some point)

  // test annotation of input
  String tmp_file;
  NEW_TMP_FILE(tmp_file);
  FeatureXMLFile().store(tmp_file, test_fm);
  TEST_EQUAL(fsc.compareFiles(tmp_file, OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_output1.featureXML")), true);

  //todo the new isotope m/z scoring should produce similar results, but still has to be tested.
  p.setValue("report_convex_hulls", "true");
  p.setValue("mz_scoring_by_elements", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 80);
  // --> this gives less features, i.e. more isotope clusters (but the input data is simulated and highly weird -- should be replaced at some point)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



