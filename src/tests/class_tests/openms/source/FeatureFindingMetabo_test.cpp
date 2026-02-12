// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Mohammed Alhigaylan$
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

// Test with ion mobility data
START_SECTION(([EXTRA] run with ion mobility data))
{
  // Load LC-IMS-MS data
  PeakMap im_input;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("../../../topp/FeatureFinderMetabo_6_input.mzML"), im_input);

  // Run MassTraceDetection with ion mobility tolerance
  MassTraceDetection mtd_im;
  Param mtd_p = mtd_im.getParameters();
  mtd_p.setValue("mass_error_ppm", 10.0);
  mtd_p.setValue("noise_threshold_int", 10.0);
  mtd_p.setValue("ion_mobility_tolerance", 0.01);
  mtd_im.setParameters(mtd_p);
  std::vector<MassTrace> im_traces;
  mtd_im.run(im_input, im_traces);
  TEST_EQUAL(im_traces.empty(), false);

  // Verify mass traces contain IM data
  bool has_im = false;
  for (const auto& mt : im_traces)
  {
    if (mt.containsIMData())
    {
      has_im = true;
      break;
    }
  }
  TEST_EQUAL(has_im, true);

  // Run ElutionPeakDetection
  ElutionPeakDetection epd_im;
  std::vector<MassTrace> im_splitted;
  epd_im.detectPeaks(im_traces, im_splitted);
  TEST_EQUAL(im_splitted.empty(), false);

  // Verify IM data survives elution peak detection splitting
  has_im = false;
  for (const auto& mt : im_splitted)
  {
    if (mt.containsIMData())
    {
      has_im = true;
      break;
    }
  }
  TEST_EQUAL(has_im, true);

  // Run FeatureFindingMetabo with IM-aware settings
  FeatureFindingMetabo ffm_im;
  Param ffm_p = ffm_im.getParameters();
  ffm_p.setValue("mz_scoring_13C", "true");
  ffm_im.setParameters(ffm_p);
  FeatureMap im_fm;
  std::vector<std::vector<OpenMS::MSChromatogram>> im_chromatograms;
  ffm_im.run(im_splitted, im_fm, im_chromatograms);
  TEST_EQUAL(im_fm.empty(), false);

  // Verify that output features contain ion mobility metadata
  bool has_im_meta = false;
  for (const auto& f : im_fm)
  {
    if (f.metaValueExists("masstrace_centroid_im"))
    {
      has_im_meta = true;
      break;
    }
  }
  TEST_EQUAL(has_im_meta, true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



