// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(CalibrationWorkflow, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((TransformationDescription doDataNormalization_(const OpenSwath::LightTargetedExperiment&, const std::vector< OpenMS::MSChromatogram >&, TransformationDescription&, std::vector< OpenSwath::SwathMap >&, double, double, const Param&, const Param&, const Param&, const bool)))
{
  CalibrationWorkflow cw;

  // Build a minimal LightTargetedExperiment with one peptide and one transition
  OpenSwath::LightTargetedExperiment targeted_exp;
  OpenSwath::LightCompound comp;
  comp.id = "pep1";
  comp.rt = 100.0;
  targeted_exp.getCompounds().push_back(comp);
  // Add a second iRT peptide so regression can be performed
  OpenSwath::LightCompound comp2;
  comp2.id = "pep2";
  comp2.rt = 200.0;
  targeted_exp.getCompounds().push_back(comp2);

  // Add a third peptide to ensure model fitting has >= 3 points
  OpenSwath::LightCompound comp3;
  comp3.id = "pep3";
  comp3.rt = 300.0;
  targeted_exp.getCompounds().push_back(comp3);

  OpenSwath::LightTransition tr;
  tr.transition_name = "t1"; // nativeID
  tr.peptide_ref = "pep1";
  tr.product_mz = 500.0;
  tr.precursor_mz = 400.0;
  targeted_exp.getTransitions().push_back(tr);

  OpenSwath::LightTransition tr2;
  tr2.transition_name = "t2"; // nativeID
  tr2.peptide_ref = "pep2";
  tr2.product_mz = 600.0;
  tr2.precursor_mz = 500.0;
  targeted_exp.getTransitions().push_back(tr2);

  OpenSwath::LightTransition tr3;
  tr3.transition_name = "t3"; // nativeID
  tr3.peptide_ref = "pep3";
  tr3.product_mz = 700.0;
  tr3.precursor_mz = 600.0;
  targeted_exp.getTransitions().push_back(tr3);

  // Create matching chromatograms for the transitions with simple peak shapes
  std::vector<MSChromatogram> chroms;

  MSChromatogram chrom;
  chrom.setNativeID("t1");
  // create a small triangular peak around RT=100 with multiple points (denser sampling)
  for (double rt = 98.0; rt <= 102.0; rt += 0.1)
  {
    ChromatogramPeak cp;
    cp.setRT(rt);
    double intensity = 0.0;
    double dist = fabs(rt - 100.0);
    if (dist < 0.0001) intensity = 2000.0;
    else if (dist < 0.6) intensity = 1200.0;
    else intensity = 150.0;
    cp.setIntensity(intensity);
    chrom.push_back(cp);
  }
  chroms.push_back(chrom);

  MSChromatogram chrom2;
  chrom2.setNativeID("t2");
  for (double rt = 198.0; rt <= 202.0; rt += 0.1)
  {
    ChromatogramPeak cp;
    cp.setRT(rt);
    double intensity = 0.0;
    double dist = fabs(rt - 200.0);
    if (dist < 0.0001) intensity = 2100.0;
    else if (dist < 0.6) intensity = 1300.0;
    else intensity = 180.0;
    cp.setIntensity(intensity);
    chrom2.push_back(cp);
  }
  chroms.push_back(chrom2);

  // third chromatogram
  MSChromatogram chrom3;
  chrom3.setNativeID("t3");
  for (double rt = 298.0; rt <= 302.0; rt += 0.1)
  {
    ChromatogramPeak cp;
    cp.setRT(rt);
    double dist = fabs(rt - 300.0);
    double intensity = 0.0;
    if (dist < 0.0001) intensity = 2200.0;
    else if (dist < 0.6) intensity = 1400.0;
    else intensity = 200.0;
    cp.setIntensity(intensity);
    chrom3.push_back(cp);
  }
  chroms.push_back(chrom3);

  TransformationDescription im_trafo;
  std::vector<OpenSwath::SwathMap> swath_maps; // empty is fine for this unit test

  // Prepare parameters: use conservative settings to avoid outlier rejection
  Param ffparam; // default feature finder params (will be populated by doDataNormalization_)
  Param irt_detection_param;
  // Populate required detection parameters (match defaults used in CalibrationWorkflow)
  irt_detection_param.setValue("outlierMethod", "none");
  irt_detection_param.setValue("alignmentMethod", "lowess");
  irt_detection_param.setValue("estimateBestPeptides", "false");
  irt_detection_param.setValue("OverallQualityCutoff", 0.0);
  irt_detection_param.setValue("NrRTBins", 100);
  irt_detection_param.setValue("MinPeptidesPerBin", 1);
  irt_detection_param.setValue("MinBinsFilled", 1);
  // Lowess / bspline model params
  irt_detection_param.setValue("lowess:span", 0.75);
  irt_detection_param.setValue("lowess:auto_span", "false");
  irt_detection_param.setValue("lowess:auto_span_min", 0.2);
  irt_detection_param.setValue("lowess:auto_span_max", 2.0);
  irt_detection_param.setValue("lowess:auto_span_grid", 5);
  irt_detection_param.setValue("b_spline:num_nodes", 20);

  Param calibration_param;
  calibration_param.setValue("mz_extraction_window", 10.0);
  calibration_param.setValue("mz_extraction_window_ppm", "false");

  // Call the routine under test
  TransformationDescription trafo_out = cw.doDataNormalization_(targeted_exp,
                                                               chroms,
                                                               im_trafo,
                                                               swath_maps,
                                                               0.0, // min_rsq
                                                               0.0, // min_coverage
                                                               ffparam,
                                                               irt_detection_param,
                                                               calibration_param,
                                                               false);

  // Expect at least one data point (we provided one chromatogram matching the peptide)
  TEST_EQUAL(trafo_out.getDataPoints().size() > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
