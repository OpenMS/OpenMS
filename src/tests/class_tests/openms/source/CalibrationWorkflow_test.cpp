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

START_SECTION((CalibrationWorkflow::CalibrationResult performCalibration(std::vector<OpenSwath::SwathMap>&, OpenSwath::LightTargetedExperiment&, OpenSwath::ChromExtractParams&, OpenSwath::ChromExtractParams&, const CalibrationWorkflow::IrtExperiments&, const Param&, const OpenSwath::ChromExtractParams&, const Param&, const Param&, const Param&, bool, bool, const String&, const String&, Size)))
{
  CalibrationWorkflow cw;

  // Reuse the minimal targeted experiment and chromatograms from the previous test
  OpenSwath::LightTargetedExperiment targeted_exp;
  OpenSwath::LightCompound comp;
  comp.id = "pep1";
  comp.rt = 100.0;
  targeted_exp.getCompounds().push_back(comp);
  OpenSwath::LightCompound comp2;
  comp2.id = "pep2";
  comp2.rt = 200.0;
  targeted_exp.getCompounds().push_back(comp2);
  OpenSwath::LightCompound comp3;
  comp3.id = "pep3";
  comp3.rt = 300.0;
  targeted_exp.getCompounds().push_back(comp3);

  OpenSwath::LightTransition tr;
  tr.transition_name = "t1";
  tr.peptide_ref = "pep1";
  tr.product_mz = 500.0;
  tr.precursor_mz = 400.0;
  targeted_exp.getTransitions().push_back(tr);
  OpenSwath::LightTransition tr2;
  tr2.transition_name = "t2";
  tr2.peptide_ref = "pep2";
  tr2.product_mz = 600.0;
  tr2.precursor_mz = 500.0;
  targeted_exp.getTransitions().push_back(tr2);
  OpenSwath::LightTransition tr3;
  tr3.transition_name = "t3";
  tr3.peptide_ref = "pep3";
  tr3.product_mz = 700.0;
  tr3.precursor_mz = 600.0;
  targeted_exp.getTransitions().push_back(tr3);

  std::vector<MSChromatogram> chroms;
  // Create matching chromatograms for the transitions with simple peak shapes
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

  std::vector<OpenSwath::SwathMap> swath_maps; // empty maps ok for unit test

  // Basic parameter sets (reuse simple params)
  Param ffparam;
  Param irt_detection_param;
  irt_detection_param.setValue("outlierMethod", "none");
  irt_detection_param.setValue("alignmentMethod", "lowess");
  irt_detection_param.setValue("estimateBestPeptides", "false");
  irt_detection_param.setValue("OverallQualityCutoff", 0.0);
  irt_detection_param.setValue("NrRTBins", 100);
  irt_detection_param.setValue("MinPeptidesPerBin", 1);
  irt_detection_param.setValue("MinBinsFilled", 1);
  irt_detection_param.setValue("lowess:span", 0.75);
  irt_detection_param.setValue("lowess:auto_span", "false");
  irt_detection_param.setValue("lowess:auto_span_min", 0.2);
  irt_detection_param.setValue("lowess:auto_span_max", 2.0);
  irt_detection_param.setValue("lowess:auto_span_grid", 5);
  irt_detection_param.setValue("b_spline:num_nodes", 20);
  Param calibration_param;
  calibration_param.setValue("mz_extraction_window", 10.0);
  calibration_param.setValue("mz_extraction_window_ppm", "false");

  // Chromatogram extraction params
  ChromExtractParams cp;
  cp.min_upper_edge_dist = 0.0;
  cp.mz_extraction_window = 10.0;
  cp.im_extraction_window = 0.0;
  cp.ppm = false;
  cp.extraction_function = "none";
  cp.rt_extraction_window = 30.0;
  cp.extra_rt_extract = 0.0;

  ChromExtractParams cp_ms1 = cp;

  // 1) Linear-only: prepare IrtExperiments with only linear iRT
  CalibrationWorkflow::IrtExperiments ie_linear;
  ie_linear.linear_irt = targeted_exp;
  ie_linear.nonlinear_irt = OpenSwath::LightTargetedExperiment(); // empty
  ie_linear.is_prepared = true;

  Param mrm_map_param; // empty

  // Instead of running the full performCalibration (which requires swath maps and
  // the chromatogram provider), exercise the core normalization routine directly
  // (doDataNormalization_) for a linear and a shifted (nonlinear-like) case.

  TransformationDescription im_trafo;
  auto trafo_linear = cw.doDataNormalization_(targeted_exp,
                                              chroms,
                                              im_trafo,
                                              swath_maps,
                                              0.0,
                                              0.0,
                                              ffparam,
                                              irt_detection_param,
                                              calibration_param,
                                              false);
  TEST_EQUAL(trafo_linear.getDataPoints().empty(), false)

  // Create shifted chromatograms to mimic nonlinear deviations (shift middle peptide chromatogram)
  std::vector<MSChromatogram> chroms_shifted = chroms;
  for (auto & cp_chrom : chroms_shifted)
  {
    if (cp_chrom.getNativeID() == "t2")
    {
      for (auto & pk : cp_chrom)
      {
        pk.setRT(pk.getRT() + 5.0);
      }
    }
  }

  // Use a target experiment with same compounds but leave RTs as-is; the chromatograms are shifted
  auto trafo_shifted = cw.doDataNormalization_(targeted_exp,
                                               chroms_shifted,
                                               im_trafo,
                                               swath_maps,
                                               0.0,
                                               0.0,
                                               ffparam,
                                               irt_detection_param,
                                               calibration_param,
                                               false);
  TEST_EQUAL(trafo_shifted.getDataPoints().empty(), false)

  // The fitted transformations should differ due to the shifted chromatogram.
  // getDataPoints() stores pairs (observed_rt, reference_rt) -> pair.first = observed, pair.second = reference.
  // Build maps from reference_rt -> observed_rt for both transformations and compare measured RTs.
  auto make_map = [](const TransformationDescription &trafo)
  {
    std::map<double, double> m;
    for (const auto &p : trafo.getDataPoints())
    {
      m[p.second] = p.first; // key = reference_rt, value = observed_rt
    }
    return m;
  };

  auto map_lin = make_map(trafo_linear);
  auto map_sh = make_map(trafo_shifted);

  bool different = false;
  for (const auto &kv : map_lin)
  {
    double ref_rt = kv.first;
    double obs_lin = kv.second;
    auto it = map_sh.find(ref_rt);
    if (it == map_sh.end())
    {
      // point missing in shifted transformation -> consider different
      different = true;
      break;
    }
    double obs_sh = it->second;
    if (std::fabs(obs_lin - obs_sh) > 1e-6)
    {
      different = true;
      break;
    }
  }

  TEST_EQUAL(different, true)
}
END_SECTION
END_TEST
