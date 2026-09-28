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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <cmath>
#include <memory>

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

START_SECTION((void setParameters(const Param&) / CalibrationResult performCalibration(...)))
{
  auto make_targeted_experiment = []()
  {
    OpenSwath::LightTargetedExperiment targeted_exp;

    OpenSwath::LightCompound comp1;
    comp1.id = "pep1";
    comp1.rt = 100.0;
    targeted_exp.getCompounds().push_back(comp1);

    OpenSwath::LightCompound comp2;
    comp2.id = "pep2";
    comp2.rt = 200.0;
    targeted_exp.getCompounds().push_back(comp2);

    OpenSwath::LightCompound comp3;
    comp3.id = "pep3";
    comp3.rt = 300.0;
    targeted_exp.getCompounds().push_back(comp3);

    OpenSwath::LightTransition tr1;
    tr1.transition_name = "t1";
    tr1.peptide_ref = "pep1";
    tr1.product_mz = 500.0;
    tr1.precursor_mz = 400.0;
    targeted_exp.getTransitions().push_back(tr1);

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

    return targeted_exp;
  };

  auto make_irt_chromatograms = []()
  {
    std::vector<MSChromatogram> chromatograms;

    auto add_chrom = [&](const std::string& native_id, double apex_rt, double precursor_mz, double product_mz)
    {
      MSChromatogram chrom;
      chrom.setNativeID(native_id);

      Precursor precursor;
      precursor.setMZ(precursor_mz);
      chrom.setPrecursor(precursor);

      Product product;
      product.setMZ(product_mz);
      chrom.setProduct(product);

      for (double rt = apex_rt - 2.0; rt <= apex_rt + 2.0; rt += 0.1)
      {
        ChromatogramPeak peak;
        peak.setRT(rt);
        const double dist = std::fabs(rt - apex_rt);
        double intensity = 200.0;
        if (dist < 0.0001) intensity = 2200.0;
        else if (dist < 0.6) intensity = 1400.0;
        peak.setIntensity(intensity);
        chrom.push_back(peak);
      }

      chromatograms.push_back(chrom);
    };

    // Keep one peptide slightly off the global linear trend so estimateWindow() yields
    // a small but positive RT estimate that can be clamped in both directions.
    add_chrom("t1", 100.0, 400.0, 500.0);
    add_chrom("t2", 201.0, 500.0, 600.0);
    add_chrom("t3", 300.0, 600.0, 700.0);

    return chromatograms;
  };

  auto make_swath_maps = [](const std::vector<MSChromatogram>& chromatograms)
  {
    auto peak_map = std::make_shared<PeakMap>();
    peak_map->setChromatograms(chromatograms);

    OpenSwath::SwathMap swath_map;
    swath_map.sptr = std::make_shared<SpectrumAccessOpenMS>(peak_map);
    swath_map.ms1 = false;

    return std::vector<OpenSwath::SwathMap>{swath_map};
  };

  auto make_irt_detection_params = []()
  {
    Param p;
    p.setValue("outlierMethod", "none");
    p.setValue("alignmentMethod", "linear");
    p.setValue("estimateBestPeptides", "false");
    p.setValue("OverallQualityCutoff", 0.0);
    p.setValue("NrRTBins", 10);
    p.setValue("MinPeptidesPerBin", 1);
    p.setValue("MinBinsFilled", 1);
    p.setValue("lowess:span", 0.75);
    p.setValue("lowess:auto_span", "false");
    p.setValue("lowess:auto_span_min", 0.2);
    p.setValue("lowess:auto_span_max", 2.0);
    p.setValue("lowess:auto_span_grid", 5);
    p.setValue("b_spline:num_nodes", 20);
    return p;
  };

  auto make_calibration_params = []()
  {
    Param p;
    p.setValue("mz_extraction_window", 10.0);
    p.setValue("mz_extraction_window_ppm", "false");
    return p;
  };

  auto make_chrom_extract_params = []()
  {
    ChromExtractParams cp;
    cp.min_upper_edge_dist = 0.0;
    cp.mz_extraction_window = 10.0;
    cp.im_extraction_window = 0.0;
    cp.ppm = false;
    cp.extraction_function = "none";
    cp.rt_extraction_window = 30.0;
    cp.extra_rt_extract = 0.0;
    return cp;
  };

  auto run_calibration = [&](double min_rt_window, double max_rt_window, double rt_padding_factor)
  {
    CalibrationWorkflow cw;
    Param workflow_params = cw.getParameters();
    workflow_params.setValue("linear:outlier_detection", "none");
    workflow_params.setValue("nonlinear:outlier_detection", "none");
    workflow_params.setValue("windows:estimate_rt", "true");
    workflow_params.setValue("windows:estimate_mz", "false");
    workflow_params.setValue("windows:estimate_im", "false");
    workflow_params.setValue("windows:rt_estimation_padding_factor", rt_padding_factor);
    workflow_params.setValue("windows:min_rt_window", min_rt_window);
    workflow_params.setValue("windows:max_rt_window", max_rt_window);
    cw.setParameters(workflow_params);

    auto targeted_exp = make_targeted_experiment();
    auto chromatograms = make_irt_chromatograms();
    auto swath_maps = make_swath_maps(chromatograms);

    Param feature_finder_param;
    feature_finder_param.setValue("use_ms1_ion_mobility", "false");

    ChromExtractParams cp = make_chrom_extract_params();
    ChromExtractParams cp_ms1 = cp;
    ChromExtractParams cp_irt = cp;

    CalibrationWorkflow::IrtExperiments irt_experiments;
    irt_experiments.linear_irt = targeted_exp;
    irt_experiments.is_prepared = true;

    const CalibrationWorkflow::CalibrationResult result = cw.performCalibration(
      swath_maps,
      targeted_exp,
      cp,
      cp_ms1,
      irt_experiments,
      feature_finder_param,
      cp_irt,
      make_irt_detection_params(),
      make_calibration_params(),
      Param(),
      false,
      false,
      "",
      "",
      0);

    return std::make_pair(cp.rt_extraction_window, result.estimated_rt_window);
  };

  const auto [floor_applied, floor_raw] = run_calibration(60.0, 600.0, 1.0);
  TEST_TRUE(floor_raw > 1e-9)
  TEST_TRUE(floor_raw < 60.0)
  TEST_TRUE(std::fabs(floor_applied - 60.0) < 1e-6)

  const auto [ceiling_applied, ceiling_raw] = run_calibration(0.0, 600.0, 10000.0);
  TEST_TRUE(ceiling_raw > 600.0)
  TEST_TRUE(std::fabs(ceiling_applied - 600.0) < 1e-6)

  const auto [unclamped_applied, unclamped_raw] = run_calibration(0.0, 0.0, 1.0);
  TEST_TRUE(unclamped_raw > 1e-9)
  TEST_TRUE(std::fabs(unclamped_applied - unclamped_raw) < 1e-6)

  const auto [equal_bound_applied, equal_bound_raw] = run_calibration(120.0, 120.0, 10000.0);
  TEST_TRUE(equal_bound_raw > 120.0)
  TEST_TRUE(std::fabs(equal_bound_applied - 120.0) < 1e-6)

  CalibrationWorkflow invalid_min_cw;
  Param invalid_min_params = invalid_min_cw.getParameters();
  invalid_min_params.setValue("windows:min_rt_window", 1e-10);
  TEST_EXCEPTION(Exception::InvalidParameter, invalid_min_cw.setParameters(invalid_min_params))

  CalibrationWorkflow invalid_max_cw;
  Param invalid_max_params = invalid_max_cw.getParameters();
  invalid_max_params.setValue("windows:max_rt_window", 1e-10);
  TEST_EXCEPTION(Exception::InvalidParameter, invalid_max_cw.setParameters(invalid_max_params))

  CalibrationWorkflow inverted_bounds_cw;
  Param inverted_bounds_params = inverted_bounds_cw.getParameters();
  inverted_bounds_params.setValue("windows:min_rt_window", 100.0);
  inverted_bounds_params.setValue("windows:max_rt_window", 60.0);
  TEST_EXCEPTION(Exception::InvalidParameter, inverted_bounds_cw.setParameters(inverted_bounds_params))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((CalibrationWorkflow::CalibrationResult performCalibration(std::vector<OpenSwath::SwathMap>&, OpenSwath::LightTargetedExperiment&, OpenSwath::ChromExtractParams&, OpenSwath::ChromExtractParams&, const CalibrationWorkflow::IrtExperiments&, const Param&, const OpenSwath::ChromExtractParams&, const Param&, const Param&, const Param&, bool, bool, const std::string&, const std::string&, Size)))
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
