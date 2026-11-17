// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIACalibration.h"

#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace OpenMS;

namespace
{
  void require(bool value, const std::string& message)
  {
    if (!value)
    {
      throw std::runtime_error("OpenDIACalibration test failed: " + message);
    }
  }

  // Known truth: run RT = 300 + 5 x library RT + 0.02 x library RT^2 (s), fragments +4 ppm, 1/K0 +0.02.
  double trueRT(double library_rt) { return 300.0 + 5.0 * library_rt + 0.02 * library_rt * library_rt; }
  constexpr double PPM_SHIFT = 4.0;
  constexpr double IM_SHIFT = 0.02;
}

int main()
{
  const std::size_t pairs = 300;
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> noise(0.0, 1.0);

  OpenSwath::LightTargetedExperiment library;
  std::vector<std::pair<std::string, std::string>> traml_ids;
  struct Truth
  {
    double library_rt, run_rt, precursor_mz, im;
    std::vector<double> fragments, intensities;
    bool present;
  };
  std::vector<Truth> truth;
  for (std::size_t p = 0; p < pairs; ++p)
  {
    Truth t;
    t.library_rt = 100.0 * static_cast<double>(p) / pairs;
    t.run_rt = trueRT(t.library_rt) + 2.0 * noise(rng);
    t.precursor_mz = 450.0 + 100.0 * uniform(rng);
    t.im = 0.8 + 0.5 * uniform(rng);
    t.present = p % 5 != 0; // 20% of targets are absent.
    for (std::size_t f = 0; f < 6; ++f)
    {
      t.fragments.push_back(200.0 + 700.0 * uniform(rng));
      t.intensities.push_back(100.0 * (6.0 - static_cast<double>(f)));
    }
    for (int decoy = 0; decoy < 2; ++decoy)
    {
      const std::string id = std::to_string(2 * p + decoy);
      OpenSwath::LightCompound compound;
      compound.id = id;
      compound.rt = t.library_rt;
      compound.drift_time = t.im;
      library.compounds.push_back(compound);
      traml_ids.emplace_back(id, (decoy ? "DECOY_pep" : "pep") + std::to_string(p));
      for (std::size_t f = 0; f < 6; ++f)
      {
        OpenSwath::LightTransition transition;
        transition.peptide_ref = id;
        transition.precursor_mz = t.precursor_mz;
        // Decoy fragments are shifted and do not co-elute.
        transition.product_mz = t.fragments[f] + (decoy ? 7.3 : 0.0);
        transition.library_intensity = t.intensities[f];
        transition.setDecoy(decoy == 1);
        library.transitions.push_back(transition);
      }
    }
    truth.push_back(std::move(t));
  }

  OpenDIA::CalibrationOptions options;
  options.panel_targets = pairs;
  OpenDIA::CalibrationCollector collector(library, traml_ids, options);
  require(collector.panelSize() == 2 * pairs, "panel size");

  // Spectra every 1 s over the run; one isolation window covering all precursors.
  const double run_start = 250.0, run_end = 300.0 + 500.0 + 200.0 + 50.0;
  for (double rt = run_start; rt <= run_end; rt += 1.0)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setRT(rt);
    Precursor precursor;
    precursor.setMZ(500.0);
    precursor.setIsolationWindowLowerOffset(60.0);
    precursor.setIsolationWindowUpperOffset(60.0);
    spectrum.setPrecursors({precursor});
    DataArrays::FloatDataArray mobility;
    IMDataConverter::setIMUnit(mobility, DriftTimeUnit::VSSC);
    for (const auto& t : truth)
    {
      if (!t.present)
      {
        continue;
      }
      const double shape = std::exp(-0.5 * std::pow((rt - t.run_rt) / 4.0, 2));
      if (shape < 1e-3)
      {
        continue;
      }
      for (std::size_t f = 0; f < 6; ++f)
      {
        const double mz = t.fragments[f] * (1.0 + (PPM_SHIFT + noise(rng)) * 1e-6);
        spectrum.emplace_back(mz, static_cast<float>(t.intensities[f] * shape * (1.0 + 0.05 * noise(rng))));
        mobility.push_back(static_cast<float>(t.im + IM_SHIFT + 0.003 * noise(rng)));
      }
    }
    // Random interference.
    for (int k = 0; k < 400; ++k)
    {
      spectrum.emplace_back(200.0 + 710.0 * uniform(rng), static_cast<float>(50.0 * uniform(rng)));
      mobility.push_back(static_cast<float>(0.8 + 0.5 * uniform(rng)));
    }
    spectrum.getFloatDataArrays().push_back(std::move(mobility));
    collector.consumeSpectrum(spectrum);
  }

  const auto result = collector.calibrate();
  require(result.em_converged, "EM converged");
  require(result.anchor_targets >= 150, "anchor count " + std::to_string(result.anchor_targets));
  require(result.anchor_decoys <= 3, "decoy anchors " + std::to_string(result.anchor_decoys));
  for (double library_rt : {10.0, 50.0, 90.0})
  {
    require(std::fabs(result.runRT(library_rt) - trueRT(library_rt)) < 5.0,
            "RT map at " + std::to_string(library_rt) + ": " + std::to_string(result.runRT(library_rt)));
  }
  require(result.rt_map.sigma > 0.5 && result.rt_map.sigma < 5.0, "RT sigma " + std::to_string(result.rt_map.sigma));
  require(std::fabs(result.ms2_offset_ppm - PPM_SHIFT) < 0.5, "ppm offset " + std::to_string(result.ms2_offset_ppm));
  require(result.ms2_sigma_ppm > 0.3 && result.ms2_sigma_ppm < 3.0, "ppm sigma " + std::to_string(result.ms2_sigma_ppm));
  // Fragments scatter N(0, 1 ppm) over an almost empty background: the fitted half width keeps the whole signal
  // (beyond 3 sigma) and stays inside the search window about the offset.
  require(result.ms2_tolerance_ppm > 3.0 && result.ms2_tolerance_ppm <= options.search_ppm - std::fabs(result.ms2_offset_ppm) + 1e-9,
          "fitted ppm tolerance " + std::to_string(result.ms2_tolerance_ppm));
  require(result.ms2_fit.background_share < 0.2, "fitted background share " + std::to_string(result.ms2_fit.background_share));
  require(result.has_im, "IM model");
  require(std::fabs(result.runIM(1.0) - (1.0 + IM_SHIFT)) < 0.005, "IM offset " + std::to_string(result.runIM(1.0)));
  require(result.im_sigma < 0.01, "IM sigma " + std::to_string(result.im_sigma));
  for (const auto& anchor : result.anchors)
  {
    require(std::stoul(anchor.precursor_id) % 2 == 0, "decoy reported as anchor");
    require(std::stoul(anchor.precursor_id) / 2 % 5 != 0, "absent target reported as anchor");
  }
  std::ostringstream json, anchors;
  OpenDIA::writeCalibrationJSON(result, options, json);
  OpenDIA::writeCalibrationAnchors(result, anchors);
  require(json.str().find("\"offset_ppm\"") != std::string::npos, "JSON output");

  // Diagnostics do not change the result and record every panel member's stage values.
  OpenDIA::CalibrationDiagnostics diagnostics;
  const auto diagnosed = collector.calibrate(&diagnostics);
  std::ostringstream json_diagnosed;
  OpenDIA::writeCalibrationJSON(diagnosed, options, json_diagnosed);
  require(json_diagnosed.str() == json.str(), "diagnostics change the calibration");
  require(diagnostics.rows.size() == collector.panelSize() && diagnostics.failure.empty(), "diagnostic rows");
  std::size_t stage3_accepted = 0;
  for (const auto& row : diagnostics.rows)
  {
    require(std::isfinite(row.em_rt) && std::isfinite(row.stage1_map_rt), "stage maps in diagnostics");
    stage3_accepted += !row.is_decoy && row.stage3_q <= options.anchor_fdr;
  }
  require(stage3_accepted == result.anchor_targets, "diagnostic stage-3 q-values match the anchors");
  std::ostringstream diagnostic_table;
  OpenDIA::writeCalibrationDiagnostics(diagnostics, diagnostic_table);
  require(diagnostic_table.str().find("stage3_q") != std::string::npos, "diagnostic TSV");

  // Windowed second pass: collect only inside +-60 s of a supplied stage-1 map, skip stages 1 and 2.
  {
    OpenDIA::CalibrationWindow window;
    for (double library_rt = 0.0; library_rt <= 100.0; library_rt += 2.0)
    {
      window.map.library_rt.push_back(library_rt);
      window.map.run_rt.push_back(trueRT(library_rt));
    }
    window.half_window_s = 60.0;
    window.half_window_spectra = result.half_window_spectra;
    window.stage1_kept = result.stage1_kept;
    window.stage1_sigma_s = result.stage1_sigma_s;
    window.peak_fwhm_s = result.peak_fwhm_s;
    OpenDIA::CalibrationCollector windowed(library, traml_ids, options, window);
    std::mt19937 rng2(11);
    for (double rt = run_start; rt <= run_end; rt += 1.0)
    {
      MSSpectrum spectrum;
      spectrum.setMSLevel(2);
      spectrum.setRT(rt);
      Precursor precursor;
      precursor.setMZ(500.0);
      precursor.setIsolationWindowLowerOffset(60.0);
      precursor.setIsolationWindowUpperOffset(60.0);
      spectrum.setPrecursors({precursor});
      for (const auto& t : truth)
      {
        const double shape = std::exp(-0.5 * std::pow((rt - t.run_rt) / 4.0, 2));
        if (!t.present || shape < 1e-3)
        {
          continue;
        }
        for (std::size_t f = 0; f < 6; ++f)
        {
          spectrum.emplace_back(t.fragments[f] * (1.0 + PPM_SHIFT * 1e-6), static_cast<float>(t.intensities[f] * shape));
        }
      }
      for (int k = 0; k < 400; ++k)
      {
        spectrum.emplace_back(200.0 + 710.0 * uniform(rng2), static_cast<float>(50.0 * uniform(rng2)));
      }
      spectrum.sortByPosition();
      windowed.consumeSpectrum(spectrum);
    }
    OpenDIA::CalibrationDiagnostics window_diagnostics;
    const auto windowed_result = windowed.calibrate(&window_diagnostics);
    require(windowed_result.windowed, "windowed flag");
    require(windowed_result.anchor_targets >= 150, "windowed anchor count " + std::to_string(windowed_result.anchor_targets));
    require(windowed_result.anchor_decoys <= 3, "windowed decoy anchors " + std::to_string(windowed_result.anchor_decoys));
    for (double library_rt : {10.0, 50.0, 90.0})
    {
      require(std::fabs(windowed_result.runRT(library_rt) - trueRT(library_rt)) < 5.0,
              "windowed RT map at " + std::to_string(library_rt) + ": " + std::to_string(windowed_result.runRT(library_rt)) + " vs " +
                std::to_string(trueRT(library_rt)) + " (window map " + std::to_string(window.map(library_rt)) + ", anchors " +
                std::to_string(windowed_result.anchor_targets) + ")");
    }
    // Each member keeps only spectra inside its window plus margin, far fewer than the whole run.
    std::size_t max_spectra = 0;
    for (const auto& row : window_diagnostics.rows)
    {
      max_spectra = std::max(max_spectra, row.spectra);
    }
    const double collected_span = 2.0 * (window.half_window_s + 30.0 + 2.0 * window.peak_fwhm_s) + 1.0;
    require(static_cast<double>(max_spectra) <= collected_span, "windowed collection range " + std::to_string(max_spectra));
    std::ostringstream windowed_json;
    OpenDIA::writeCalibrationJSON(windowed_result, options, windowed_json);
    require(windowed_json.str().find("\"windowed_second_pass\": true") != std::string::npos, "windowed JSON marker");
    require(json.str().find("windowed_second_pass") == std::string::npos, "default JSON has no windowed marker");
  }

  // Windowed pass whose stage-1 map is supported only for library RT >= 40 and is 150 s late below it (like an
  // extrapolation past the stage-1 anchors). Members below the support are collected and searched to the run start,
  // so anchors and the final map recover the true early RT.
  {
    OpenDIA::CalibrationWindow window;
    for (double library_rt = 0.0; library_rt <= 100.0; library_rt += 2.0)
    {
      window.map.library_rt.push_back(library_rt);
      window.map.run_rt.push_back(trueRT(library_rt) + (library_rt < 40.0 ? 150.0 * (40.0 - library_rt) / 40.0 : 0.0));
    }
    window.half_window_s = 60.0;
    window.half_window_spectra = result.half_window_spectra;
    window.stage1_kept = result.stage1_kept;
    window.stage1_sigma_s = result.stage1_sigma_s;
    window.peak_fwhm_s = result.peak_fwhm_s;
    const auto run_windowed = [&](const OpenDIA::CalibrationWindow& w)
    {
      OpenDIA::CalibrationCollector windowed(library, traml_ids, options, w);
      std::mt19937 rng3(13);
      for (double rt = run_start; rt <= run_end; rt += 1.0)
      {
        MSSpectrum spectrum;
        spectrum.setMSLevel(2);
        spectrum.setRT(rt);
        Precursor precursor;
        precursor.setMZ(500.0);
        precursor.setIsolationWindowLowerOffset(60.0);
        precursor.setIsolationWindowUpperOffset(60.0);
        spectrum.setPrecursors({precursor});
        for (const auto& t : truth)
        {
          const double shape = std::exp(-0.5 * std::pow((rt - t.run_rt) / 4.0, 2));
          if (!t.present || shape < 1e-3)
          {
            continue;
          }
          for (std::size_t f = 0; f < 6; ++f)
          {
            spectrum.emplace_back(t.fragments[f] * (1.0 + PPM_SHIFT * 1e-6), static_cast<float>(t.intensities[f] * shape));
          }
        }
        for (int k = 0; k < 400; ++k)
        {
          spectrum.emplace_back(200.0 + 710.0 * uniform(rng3), static_cast<float>(50.0 * uniform(rng3)));
        }
        spectrum.sortByPosition();
        windowed.consumeSpectrum(spectrum);
      }
      return windowed.calibrate();
    };
    // Without a support range the late early-gradient centres hide those precursors (the old behaviour).
    const auto unsupported = run_windowed(window);
    std::size_t early_unsupported = 0;
    for (const auto& anchor : unsupported.anchors)
    {
      early_unsupported += anchor.used && anchor.library_rt < 20.0;
    }
    window.support_low = 40.0;
    const auto supported = run_windowed(window);
    std::size_t early_supported = 0;
    for (const auto& anchor : supported.anchors)
    {
      early_supported += anchor.used && anchor.library_rt < 20.0;
    }
    require(supported.stage3_extended_targets > 0, "no target extended beyond the support");
    require(early_supported >= 30 && early_supported > early_unsupported + 20,
            "early anchors with support " + std::to_string(early_supported) + " vs without " + std::to_string(early_unsupported));
    require(supported.anchor_decoys <= 3, "decoy anchors with extended windows " + std::to_string(supported.anchor_decoys));
    require(std::fabs(supported.runRT(5.0) - trueRT(5.0)) < 5.0,
            "early RT map with support " + std::to_string(supported.runRT(5.0)) + " vs " + std::to_string(trueRT(5.0)));
  }

  // One-pass calibration on a gradient that bends at the start (run RT rises 1.6x faster below library RT 20; HeLa's
  // stage-1 anchors bend by about this much, PRECHECK_RTCENTRE). Stage 3 is centred on the stage-1 map, which follows
  // its anchors, so the final map is right at the start (the smooth stage-2 EM map must not set the window).
  {
    const auto bentRT = [](double library_rt) { return library_rt < 20.0 ? 200.0 + 8.0 * library_rt : 360.0 + 5.0 * (library_rt - 20.0); };
    std::vector<Truth> bent = truth;
    std::mt19937 rng4(17);
    for (auto& t : bent) t.run_rt = bentRT(t.library_rt) + 2.0 * noise(rng4);
    OpenDIA::CalibrationCollector one_pass(library, traml_ids, options);
    for (double rt = 50.0; rt <= 1050.0; rt += 1.0)
    {
      MSSpectrum spectrum;
      spectrum.setMSLevel(2);
      spectrum.setRT(rt);
      Precursor precursor;
      precursor.setMZ(500.0);
      precursor.setIsolationWindowLowerOffset(60.0);
      precursor.setIsolationWindowUpperOffset(60.0);
      spectrum.setPrecursors({precursor});
      for (const auto& t : bent)
      {
        const double shape = std::exp(-0.5 * std::pow((rt - t.run_rt) / 4.0, 2));
        if (!t.present || shape < 1e-3) continue;
        for (std::size_t f = 0; f < 6; ++f) spectrum.emplace_back(t.fragments[f] * (1.0 + PPM_SHIFT * 1e-6), static_cast<float>(t.intensities[f] * shape));
      }
      for (int k = 0; k < 400; ++k) spectrum.emplace_back(200.0 + 710.0 * uniform(rng4), static_cast<float>(50.0 * uniform(rng4)));
      spectrum.sortByPosition();
      one_pass.consumeSpectrum(spectrum);
    }
    OpenDIA::CalibrationDiagnostics bent_diagnostics;
    const auto bent_result = one_pass.calibrate(&bent_diagnostics);
    require(!bent_result.windowed, "bent gradient calibrates in one pass");
    // The centre of every stage-3 window is the stage-1 map (the EM map is reported, not used for the window).
    for (const auto& row : bent_diagnostics.rows)
    {
      require(row.stage3_centre_rt == row.stage1_map_rt, "one-pass stage 3 is centred on the stage-1 map");
    }
    std::size_t early = 0;
    for (const auto& anchor : bent_result.anchors) early += anchor.used && anchor.library_rt < 10.0;
    require(early >= 15, "early anchors on a bent gradient " + std::to_string(early));
    for (double library_rt : {3.0, 8.0, 50.0})
    {
      require(std::fabs(bent_result.runRT(library_rt) - bentRT(library_rt)) < 6.0,
              "bent RT map at " + std::to_string(library_rt) + ": " + std::to_string(bent_result.runRT(library_rt)) + " vs " +
                std::to_string(bentRT(library_rt)));
    }
  }

  // Out-of-order spectra are rejected rather than silently mis-ordered.
  OpenDIA::CalibrationCollector unordered(library, traml_ids, options);
  bool threw = false;
  for (double rt : {20.0, 10.0})
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setRT(rt);
    Precursor precursor;
    precursor.setMZ(500.0);
    precursor.setIsolationWindowLowerOffset(60.0);
    precursor.setIsolationWindowUpperOffset(60.0);
    spectrum.setPrecursors({precursor});
    try
    {
      unordered.consumeSpectrum(spectrum);
    }
    catch (const std::runtime_error&)
    {
      threw = true;
    }
  }
  require(threw, "out-of-order rejection");
  std::cout << "OpenDIACalibration tests passed: anchors=" << result.anchor_targets << " rt_sigma=" << result.rt_map.sigma
            << " ppm=" << result.ms2_offset_ppm << " im=" << result.runIM(1.0) - 1.0 << '\n';
  return 0;
}
