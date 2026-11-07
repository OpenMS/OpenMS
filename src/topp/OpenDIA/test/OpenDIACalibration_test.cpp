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
