// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAStage2.h"

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS;
using OpenDIA::TraceCollector;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  constexpr std::size_t RT = TraceCollector::RT_BINS;
  constexpr std::size_t IM = TraceCollector::IM_BINS;

  /// Tensor with the given per-fragment RT profiles in the central IM bin.
  std::vector<float> tensor(const std::vector<std::vector<double>>& profiles)
  {
    std::vector<float> t(TraceCollector::SAMPLE_SIZE, 0.0f);
    for (std::size_t c = 0; c < profiles.size(); ++c)
    {
      for (std::size_t b = 0; b < RT; ++b) t[(c * RT + b) * IM + IM / 2] = static_cast<float>(profiles[c][b]);
    }
    return t;
  }

  std::vector<double> gaussian(double scale, double centre_bin, double width_bins)
  {
    std::vector<double> v(RT);
    for (std::size_t b = 0; b < RT; ++b) v[b] = scale * std::exp(-0.5 * std::pow((b - centre_bin) / width_bins, 2));
    return v;
  }

  OpenDIA::TraceEvidence fullEvidence()
  {
    OpenDIA::TraceEvidence e;
    e.ms1_mask = e.ms2_mask = 0xFFFFFFFFu;
    return e;
  }
}

int main()
{
  OpenDIA::Stage2ExtractionOptions options;  // 15 s half width, 4 s core: bins 11..20 are core
  const std::array<double, 6> library{100, 80, 60, 40, 20, 10};
  const std::array<double, 6> mz{300, 400, 500, 600, 700, 800};

  // 1. A real peptide: co-eluting fragments in library ratios, MS1 at the right mass, isotope present.
  {
    std::vector<std::vector<double>> p;
    for (double l : library) p.push_back(gaussian(l, 15.5, 2.0));
    const auto t = tensor(p);
    auto e = fullEvidence();
    const auto ms1 = gaussian(1000, 15.5, 2.0);
    for (std::size_t b = 0; b < RT; ++b) e.ms1[b] = static_cast<float>(ms1[b]);
    e.ms1_weight = 10.0; e.ms1_ppm_sum = 10.0 * 1.0;   // +1 ppm
    e.ms2_weight = 5.0; e.ms2_ppm_sum = 5.0 * -0.5;    // -0.5 ppm
    e.ms1_isotope = 3.0f;
    options.ms1_offset_ppm = 1.0;
    options.ms2_offset_ppm = -0.5;
    const auto f = OpenDIA::stage2Features(t.data(), library, mz, e, options);
    require(f[0] == 0.0, "feature 0 is left for the caller");
    require(f[1] == 1.0, "all fragments carry signal");
    require(f[2] > 0.99, "co-eluting fragments correlate");
    require(f[3] > 0.99, "library ratios match");
    require(std::abs(f[4]) < 1e-9 && std::abs(f[5]) < 1e-9, "mass errors relative to the run offset are zero");
    require(f[6] > 0.99, "MS1 co-elutes with the fragments");
    require(f[7] == 1.0, "isotope present");
  }
  // 2. Unrelated signal: each fragment peaks somewhere else, reversed library ratios, no MS1.
  {
    std::vector<std::vector<double>> p;
    // Core bins are 12..19; keep every displaced peak inside them so only ratios and timing differ.
    for (std::size_t c = 0; c < 6; ++c) p.push_back(gaussian(library[5 - c], 13.0 + 1.0 * c, 0.6));
    const auto t = tensor(p);
    auto e = fullEvidence();
    const auto f = OpenDIA::stage2Features(t.data(), library, mz, e, options);
    require(f[2] < 0.3, "non-co-eluting fragments do not correlate");
    require(f[3] < -0.9, "reversed library ratios anti-correlate");
    require(f[4] == options.ppm_cap && f[5] == options.ppm_cap, "missing signal gets the ppm cap");
    require(f[7] == 0.0, "no isotope without MS1");
  }
  // 3. Bins without an MS2 scan are excluded rather than read as zero signal.
  {
    std::vector<std::vector<double>> p(2, gaussian(10, 15.5, 2.0));
    const auto t = tensor(p);
    OpenDIA::TraceEvidence e;
    e.ms2_mask = (1u << 14) | (1u << 15);  // only two core bins scanned
    const auto f = OpenDIA::stage2Features(t.data(), library, {300, 400, 0, 0, 0, 0}, e, options);
    require(f[2] == 0.0, "fewer than three scanned bins give no co-elution");
    require(f[1] == 1.0, "two valid fragments, both with signal");
  }

  // 7. Search peak-group features.
  {
    OpenDIA::SearchPeakFeatures peak;
    peak.transition_support = 4;
    peak.matched_scans = 9;
    peak.peak_width_seconds = 20;
    peak.abs_rt_error_seconds = 12.5;
    const auto f = OpenDIA::stage2SearchFeatures(peak);
    require(f[0] == 4.0 && std::abs(f[1] - std::log1p(9.0)) < 1e-12 && std::abs(f[2] - std::log1p(20.0)) < 1e-12 && f[3] == 12.5,
            "support, log scans, log width, |RT error|");
  }

  // 6. TraceCollector evidence: MS1 trace, mass errors and isotope; tensor unchanged by enabling evidence.
  {
    OpenDIA::Library lib;
    OpenDIA::LibraryCompound compound;
    compound.id = "1";
    compound.charge = 2;
    lib.compounds.push_back(compound);
    for (int k = 0; k < 3; ++k)
    {
      OpenDIA::LibraryTransition t;
      t.transition_name = "t" + std::to_string(k);
      t.peptide_ref = "1";
      t.precursor_mz = 500.0;
      t.product_mz = 300.0 + 100.0 * k;
      t.library_intensity = 10.0 - k;
      lib.transitions.push_back(t);
    }
    const auto run = [&](bool evidence)
    {
      TraceCollector collector({{"s1", "1", "g1", -1, 100.0}}, lib);
      if (evidence) collector.enableEvidence(20.0, 4.0);
      for (double rt = 90.0; rt <= 110.0; rt += 1.0)
      {
        MSSpectrum ms1;
        ms1.setMSLevel(1);
        ms1.setRT(rt);
        ms1.push_back(Peak1D(500.0 * (1 + 2e-6), 100.0));               // +2 ppm
        ms1.push_back(Peak1D(500.0 + 1.0033548 / 2.0, 40.0));           // first isotope
        collector.consumeSpectrum(ms1);
        MSSpectrum ms2;
        ms2.setMSLevel(2);
        ms2.setRT(rt + 0.5);
        Precursor precursor;
        precursor.setMZ(500.0);
        precursor.setIsolationWindowLowerOffset(5.0);
        precursor.setIsolationWindowUpperOffset(5.0);
        ms2.setPrecursors({precursor});
        ms2.push_back(Peak1D(300.0, 5.0));
        ms2.push_back(Peak1D(400.0 * (1 - 1e-6), 5.0));               // -1 ppm
        collector.consumeSpectrum(ms2);
      }
      return collector;
    };
    const auto plain = run(false);
    const auto with = run(true);
    require(plain.intensities() == with.intensities(), "enabling evidence leaves the tensor unchanged");
    const auto& e = with.evidence(0);
    require(e.ms1_mask != 0 && e.ms2_mask != 0, "both MS levels marked");
    require(std::abs(e.ms1_ppm_sum / e.ms1_weight - 2.0) < 1e-3, "MS1 error +2 ppm");
    require(std::abs(e.ms2_ppm_sum / e.ms2_weight + 0.5) < 1e-3, "fragment error is the intensity mean of 0 and -1 ppm");
    require(e.ms1_isotope > 0.0f, "first isotope recorded");
  }

  // 9. Boundary features: picked on scan bins only; a narrow co-eluting peak in library ratios scores high, and a
  //    neighbouring peak outside the picked boundaries does not dilute co-elution or library correlation.
  {
    std::vector<std::vector<double>> p;
    for (std::size_t c = 0; c < 6; ++c)
    {
      auto v = gaussian(library[c], 16.0, 1.5);
      // A second species at bin 26 with a reversed fragment pattern.
      const auto other = gaussian(library[5 - c], 26.0, 1.5);
      for (std::size_t b = 0; b < RT; ++b) v[b] += other[b];
      p.push_back(v);
    }
    const auto t = tensor(p);
    auto e = fullEvidence();
    const auto f = OpenDIA::stage2BoundaryFeatures(t.data(), library, mz, e, options);
    // FWHM of sigma 1.5 bins after one smoothing pass: 2.355 * sqrt(2.25 + 0.5) bins * 0.9375 s = 3.66 s.
    require(std::abs(std::expm1(f[0]) - 3.66) < 0.3, "boundary FWHM off");
    require(std::expm1(f[1]) < 12.0, "boundary width should stop before the neighbour");
    require(f[2] > 0.95 && f[3] > 0.95 && f[4] == 1.0, "co-eluting library-ratio peak should score high");
    // Every other bin without a scan: the picker sees only the scan bins, so the peak is still found and narrow.
    e.ms2_mask = 0x55555555u;
    const auto g = OpenDIA::stage2BoundaryFeatures(t.data(), library, mz, e, options);
    require(std::expm1(g[1]) < 12.0 && g[3] > 0.9, "sparse scan bins break the boundary features");
    // Fewer than 3 scan bins or no signal: all zero.
    e.ms2_mask = 0x3u;
    const auto z = OpenDIA::stage2BoundaryFeatures(t.data(), library, mz, e, options);
    require(z[0] == 0.0 && z[2] == 0.0 && z[4] == 0.0, "too few scan bins must give zero features");
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    require(OpenDIA::stage2BoundaryFeatures(empty.data(), library, mz, fullEvidence(), options)[4] == 0.0, "empty window gives zero features");
  }

  // 10. b/y ion features: counts copied from the evidence.
  {
    auto e = fullEvidence();
    e.byion_b = 7;
    e.byion_y = 9;
    const auto f = OpenDIA::stage2ByIonFeatures(e);
    require(f[0] == 7.0 && f[1] == 9.0, "byion features copy the counts");
  }

  // 12. Experimental evidence block (-extra_scores): copies of the extraction features, core areas and mobility.
  {
    std::vector<std::vector<double>> p;
    for (std::size_t c = 0; c < 6; ++c) p.push_back(gaussian(library[c], 16.0, 1.5));
    auto t = tensor(p);
    auto e = fullEvidence();
    const auto base = OpenDIA::stage2Features(t.data(), library, mz, e, options);
    const auto x = OpenDIA::stage2ExtraFeatures(t.data(), library, mz, e, options);
    require(OpenDIA::stage2ExtraFeatureNames().size() == OpenDIA::STAGE2_EXTRA_FEATURES, "extra names");
    for (std::size_t k = 0; k < 7; ++k) require(x[k] == base[k + 1], "extra block copies extraction features 1..7");
    require(x[8] > 0.0, "core MS2 area");
    // All signal in the central mobility bin (IM / 2): centroid offset 0.5 bin, central share 1, profiles identical.
    require(std::abs(x[9] - 0.5) < 1e-9 && std::abs(x[11] - 1.0) < 1e-9, "central mobility signal");
    // Move every fragment to mobility bin 0: offset 3.5 bins, central share 0.
    for (std::size_t c = 0; c < 6; ++c)
      for (std::size_t b = 0; b < RT; ++b) std::swap(t[(c * RT + b) * IM + IM / 2], t[(c * RT + b) * IM]);
    const auto y = OpenDIA::stage2ExtraFeatures(t.data(), library, mz, e, options);
    require(std::abs(y[9] - 3.5) < 1e-9 && y[11] == 0.0, "off-centre mobility signal");
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    const auto z = OpenDIA::stage2ExtraFeatures(empty.data(), library, mz, e, options);
    require(z[8] == 0.0 && z[9] == IM / 2.0 && z[10] == 0.0 && z[11] == 0.0, "empty window");
  }

  // 15. PRECHECK_MF matched-filter features: names, no signal, coherent peak, single interferer, robust z by hand.
  {
    const auto& names = OpenDIA::stage2MatchedFilterFeatureNames();
    require(names.size() == 9 && names[0] == "mt_z_med" && names[3] == "mt_n_z3" && names[8] == "mt_llr", "matched-filter names are stable");
    const auto e = fullEvidence();
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    for (double v : OpenDIA::stage2MatchedFilterFeatures(empty.data(), library, mz, e, options))
      require(v == 0.0, "no signal gives zero matched-filter features");

    // Coherent library-shaped peak at the centre bin on a flat background.
    std::vector<std::vector<double>> p;
    for (double l : library)
    {
      auto g = gaussian(50.0 * std::sqrt(l), 16.0, 1.5);
      for (double& v : g) v += 1.0;
      p.push_back(g);
    }
    const auto t = tensor(p);
    const auto f = OpenDIA::stage2MatchedFilterFeatures(t.data(), library, mz, e, options);
    require(f[3] == 6.0, "coherent peak: all six fragments z >= 3");
    require(f[4] > 0.9 && f[6] > 0.9, "coherent peak: high leave-one-out correlation and cosine");

    // Only one (weak library) fragment carries signal: dominance > 0 and fewer z >= 3 fragments.
    std::vector<std::vector<double>> q(6, std::vector<double>(RT, 1.0));
    q[5] = gaussian(500.0, 16.0, 1.5);
    for (double& v : q[5]) v += 1.0;
    const auto ti = tensor(q);
    const auto g = OpenDIA::stage2MatchedFilterFeatures(ti.data(), library, mz, e, options);
    require(g[7] > 0.0 && g[3] < f[3], "single interferer: dominance positive and fewer z >= 3");

    // By hand: box sums bg {1,3,3,3,4,4,3}, median 3, MAD 0, apex 20, so z = 17 / (sqrt(3) + 1).
    const std::vector<double> x{1, 0, 2, 1, 0, 6, 9, 5, 0, 1, 3, 0};
    require(std::abs(OpenDIA::stage2RobustZ(x, 6, 5, 8) - 17.0 / (std::sqrt(3.0) + 1.0)) < 1e-9, "robust z by hand");
  }

  {
    // Sibling features. Group 0: target peptide with 3 library precursors (precursor 5 has no window), p0 and p1 co-elute,
    // p1 also has a far window. Group 1: isolated p2. Group 2: decoy peptide p3, p4 co-eluting.
    const float NK = -1e30f;
    const std::vector<std::uint32_t> group{0, 0, 0, 1, 2, 2}, prec{0, 1, 1, 2, 3, 4}, counts{3, 1, 2};
    const std::vector<double> rt{100.0, 101.5, 200.0, 100.0, 50.0, 52.0};
    const std::vector<float> score{5.0f, NK, 7.0f, 3.0f, 2.0f, 4.0f};
    const auto s = OpenDIA::stage2SiblingFeatures(group, prec, rt, score, counts, NK);
    const std::vector<std::vector<float>> expect{{2, 1, 2, 1.5f}, {2, 1, 5, 1.5f}, {2, 0, 2, 3}, {0, 0, 2, 3}, {1, 1, 4, 2}, {1, 1, 2, 2}};
    require(s.size() == 6 * OpenDIA::STAGE2_SIBLING_FEATURES, "sibling feature size");
    for (std::size_t i = 0; i < 6; ++i)
      for (std::size_t k = 0; k < 4; ++k) require(std::abs(s[i * 4 + k] - expect[i][k]) < 1e-6f, "sibling feature value");
  }

  std::cout << "OpenDIAStage2_test passed\n";
  return 0;
}
