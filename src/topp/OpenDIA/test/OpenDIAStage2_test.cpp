// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAStage2.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
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

  // 4. Cross-fitted LDA: a weak feature 0 plus an informative feature 2 separates better than feature 0 alone.
  {
    std::mt19937 rng(7);
    std::normal_distribution<double> noise(0.0, 1.0);
    std::vector<OpenDIA::Stage2Row> rows;
    for (int k = 0; k < 3000; ++k)
    {
      const bool decoy = k % 2 == 1;
      const bool real = !decoy && k % 4 == 0;  // half the targets are true
      OpenDIA::Stage2Row row;
      row.features.assign(OpenDIA::STAGE2_FEATURES, 0.0);
      row.precursor_id = "p" + std::to_string(k);
      row.is_decoy = decoy;
      row.fold = k % 3;
      row.features[0] = noise(rng) + (real ? 1.5 : 0.0);
      row.features[2] = noise(rng) + (real ? 3.0 : 0.0);
      for (std::size_t p : {1u, 3u, 4u, 5u, 6u, 7u}) row.features[p] = noise(rng);
      rows.push_back(row);
    }
    std::vector<OpenDIA::Stage2FoldReport> report;
    const auto score = OpenDIA::crossfitStage2LDA(rows, 3, OpenDIA::Stage2Options{}, report);
    require(report.size() == 3 && score.size() == rows.size(), "one report per fold, one score per row");
    for (const auto& fold : report)
    {
      require(fold.weights[2] > fold.weights[0] && fold.weights[2] > 0.0, "informative feature gets the largest weight");
      require(fold.weights[2] > 5.0 * std::abs(fold.weights[5]), "noise features get small weights");
    }
    // AUC of real targets against decoys, per fold (scores are fold-specific).
    for (int fold = 0; fold < 3; ++fold)
    {
      std::size_t wins = 0, pairs = 0, wins0 = 0;
      for (std::size_t i = 0; i < rows.size(); ++i)
      {
        if (rows[i].fold != fold || rows[i].is_decoy || std::stoi(rows[i].precursor_id.substr(1)) % 4 != 0) continue;
        for (std::size_t j = 0; j < rows.size(); ++j)
        {
          if (rows[j].fold != fold || !rows[j].is_decoy) continue;
          ++pairs;
          wins += score[i] > score[j];
          wins0 += rows[i].features[0] > rows[j].features[0];
        }
      }
      require(static_cast<double>(wins) / pairs > 0.95 && wins > wins0, "LDA beats feature 0 alone on held-out rows");
    }
  }
  // 5. Fail closed: a fold without enough positives, and non-finite features.
  {
    std::vector<OpenDIA::Stage2Row> rows(30);
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
      rows[i].features.assign(OpenDIA::STAGE2_FEATURES, 0.0);
      rows[i].precursor_id = "q" + std::to_string(i);
      rows[i].is_decoy = i % 2;
      rows[i].fold = static_cast<int>(i % 3);
    }
    std::vector<OpenDIA::Stage2FoldReport> report;
    bool threw = false;
    try { OpenDIA::crossfitStage2LDA(rows, 3, OpenDIA::Stage2Options{}, report); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "too few positives is rejected");
    rows[0].features[3] = std::nan("");
    threw = false;
    try { OpenDIA::crossfitStage2LDA(rows, 3, OpenDIA::Stage2Options{}, report); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "non-finite features are rejected");
    rows[0].features.assign(OpenDIA::STAGE2_FEATURES, 0.0);
    rows[1].features.resize(3);
    threw = false;
    try { OpenDIA::crossfitStage2LDA(rows, 3, OpenDIA::Stage2Options{}, report); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "rows with differing feature counts are rejected");
  }

  // 7. Feature sets and search peak-group features.
  {
    require(OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("extraction")).size() == OpenDIA::STAGE2_FEATURES, "extraction set is the original 8");
    const auto search = OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search"));
    require(search.size() == 5 && search[0] == "cnn_z" && search[4] == "abs_library_rt_error", "search set = cnn_z + 4");
    require(OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("extraction_search")).size() == 12, "combined set = 8 + 4");
    const auto all = OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search+isotope"));
    require(all.size() == 8 && all[5] == "ms1_isotope_correlation" && all[7] == "ms1_isotope_coelution", "search+isotope = cnn_z + 4 + 3");
    bool threw = false;
    try { OpenDIA::parseStage2FeatureSet("other"); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "unknown feature set rejected");
    threw = false;
    try { OpenDIA::parseStage2FeatureSet("search+search"); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "repeated group rejected");

    const std::string path = "OpenDIAStage2_test_search.tsv";
    {
      std::ofstream out(path);
      out << "candidate_row\tprecursor_id\tis_decoy\trt_apex\tpeak_width_seconds\ttransition_support\tmatched_scans\thas_library_rt\tapex_rt_error_seconds\n"
          << "0\t11\t0\t100.5\t20\t4\t9\t1\t-12.5\n"
          << "1\t12\t1\t200\t0\t1\t0\t0\t33\n";
    }
    const auto peaks = OpenDIA::readSearchPeakFeatures(path);
    require(peaks.size() == 2 && peaks.at("0").precursor_id == "11" && peaks.at("0").rt_apex == 100.5, "rows keyed by candidate_row");
    const auto f = OpenDIA::stage2SearchFeatures(peaks.at("0"));
    require(f[0] == 4.0 && std::abs(f[1] - std::log1p(9.0)) < 1e-12 && std::abs(f[2] - std::log1p(20.0)) < 1e-12 && f[3] == 12.5,
            "support, log scans, log width, |RT error|");
    require(OpenDIA::stage2SearchFeatures(peaks.at("1"))[3] == 0.0, "no library RT gives zero RT error");
    {
      std::ofstream out(path);
      out << "candidate_row\tprecursor_id\trt_apex\tpeak_width_seconds\ttransition_support\tmatched_scans\thas_library_rt\tapex_rt_error_seconds\n"
          << "0\t11\t100\t20\t4\t9\t1\tnan\n";
    }
    threw = false;
    try { OpenDIA::readSearchPeakFeatures(path); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "non-finite search value rejected");
    {
      std::ofstream out(path);
      out << "candidate_row\tprecursor_id\trt_apex\n0\t11\t100\n";
    }
    threw = false;
    try { OpenDIA::readSearchPeakFeatures(path); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "missing search column rejected");
    std::remove(path.c_str());
  }

  // 8. Isotope-envelope features: averagine-shaped co-eluting isotopes score high, a dominant M-1 peak is flagged.
  {
    OpenDIA::Stage2ExtractionOptions options;
    const double mz = 600.0;
    const int charge = 2;
    const auto& averagine = DIAHelpers::averagineIsotopeDistribution(mz * charge, 4);
    OpenDIA::TraceEvidence e;
    const auto shape = gaussian(1.0, 15.5, 2.0);
    for (std::size_t b = 12; b <= 19; ++b)
    {
      e.ms1_mask |= 1u << b;
      e.ms1[b] = static_cast<float>(1000.0 * shape[b] * averagine[0].getIntensity());
      e.ms1_isotope_traces[0][b] = static_cast<float>(1000.0 * shape[b] * averagine[1].getIntensity());
      e.ms1_isotope_traces[1][b] = static_cast<float>(1000.0 * shape[b] * averagine[2].getIntensity());
    }
    for (std::size_t k = 0; k < 4; ++k) e.ms1_isotope_core[k + 1] = static_cast<float>(1000.0 * averagine[k].getIntensity());
    auto f = OpenDIA::stage2IsotopeFeatures(e, mz, charge, options);
    require(f[0] > 0.99, "averagine-shaped envelope correlates");
    require(f[1] == 0.0, "no M-1 signal, no overlap");
    require(f[2] > 0.99, "co-eluting isotope traces correlate");
    e.ms1_isotope_core[0] = 3000.0f * static_cast<float>(averagine[0].getIntensity());
    std::reverse(e.ms1_isotope_core.begin() + 1, e.ms1_isotope_core.end());
    f = OpenDIA::stage2IsotopeFeatures(e, mz, charge, options);
    require(f[0] < 0.0, "reversed envelope anti-correlates");
    require(f[1] > 1.0, "larger M-1 peak gives overlap above 1");
    require(OpenDIA::stage2IsotopeFeatures(e, mz, 0, options)[0] == 0.0, "unknown charge gives zero features");
  }

  // 6. TraceCollector evidence: MS1 trace, mass errors and isotope; tensor unchanged by enabling evidence.
  {
    OpenSwath::LightTargetedExperiment lib;
    OpenSwath::LightCompound compound;
    compound.id = "1";
    compound.charge = 2;
    lib.compounds.push_back(compound);
    for (int k = 0; k < 3; ++k)
    {
      OpenSwath::LightTransition t;
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
    require(std::abs(e.ms1_isotope_core[2] - e.ms1_isotope) < 1e-3f, "core M+1 equals the first-isotope sum");
    require(e.ms1_isotope_core[1] > 0.0f && e.ms1_isotope_core[0] == 0.0f && e.ms1_isotope_core[3] == 0.0f, "M+0 present, M-1 and M+2 absent");
    float m1_trace = 0.0f;
    for (float v : e.ms1_isotope_traces[0]) m1_trace += v;
    require(m1_trace >= e.ms1_isotope, "M+1 trace covers at least the core");
  }

  // 9. Boundary features: picked on scan bins only; a narrow co-eluting peak in library ratios scores high, and a
  //    neighbouring peak outside the picked boundaries does not dilute co-elution or library correlation.
  {
    require(OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search+boundary")).size() == 10, "search+boundary = cnn_z + 4 + 5");
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

  // 10. b/y ion group: appended last, counts copied from the evidence.
  {
    const auto names = OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search+boundary+byion"));
    require(names.size() == 12 && names[10] == "byion_b" && names[11] == "byion_y", "byion columns come last");
    auto e = fullEvidence();
    e.byion_b = 7;
    e.byion_y = 9;
    const auto f = OpenDIA::stage2ByIonFeatures(e);
    require(f[0] == 7.0 && f[1] == 9.0, "byion features copy the counts");
    bool repeated = false;
    try { OpenDIA::parseStage2FeatureSet("byion+byion"); } catch (const std::invalid_argument&) { repeated = true; }
    require(repeated, "repeated byion group accepted");
  }

  // 11. Shape group: a co-eluting library-ratio peak scores near-perfect; a shifted fragment raises the lag; empty is zero.
  {
    const auto names = OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search+byion+shape"));
    require(names.size() == 20 && names[7] == "shape_xcorr_coelution" && names.back() == "shape_min_fragment_share", "shape columns come last");
    std::vector<std::vector<double>> p;
    for (double l : library) p.push_back(gaussian(l, 16.0, 2.0));  // area proportional to library intensity
    for (auto& v : p) for (auto& x : v) x += 0.01;                         // small flat background, so the window median is the noise
    const auto t = tensor(p);
    const auto f = OpenDIA::stage2ShapeFeatures(t.data(), library, mz, fullEvidence());
    require(f[0] < 0.01 && f[1] > 0.99 && f[4] > 0.99, "co-eluting peak: zero lag and shape ~1");
    require(f[5] > 0.99 && f[8] > 0.9 && f[6] < 0.05 && f[7] < 0.02, "library-ratio peak: dot product ~1, small distances");
    require(f[9] > 3.0 && f[11] == 6.0 && f[10] > 0.5 && f[12] > 0.0, "signal features of a strong peak");
    auto shifted = p;
    shifted[5] = gaussian(library[5], 18.0, 2.0);
    for (auto& x : shifted[5]) x += 0.01;
    const auto t2 = tensor(shifted);
    const auto g = OpenDIA::stage2ShapeFeatures(t2.data(), library, mz, fullEvidence());
    require(g[0] > f[0] && g[4] < f[4], "a shifted fragment must raise the lag and lower the worst shape");
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    const auto z = OpenDIA::stage2ShapeFeatures(empty.data(), library, mz, fullEvidence());
    require(z[1] == 0.0 && z[5] == 0.0 && z[11] == 0.0, "empty window gives zero shape features");
  }

  // 12. Refine group: clean co-elution is all quantitative; one interfered fragment is flagged and the median keeps the rest.
  {
    const auto names = OpenDIA::stage2FeatureSetNames(OpenDIA::parseStage2FeatureSet("search+refine"));
    require(names.size() == 11 && names[5] == "refine_quantitative" && names.back() == "refine_best_reference_r", "refine columns come last");
    std::vector<std::vector<double>> p;
    for (double l : library) p.push_back(gaussian(l, 16.0, 2.0));
    const auto t = tensor(p);
    const auto f = OpenDIA::stage2RefineFeatures(t.data(), library, mz, fullEvidence());
    require(f[0] == 6.0 && f[1] == 6.0 && f[2] > 0.99 && std::abs(f[3] - 1.0) < 1e-12 && std::abs(f[4] - 1.0) < 1e-12 && f[5] > 0.99,
            "co-eluting peak: all six fragments quantitative");
    // Fragment 1 (second most intense) carries a co-eluting peak plus a stronger interfering peak four bins later.
    auto interfered = p;
    const auto other = gaussian(2.0 * library[1], 20.0, 1.5);
    for (std::size_t b = 0; b < RT; ++b) interfered[1][b] += other[b];
    const auto t2 = tensor(interfered);
    const auto r = OpenDIA::refineFragmentCorrelations(t2.data(), mz, fullEvidence());
    require(r[1] < 0.75, "interfered fragment falls below the usable cut");
    require(r[0] >= 0.9 && r[2] >= 0.9 && r[3] >= 0.9 && r[4] >= 0.9 && r[5] >= 0.9, "clean fragments stay quantitative");
    const auto g = OpenDIA::stage2RefineFeatures(t2.data(), library, mz, fullEvidence());
    require(g[0] == 5.0 && g[1] == 5.0 && g[3] < 1.0 && g[4] < 0.75 && g[5] < f[5], "one interfered fragment: 5 quantitative, lower shares");
    // Six unrelated peaks (no shared elution): few or no fragments match the median shape.
    std::vector<std::vector<double>> scattered;
    for (std::size_t c = 0; c < 6; ++c) scattered.push_back(gaussian(library[c], 8.0 + 3.0 * c, 1.2));
    const auto t3 = tensor(scattered);
    const auto s = OpenDIA::stage2RefineFeatures(t3.data(), library, mz, fullEvidence());
    require(s[0] < f[0] && s[2] < 0.75 && s[5] < 0.75, "scattered fragments: low refinement scores");
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    const auto z = OpenDIA::stage2RefineFeatures(empty.data(), library, mz, fullEvidence());
    require(z[0] == 0.0 && z[2] == 0.0 && z[5] == 0.0, "empty window gives zero refine features");
    auto e = fullEvidence();
    e.ms2_mask = 0b11u;
    require(OpenDIA::stage2RefineFeatures(t.data(), library, mz, e)[0] == 0.0, "fewer than 5 scan bins gives zero refine features");
    // Mask: only the interfered channel is zeroed; a clean window and an unrefinable window are untouched.
    auto masked = t2;
    require(OpenDIA::maskInterferedFragments(masked.data(), mz, fullEvidence(), 0.75) == 1, "one channel masked");
    const std::size_t block = RT * IM;
    for (std::size_t c = 0; c < 6; ++c)
    {
      const bool zero = std::all_of(masked.begin() + c * block, masked.begin() + (c + 1) * block, [](float v) { return v == 0.0f; });
      const bool same = std::equal(masked.begin() + c * block, masked.begin() + (c + 1) * block, t2.begin() + c * block);
      require(c == 1 ? zero : same, "mask zeroes only the interfered channel");
    }
    auto clean = t;
    require(OpenDIA::maskInterferedFragments(clean.data(), mz, fullEvidence(), 0.75) == 0 && clean == t, "clean window unchanged");
    auto unrefinable = t2;
    require(OpenDIA::maskInterferedFragments(unrefinable.data(), mz, e, 0.75) == 0 && unrefinable == t2, "unrefinable window unchanged");
  }

  std::cout << "OpenDIAStage2_test passed\n";
  return 0;
}
