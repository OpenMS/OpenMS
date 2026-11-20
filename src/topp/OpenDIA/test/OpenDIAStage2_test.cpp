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
    bool removed = false;
    try { OpenDIA::parseStage2FeatureSet("search+isotope"); } catch (const std::invalid_argument&) { removed = true; }
    require(removed, "removed isotope group rejected");
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

  // 11. Peak-group scores (PLAN_PGSCORES).
  {
    const std::array<std::int8_t, 6> charge{1, 1, 1, 1, 1, 1};
    OpenDIA::PeakGroupEvidence pg;
    for (std::size_t c = 0; c < 6; ++c)
    {
      pg.mono[c] = static_cast<float>(library[c]);
      pg.ppm_sum[c] = static_cast<float>(library[c] * 2.0);  // +2 ppm each
      // Averagine-like decreasing isotopes and nothing before the monoisotope.
      for (std::size_t i = 0; i < 4; ++i) pg.isotope[c][i] = static_cast<float>(library[c] * std::pow(0.4, i + 1));
    }
    std::vector<std::vector<double>> p, shifted;
    for (std::size_t c = 0; c < 6; ++c)
    {
      p.push_back(gaussian(library[c], 16.0, 2.0));
      shifted.push_back(gaussian(library[c], 16.0 + (c % 2 == 0 ? -2.0 : 2.0), 2.0));
    }
    const auto good = OpenDIA::stage2PeakGroupFeatures(tensor(p).data(), library, mz, charge, fullEvidence(), pg, 10.0, 1.0);
    require(OpenDIA::stage2PeakGroupFeatureNames().size() == good.size(), "names match features");
    require(good[0] == 0.0 && good[1] == 0.0, "co-eluting traces have zero best lag");
    require(good[2] > 0.95 && good[3] > 0.95, "identical shapes correlate");
    require(good[4] > std::log1p(1.0), "a peak over a flat background has S/N above 1");
    require(good[5] > 0.9 && good[5] <= 1.0, "nearly all signal inside the peak");
    require(std::abs(good[6] - 0.2) < 1e-6 && std::abs(good[7] - 0.2) < 1e-6, "2 ppm of 10 ppm tolerance");
    require(good[8] > 0.8 && good[9] == 0.0, "decreasing isotopes correlate with a small-fragment averagine, no overlap");
    const auto bad = OpenDIA::stage2PeakGroupFeatures(tensor(shifted).data(), library, mz, charge, fullEvidence(), pg, 10.0, 1.0);
    require(bad[0] > good[0] && bad[2] < good[2], "shifted fragments have a best lag and a lower best shape");
    // A larger peak one 13C spacing below at charges 1 and 2 counts twice per fragment.
    auto overlap = pg;
    for (std::size_t c = 0; c < 6; ++c) overlap.before[c] = {1e6f, 1e6f, 0.0f, 0.0f};
    const auto o = OpenDIA::stage2PeakGroupFeatures(tensor(p).data(), library, mz, charge, fullEvidence(), overlap, 10.0, 1.0);
    require(std::abs(o[9] - 2.0) < 1e-6, "two larger preceding peaks per fragment");
    // No signal: lags 10, shapes 0, mass deviation 1.
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    const auto z = OpenDIA::stage2PeakGroupFeatures(empty.data(), library, mz, charge, fullEvidence(), OpenDIA::PeakGroupEvidence{}, 10.0, 1.0);
    require(z[0] == 10.0 && z[2] == 0.0 && std::abs(z[6] - 1.0) < 1e-12 && z[8] == 0.0, "empty window defaults");

    // Collection: core-window monoisotope, ppm and isotope sums per channel; the tensor is unchanged.
    OpenDIA::Library lib;
    OpenDIA::LibraryCompound compound;
    compound.id = "1";
    compound.charge = 2;
    lib.compounds.push_back(compound);
    for (int k = 0; k < 2; ++k)
    {
      OpenDIA::LibraryTransition t;
      t.transition_name = "t" + std::to_string(k);
      t.peptide_ref = "1";
      t.precursor_mz = 500.0;
      t.product_mz = 300.0 + 100.0 * k;
      t.library_intensity = 10.0 - k;
      t.fragment_charge = static_cast<std::int8_t>(k + 1);
      lib.transitions.push_back(t);
    }
    const auto collect = [&](bool pg_on)
    {
      TraceCollector collector({{"s1", "1", "g1", -1, 100.0}}, lib);
      collector.enableEvidence(20.0, 4.0);
      if (pg_on) collector.enablePeakGroupEvidence(10.0);
      for (double rt = 90.0; rt <= 110.0; rt += 1.0)
      {
        MSSpectrum ms2;
        ms2.setMSLevel(2);
        ms2.setRT(rt + 0.5);
        Precursor precursor;
        precursor.setMZ(500.0);
        precursor.setIsolationWindowLowerOffset(5.0);
        precursor.setIsolationWindowUpperOffset(5.0);
        ms2.setPrecursors({precursor});
        ms2.push_back(Peak1D(300.0 - 1.0033548, 50.0));            // before fragment 0 (charge 1)
        ms2.push_back(Peak1D(300.0 * (1 + 3e-6), 10.0));           // +3 ppm
        ms2.push_back(Peak1D(300.0 + 1.0033548, 4.0));             // M+1 of the singly charged fragment
        ms2.push_back(Peak1D(400.0, 10.0));
        ms2.push_back(Peak1D(400.0 + 1.0033548 / 2.0, 6.0));       // M+1 of the doubly charged fragment
        collector.consumeSpectrum(ms2);
      }
      return collector;
    };
    const auto off = collect(false);
    const auto on = collect(true);
    require(off.intensities() == on.intensities(), "peak-group evidence leaves the tensor unchanged");
    const auto& e = on.peakGroupEvidence(0);
    // Spectra at 96.5 .. 103.5 lie within the 4 s core: 8 spectra.
    require(std::abs(e.mono[0] - 80.0f) < 1e-3 && std::abs(e.ppm_sum[0] / e.mono[0] - 3.0f) < 1e-2, "core monoisotope and ppm");
    require(std::abs(e.isotope[0][0] - 32.0f) < 1e-3 && std::abs(e.isotope[1][0] - 48.0f) < 1e-3, "isotopes at the fragment charge");
    require(std::abs(e.before[0][0] - 400.0f) < 1e-3 && e.before[1][0] == 0.0f, "preceding peak at charge 1");
    require(on.fragmentCharges(0)[0] == 1 && on.fragmentCharges(0)[1] == 2, "fragment charges kept");
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

  // 13. Own-mobility-apex block (-mobility_apex): from the core-window 1/K0 histogram.
  {
    constexpr std::size_t B = TraceCollector::MOBILITY_HISTOGRAM_BINS;
    constexpr double W = TraceCollector::MOBILITY_HISTOGRAM_BIN;
    require(OpenDIA::stage2MobilityApexFeatureNames().size() == OpenDIA::STAGE2_MOBILITY_APEX_FEATURES, "mobility apex names");
    const auto none = OpenDIA::stage2MobilityApexFeatures(nullptr, library, mz);
    require(std::abs(none[0] - B * W / 2.0) < 1e-12 && none[1] == 0.03 && none[2] == 0.0 && none[3] == 0.0 && none[4] == 0.0, "no histogram");
    // All fragments at bin 20 with the library pattern: apex offset = |centre(20)|, no spread, cosine 1, 6 supporting,
    // nothing in the side band.
    std::vector<float> h(TraceCollector::MOBILITY_HISTOGRAM_SIZE, 0.0f);
    for (std::size_t c = 0; c < 6; ++c) h[c * B + 20] = static_cast<float>(std::sqrt(library[c]) * 100.0);
    const auto a = OpenDIA::stage2MobilityApexFeatures(h.data(), library, mz);
    require(std::abs(a[0] - std::abs(20.5 * W - B * W / 2.0)) < 1e-12, "apex offset");
    require(a[1] < 1e-12 && std::abs(a[2] - 1.0) < 1e-9 && a[3] == 6.0, "coherent apex");
    double inside = 0.0;
    for (std::size_t c = 0; c < 6; ++c) inside += h[c * B + 20];
    require(std::abs(a[4] - std::log1p(inside)) < 1e-6, "contrast without side band");
    // One fragment moved 4 bins away: spread > 0, support 5, the side band still empty (4 < 5).
    std::swap(h[5 * B + 20], h[5 * B + 24]);
    const auto b = OpenDIA::stage2MobilityApexFeatures(h.data(), library, mz);
    require(b[1] > 0.0 && b[3] == 5.0 && b[2] < 1.0, "one fragment off the apex");
    // Padding channels are ignored.
    auto mz_pad = mz;
    mz_pad[5] = 0.0;
    const auto c = OpenDIA::stage2MobilityApexFeatures(h.data(), library, mz_pad);
    require(c[1] < 1e-12 && c[3] == 5.0 && std::abs(c[2] - 1.0) < 1e-9, "padding ignored");
  }

  // 14. PRECHECK_FAINT pooled tensor features: known apex cosine/local signal and no-signal behavior.
  {
    std::vector<std::vector<double>> p(6, std::vector<double>(RT, 0.0));
    for (std::size_t c = 0; c < 6; ++c)
      for (std::size_t b = 15; b <= 17; ++b) p[c][b] = std::sqrt(library[c]);
    const auto t = tensor(p);
    const auto e = fullEvidence();
    const auto f = OpenDIA::stage2FaintFeatures(t.data(), library, mz, e, options);
    require(std::abs(f[0] - 1.0) < 1e-12, "apex cosine is one for library-shaped signal");
    require(std::isfinite(f[1]) && f[1] >= 0.0, "local signal feature is finite and nonnegative for a central peak");
    const std::vector<float> empty(TraceCollector::SAMPLE_SIZE, 0.0f);
    const auto z = OpenDIA::stage2FaintFeatures(empty.data(), library, mz, e, options);
    require(z[0] == 0.0 && z[1] == 0.0, "no signal gives zero faint features");
    require(OpenDIA::stage2FaintFeatureNames()[0] == "p_cos_apex" && OpenDIA::stage2FaintFeatureNames()[1] == "p_local_sn",
            "faint feature names are stable");
  }

  std::cout << "OpenDIAStage2_test passed\n";
  return 0;
}
