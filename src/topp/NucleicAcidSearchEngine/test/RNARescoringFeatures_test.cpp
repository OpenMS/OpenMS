// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNARescoringFeatures.h"

#include <OpenMS/CONCEPT/ClassTest.h>

#include <algorithm>

using namespace OpenMS;

namespace
{
PeptideHit::PeakAnnotation annotation(const std::string& name, double mz, double intensity)
{
  PeptideHit::PeakAnnotation result;
  result.annotation = name;
  result.mz = mz;
  result.intensity = intensity;
  result.charge = -1;
  return result;
}

double feature(const RNARescoringFeatures& features, const std::vector<double>& row, const std::string& name)
{
  const auto& names = features.featureNames();
  const auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end()) { throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name); }
  return row[found - names.begin()];
}
} // namespace

START_TEST(RNARescoringFeatures, "$Id$")

RNARescoringFeatures::Options options;
options.ion_series = {"a-B", "c", "w", "y"};
options.min_charge = 1;
options.max_charge = 3;

START_SECTION(RNARescoringFeatures(Options options))
{
  RNARescoringFeatures::Options bad = options;
  bad.ion_series = {"q"};
  TEST_EXCEPTION(Exception::InvalidParameter, RNARescoringFeatures{bad});
  bad = options;
  bad.min_charge = 0;
  TEST_EXCEPTION(Exception::InvalidParameter, RNARescoringFeatures{bad});
  bad = options;
  bad.min_charge = 4;
  TEST_EXCEPTION(Exception::InvalidParameter, RNARescoringFeatures{bad});
}
END_SECTION

START_SECTION(const StringList& featureNames() const)
{
  RNARescoringFeatures basic(options);
  TEST_EQUAL(ListUtils::concatenate(basic.featureNames(), ","),
             "native_score,delta_next,delta_best,abs_ppm_error,isotope_offset,charge_1,charge_2,charge_3,length,"
             "modifications,adduct,frac_a-B,frac_c,frac_w,frac_y,matched_ions,matched_intensity,longest_series,coverage")
  RNARescoringFeatures::Options all = options;
  all.route = all.graph = all.diagnostics = all.flash = true;
  RNARescoringFeatures extended(all);
  TEST_EQUAL(extended.featureNames().size(), basic.featureNames().size() + 4 + 1 + 1 + 3)
  TEST_EQUAL(extended.featureNames().back(), "flash_snr")
}
END_SECTION

START_SECTION(parseAnnotation)
{
  auto parsed = RNARescoringFeatures::parseAnnotation("a2-B");
  TEST_EQUAL(parsed.has_value(), true)
  TEST_EQUAL(parsed->first, "a-B")
  TEST_EQUAL(parsed->second, 2)
  parsed = RNARescoringFeatures::parseAnnotation("w13");
  TEST_EQUAL(parsed->first, "w")
  TEST_EQUAL(parsed->second, 13)
  parsed = RNARescoringFeatures::parseAnnotation("y3-");
  TEST_EQUAL(parsed->first, "y")
  TEST_EQUAL(parsed->second, 3)
  TEST_EQUAL(RNARescoringFeatures::parseAnnotation("").has_value(), false)
  TEST_EQUAL(RNARescoringFeatures::parseAnnotation("w").has_value(), false)
  TEST_EQUAL(RNARescoringFeatures::parseAnnotation("c2-B").has_value(), false)
  TEST_EQUAL(RNARescoringFeatures::parseAnnotation("a3-H2O").has_value(), false)
  TEST_EQUAL(RNARescoringFeatures::parseAnnotation("M-H2O").has_value(), false)
}
END_SECTION

START_SECTION(computeSpectrum)
{
  RNARescoringFeatures::Options all = options;
  all.route = all.graph = all.diagnostics = all.flash = true;
  RNARescoringFeatures features(all);

  // 5 residues, 4 backbone rungs; y2 shares its peak with c2
  const std::vector<PeptideHit::PeakAnnotation> annotations = {
    annotation("c1", 100.0, 1.0), annotation("c2", 200.0, 2.0), annotation("w1", 300.0, 3.0),
    annotation("a2-B", 400.0, 4.0), annotation("y2", 200.0, 2.0), annotation("x1", 500.0, 5.0),
    annotation("c9", 600.0, 6.0)};

  RNARescoringFeatures::HitInput best;
  best.native_score = 10.0;
  best.sequence = NASequence::fromString("[m1A]UCGG");
  best.precursor_error_ppm = -2.5;
  best.isotope_offset = 1;
  best.charge = -2;
  best.annotations = &annotations;
  best.route_matched_intensity = 3.0;
  best.route_local_intensity = 6.0;
  best.graph_tag_length = 3;

  RNARescoringFeatures::HitInput same_sequence = best; // other charge hypothesis: not a competitor
  same_sequence.native_score = 8.0;
  same_sequence.charge = -3;
  same_sequence.annotations = nullptr;

  RNARescoringFeatures::HitInput runner_up;
  runner_up.native_score = 7.0;
  runner_up.sequence = NASequence::fromString("AUC[m5C]G");
  runner_up.charge = -5; // outside the one-hot range
  runner_up.adduct = true;

  RNARescoringFeatures::SpectrumInput spectrum;
  spectrum.total_intensity = 20.0;
  spectrum.diagnostic_allowed = {"m1A"};

  const auto rows = features.computeSpectrum({best, same_sequence, runner_up}, spectrum);
  TEST_EQUAL(rows.size(), 3)
  for (const auto& row : rows) { TEST_EQUAL(row.size(), features.featureNames().size()) }

  const auto& row = rows[0];
  TEST_REAL_SIMILAR(feature(features, row, "native_score"), 10.0)
  TEST_REAL_SIMILAR(feature(features, row, "delta_next"), 3.0)
  TEST_REAL_SIMILAR(feature(features, row, "delta_best"), 0.0)
  TEST_REAL_SIMILAR(feature(features, row, "abs_ppm_error"), 2.5)
  TEST_REAL_SIMILAR(feature(features, row, "isotope_offset"), 1.0)
  TEST_REAL_SIMILAR(feature(features, row, "charge_1"), 0.0)
  TEST_REAL_SIMILAR(feature(features, row, "charge_2"), 1.0)
  TEST_REAL_SIMILAR(feature(features, row, "length"), 5.0)
  TEST_REAL_SIMILAR(feature(features, row, "modifications"), 1.0)
  TEST_REAL_SIMILAR(feature(features, row, "adduct"), 0.0)
  TEST_REAL_SIMILAR(feature(features, row, "frac_c"), 0.5)
  TEST_REAL_SIMILAR(feature(features, row, "frac_w"), 0.25)
  TEST_REAL_SIMILAR(feature(features, row, "frac_a-B"), 0.25)
  TEST_REAL_SIMILAR(feature(features, row, "frac_y"), 0.25)
  // x1 (series not configured) and c9 (beyond the last rung) are ignored
  TEST_REAL_SIMILAR(feature(features, row, "matched_ions"), 5.0 / 16.0)
  TEST_REAL_SIMILAR(feature(features, row, "matched_intensity"), (1.0 + 2.0 + 3.0 + 4.0) / 20.0)
  TEST_REAL_SIMILAR(feature(features, row, "longest_series"), 0.5)
  // c1 -> 1, c2/a2-B -> 2, y2 -> 3, w1 -> 4
  TEST_REAL_SIMILAR(feature(features, row, "coverage"), 1.0)
  TEST_REAL_SIMILAR(feature(features, row, "route_intensity_ratio"), 0.5)
  TEST_REAL_SIMILAR(feature(features, row, "graph_tag_length"), 3.0)
  TEST_REAL_SIMILAR(feature(features, row, "diagnostic_supported_mods"), 1.0)
  TEST_REAL_SIMILAR(feature(features, row, "flash_corrected"), 0.0)

  TEST_REAL_SIMILAR(feature(features, rows[1], "delta_next"), 1.0)
  TEST_REAL_SIMILAR(feature(features, rows[1], "delta_best"), -2.0)
  TEST_REAL_SIMILAR(feature(features, rows[1], "charge_3"), 1.0)
  TEST_REAL_SIMILAR(feature(features, rows[1], "matched_ions"), 0.0)

  const auto& last = rows[2];
  TEST_REAL_SIMILAR(feature(features, last, "delta_next"), 7.0)
  TEST_REAL_SIMILAR(feature(features, last, "delta_best"), -3.0)
  TEST_REAL_SIMILAR(feature(features, last, "charge_1") + feature(features, last, "charge_2") + feature(features, last, "charge_3"), 0.0)
  TEST_REAL_SIMILAR(feature(features, last, "adduct"), 1.0)
  TEST_REAL_SIMILAR(feature(features, last, "diagnostic_supported_mods"), 0.0)
  TEST_REAL_SIMILAR(feature(features, last, "route_intensity_ratio"), 0.0)

  // Ties between different sequences leave no margin; FLASH context is per spectrum
  RNARescoringFeatures::HitInput tie = runner_up;
  tie.native_score = 10.0;
  spectrum.flash_cosine = 0.9;
  spectrum.flash_snr = 2.0;
  const auto tied = features.computeSpectrum({best, tie}, spectrum);
  TEST_REAL_SIMILAR(feature(features, tied[0], "delta_next"), 0.0)
  TEST_REAL_SIMILAR(feature(features, tied[1], "delta_next"), 0.0)
  TEST_REAL_SIMILAR(feature(features, tied[1], "flash_corrected"), 1.0)
  TEST_REAL_SIMILAR(feature(features, tied[1], "flash_cosine"), 0.9)
  TEST_REAL_SIMILAR(feature(features, tied[1], "flash_snr"), 2.0)

  TEST_EQUAL(features.computeSpectrum({}, spectrum).empty(), true)
}
END_SECTION

START_SECTION(computeRow)
{
  RNARescoringFeatures::Options all = options;
  all.route = all.graph = all.diagnostics = all.flash = true;
  RNARescoringFeatures features(all);
  const std::vector<PeptideHit::PeakAnnotation> annotations = {annotation("c1", 100.0, 1.0), annotation("w1", 300.0, 3.0)};

  RNARescoringFeatures::HitInput best;
  best.native_score = 10.0;
  best.sequence = NASequence::fromString("[m1A]UCGG");
  best.charge = -2;
  best.annotations = &annotations;
  RNARescoringFeatures::HitInput same_sequence = best;
  same_sequence.native_score = 8.0;
  same_sequence.charge = -3;
  RNARescoringFeatures::HitInput runner_up;
  runner_up.native_score = 7.0;
  runner_up.sequence = NASequence::fromString("AUC[m5C]G");
  runner_up.charge = -1;
  RNARescoringFeatures::HitInput tie = runner_up;
  tie.native_score = 10.0;

  RNARescoringFeatures::SpectrumInput spectrum;
  spectrum.total_intensity = 20.0;
  spectrum.flash_cosine = 0.5;
  // Every row equals the corresponding row of computeSpectrum, in any hit order
  for (const auto& hits : std::vector<std::vector<RNARescoringFeatures::HitInput>>{
         {best, same_sequence, runner_up}, {runner_up, same_sequence, best}, {best, tie}, {best}})
  {
    const auto rows = features.computeSpectrum(hits, spectrum);
    for (Size index = 0; index < hits.size(); ++index)
    {
      const auto row = features.computeRow(hits, index, spectrum);
      TEST_EQUAL(row.size(), rows[index].size())
      for (Size column = 0; column < row.size(); ++column) { TEST_REAL_SIMILAR(row[column], rows[index][column]) }
    }
  }
  TEST_REAL_SIMILAR(feature(features, features.computeRow({best, same_sequence, runner_up}, 0, spectrum), "delta_next"), 3.0)
  TEST_REAL_SIMILAR(feature(features, features.computeRow({best}, 0, spectrum), "delta_next"), 10.0)
}
END_SECTION

END_TEST
