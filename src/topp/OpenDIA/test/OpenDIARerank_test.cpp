// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIARerank.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS::OpenDIA;

namespace
{
  void require_(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  std::vector<CandidateFeatures> makeRun_(const std::string& prefix, std::mt19937& rng)
  {
    constexpr OpenMS::Size assays_per_class = 120;
    std::normal_distribution<double> noise(0.0, 0.45);
    std::vector<CandidateFeatures> rows;
    rows.reserve(4 * assays_per_class);
    OpenMS::Size row = 0;
    for (int decoy = 0; decoy < 2; ++decoy)
    {
      for (OpenMS::Size i = 0; i < assays_per_class; ++i)
      {
        const double base = std::exp((decoy ? 1.2 : 2.0) + noise(rng));
        for (int alternative = 0; alternative < 2; ++alternative)
        {
          CandidateFeatures candidate;
          candidate.row = row++;
          candidate.precursor_id = prefix + (decoy ? "decoy_" : "target_") + std::to_string(i);
          candidate.precursor_mz = 400.0 + static_cast<double>(i) * 0.01 + (decoy ? 100.0 : 0.0);
          candidate.is_decoy = decoy != 0;
          candidate.area = base * (alternative ? 0.4 : 1.0);
          candidate.apex_intensity = candidate.area * (0.7 + 0.1 * static_cast<double>(i % 3));
          candidate.peak_width_seconds = 3.0 + static_cast<double>(i % 8);
          candidate.transition_support = 3 + i % 5;
          candidate.matched_scans = 5 + i % 7;
          candidate.has_library_rt = i % 4 != 0;
          candidate.apex_rt_error_seconds = candidate.has_library_rt ? noise(rng) * 10.0 : 0.0;
          candidate.has_library_im = i % 3 != 0;
          candidate.has_im_match_evidence = i % 5 != 0;
          candidate.mean_abs_im_error_vssc = candidate.has_im_match_evidence
                                                ? (decoy ? 0.02 : 0.006) + std::abs(noise(rng)) * 0.001
                                                : 0.0;
          rows.push_back(std::move(candidate));
        }
      }
    }
    return rows;
  }
}

int main()
{
  try
  {
    constexpr char feature_table[] =
      "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\trt_start\trt_apex\trt_end\t"
      "peak_width_seconds\tarea\tapex_intensity\ttransition_support\tmatched_scans\t"
      "has_library_rt\tapex_rt_error_seconds\thas_library_im\n"
      "0\tassayA\t500\t0\t1\t2\t3\t2\t10\t7\t4\t6\t1\t-0.5\t1\n";
    std::istringstream parsed_input(feature_table);
    const auto parsed = readCandidateFeatures(parsed_input);
    require_(parsed.size() == 1 && parsed[0].precursor_id == "assayA" && parsed[0].area == 10.0,
             "candidate feature TSV parser changed column mapping");
    require_(!parsed[0].has_im_match_evidence && parsed[0].mean_abs_im_error_vssc == 0.0,
             "legacy candidate feature TSV did not default missing mobility evidence");

    const std::string mobility_table =
      "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\trt_start\trt_apex\trt_end\t"
      "peak_width_seconds\tarea\tapex_intensity\ttransition_support\tmatched_scans\t"
      "has_library_rt\tapex_rt_error_seconds\thas_library_im\thas_im_match_evidence\tmean_abs_im_error_vssc\n"
      "0\tassayA\t500\t0\t1\t2\t3\t2\t10\t7\t4\t6\t1\t-0.5\t1\t1\t0.0125\n";
    std::istringstream mobility_input(mobility_table);
    const auto mobility_parsed = readCandidateFeatures(mobility_input);
    require_(mobility_parsed.size() == 1 && mobility_parsed[0].has_im_match_evidence &&
             std::abs(mobility_parsed[0].mean_abs_im_error_vssc - 0.0125) < 1e-12,
             "mobility evidence columns were not parsed");

    const std::string negative_count_table =
      "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\trt_start\trt_apex\trt_end\t"
      "peak_width_seconds\tarea\tapex_intensity\ttransition_support\tmatched_scans\t"
      "has_library_rt\tapex_rt_error_seconds\thas_library_im\n"
      "0\tassayA\t500\t0\t1\t2\t3\t2\t10\t7\t-1\t6\t1\t-0.5\t1\n";
    bool rejected_negative_count = false;
    try
    {
      std::istringstream invalid_input(negative_count_table);
      (void)readCandidateFeatures(invalid_input);
    }
    catch (const std::exception&)
    {
      rejected_negative_count = true;
    }
    require_(rejected_negative_count, "negative transition support was accepted");

    CandidateFeatures invalid_feature = parsed.front();
    invalid_feature.area = std::numeric_limits<double>::quiet_NaN();
    bool rejected_non_finite_feature = false;
    try
    {
      (void)rescoreCandidates({invalid_feature}, parsed, 17, 1);
    }
    catch (const std::exception&)
    {
      rejected_non_finite_feature = true;
    }
    require_(rejected_non_finite_feature, "non-finite model feature was accepted");

    std::mt19937 rng(17);
    const auto training = makeRun_("train_", rng);
    const auto scoring = makeRun_("heldout_", rng);
    const auto results = rescoreCandidates(training, scoring, 17, 1);
    require_(results.size() == scoring.size(), "scorer did not return one row per held-out candidate");

    OpenMS::Size winners = 0;
    OpenMS::Size target_winners = 0;
    OpenMS::Size decoy_winners = 0;
    for (const auto& result : results)
    {
      require_(std::isfinite(result.score), "Percolator returned a non-finite candidate score");
      if (result.assay_best)
      {
        ++winners;
        result.candidate.is_decoy ? ++decoy_winners : ++target_winners;
        require_(result.assay_q_value >= 0.0 && result.assay_q_value <= 1.0,
                 "assay-level q-value is outside [0, 1]");
      }
    }
    require_(winners == 240 && target_winners == 120 && decoy_winners == 120,
             "multiple candidate windows were not reduced to one target/decoy winner per assay");

    std::ostringstream output;
    writeScoredCandidates(output, results);
    require_(output.str().find("assay_best\tassay_q_value\n") != std::string::npos,
             "scored output is missing assay-level fields");
    std::cout << "candidate rescoring pilot checks passed\n";
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
