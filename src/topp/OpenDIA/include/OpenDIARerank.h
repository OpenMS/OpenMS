// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <iosfwd>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  struct CandidateFeatures
  {
    Size row = 0;
    std::string precursor_id;
    double precursor_mz = 0.0;
    bool is_decoy = false;
    double peak_width_seconds = 0.0;
    double area = 0.0;
    double apex_intensity = 0.0;
    Size transition_support = 0;
    Size matched_scans = 0;
    bool has_library_rt = false;
    double apex_rt_error_seconds = 0.0;
    bool has_library_im = false;
  };

  struct ScoredCandidate
  {
    CandidateFeatures candidate;
    double score = 0.0;
    bool assay_best = false;
    double assay_q_value = 1.0;
  };

  /// Read the candidate feature TSV written by OpenDIA -features_out.
  std::vector<CandidateFeatures> readCandidateFeatures(std::istream& input);

  /// Train on the area-best candidate per assay, score every held-out candidate, and compute assay-level q-values.
  std::vector<ScoredCandidate> rescoreCandidates(const std::vector<CandidateFeatures>& training,
                                                  const std::vector<CandidateFeatures>& scoring,
                                                  int seed = 1,
                                                  int threads = 1);

  /// Write candidate scores and assay-level winners/q-values in input row order.
  void writeScoredCandidates(std::ostream& output, const std::vector<ScoredCandidate>& candidates);
}
