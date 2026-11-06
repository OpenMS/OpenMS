// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cstddef>
#include <iosfwd>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  struct TransitionQuery
  {
    std::string_view precursor_id;
    double precursor_mz = 0.0;
    double fragment_mz = 0.0;
    double library_rt = 0.0;
    bool has_library_rt = false;
  };

  struct SearchOptions
  {
    double fragment_tolerance_da = 0.02;
    double rt_window_seconds = 0.0;
    double peak_gap_seconds = 30.0;
    std::size_t min_transitions = 1;
  };

  /**
    @brief Stream DIA spectra and report fragment-group peak candidates.

    Queries are sorted once by product m/z. Each centroided MS2 peak visits only
    nearby library queries. Candidate state is constant-size per precursor and
    completed groups are written immediately, so raw scans and chromatograms are
    never retained.
  */
  class StreamSearch final : public Interfaces::IMSDataConsumer
  {
  public:
    StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;

    /// Flush active groups at end of input.
    void finish();

  private:
    struct Candidate
    {
      std::string id;
      double precursor_mz = 0.0;
      double library_rt = 0.0;
      bool has_library_rt = false;
      bool active = false;
      double start_rt = 0.0;
      double apex_rt = 0.0;
      double last_rt = 0.0;
      double last_intensity = 0.0;
      double apex_intensity = 0.0;
      double area = 0.0;
      std::size_t points = 0;
      std::size_t max_transition_support = 0;
    };

    struct Query
    {
      double fragment_mz = 0.0;
      std::size_t candidate = 0;
    };

    void consumeSpectrum_(const MSSpectrum& spectrum);
    void updateCandidate_(std::size_t candidate, double rt, double intensity, std::size_t transition_support);
    void writeCandidate_(const Candidate& candidate);
    bool queryMatches_(const Query& query, double rt) const;

    SearchOptions options_;
    std::ostream& output_;
    std::vector<Candidate> candidates_;
    std::vector<Query> queries_;
    std::vector<double> query_intensity_;
    std::vector<std::size_t> query_generation_;
    std::vector<std::size_t> query_touched_;
    std::vector<double> candidate_scan_intensity_;
    std::vector<std::size_t> candidate_scan_support_;
    std::vector<std::size_t> candidate_generation_;
    std::vector<std::size_t> candidate_touched_;
    std::vector<std::pair<double, double>> isolation_windows_;
    std::size_t generation_ = 0;
    bool finished_ = false;
  };
}
