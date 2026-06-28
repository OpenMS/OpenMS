// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "FeatureTypes.h"
#include "GridWithStorage.h"
#include "MzDiff.h"
#include "Run.h"
#include "RunStatistics.h"
#include "Window.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

namespace OpenMS::PipEcho
{

/******************************************************************************/
// Handy aliases.
using DonorMap = GridWithStorage<Donor>;
using AcceptorMap = GridWithStorage<Acceptor>;
using RunMap = std::map<std::string, Run>;

/// [LOCAL-RT] Per-run-pair local RT calibration; defined in Impl.cpp.
struct LocalRtModel;

/******************************************************************************/
/**
 * Internal representation of the PIP-ECHO algorithm.
 */
class Impl
{
public:
  /// Type used when searching for a matching donor.
  using match_t = std::optional<std::pair<Score, std::shared_ptr<Acceptor>>>;

  /// Construct a new implementation object.
  Impl(const Param& params, const std::pair<double, double>& mz_range);

  /// Separate donors from acceptors.
  void partition_features(const std::vector<FeatureMap>&, RunMap&);

  /// Match identified features with unidentified features.
  void link_donors_and_acceptors(const RunStatistics&,
                                 const DonorMap&,
                                 AcceptorMap&);

  /// Search for a matching acceptor.
  match_t find_acceptor_for(const RunStatistics&,
                            const AcceptorMap&,
                            const Donor&,
                            const Window&,
                            const std::optional<double> = {});

  /// Find a random Donor that is dissimilar to a given Donor.
  std::optional<std::shared_ptr<const Donor>>
  find_random_donor(const DonorMap&, const Donor&, const Window&) const;

  /// [LOCAL-RT] Build the local RT calibration for the current (donor_run,
  /// acceptor_run) pair from peptides identified in BOTH runs. No-op unless the
  /// 'local_rt:enabled' parameter is set.
  void set_local_rt_model(const DonorMap& donor_donors,
                          const DonorMap& acceptor_donors);

  /// [LOCAL-RT auto] Estimate global window scales (cap/floor/locality/fallback)
  /// from the data before matching. No-op unless 'local_rt:auto' is enabled.
  void estimate_auto_params(const RunMap& runs);

  /// [LOCAL-RT] Whether the local adaptive RT window is enabled.
  bool local_rt_enabled() const { return use_local_rt_; }

  /// [LOCAL-RT] Decoy count needed to resolve the requested FDR (the adaptive-
  /// widening target): max(min_decoys, ceil(1/fdr)).
  std::size_t target_decoy_count() const
  {
    if (mbr_fdr >= 1.0) { return 0; }  // FDR control disabled (keep all) -> no decoy requirement
    const double f = mbr_fdr > 0.0 ? mbr_fdr : 1.0;
    return std::max<std::size_t>(min_decoys,
                                 static_cast<std::size_t>(std::ceil(1.0 / f)));
  }

  /// [LOCAL-RT] Begin a match pass at widening multiplier w: reseed the decoy RNG
  /// and reset the per-pass diagnostic counters, so each (re-)match pass is
  /// independent of discarded passes and the final log/result describe the
  /// retained pass (Codex review).
  void set_widen_factor(double w) { lrt_widen_ = w; rng_.seed(rng_seed_); lrt_supported_ = 0; lrt_fallback_ = 0; }

  /// [LOCAL-RT] Widening multiplier at which even the SMALLEST local window
  /// (the floor) saturates the global ceiling -- widening further cannot add
  /// decoys, so stop. (Based on floor, not cap: the data-driven c*sigma window
  /// is usually far below the cap, so the cap is the wrong saturation basis.)
  double widen_ceiling() const { return rt_sec_max_window / std::max(lrt_floor_, 1.0); }

  /// Fill in the final ConsensusMap.
  void generate_consensus_map(RunMap&, ConsensusMap&);

public:
  /// Max allowed m/z difference between donor and acceptor.
  MzDiff mz_max_diff;

  /// Grid center for the m/z dimension.  This is used to decide
  /// which features are close to one another.
  double mz_grid_center;

  /// Max allowed RT difference between donor and acceptor.
  double rt_sec_max_window;

  /// The FDR to use for the MBR process.
  double mbr_fdr;

  /// Minimum number of MBR decoy transfers required before transferred
  /// features are kept (the FDR-estimability floor; see TransferFDRModel::run).
  std::size_t min_decoys;

  /// Upper bound on the SVM training-set size per cross-validation fold
  /// (0 = unlimited). Only the model fit is subsampled; prediction and the
  /// FDR/q-values still use every transfer (see TransferFDRModel::round).
  std::size_t max_training_points;

private:
  std::string path_from_feature_map(const FeatureMap&);
  bool is_donor_feature(const Feature&);
  Run& get_run_from_file_name(RunMap&, const std::string&);

  std::optional<Window> initial_window(const Donor& donor);
  std::optional<Window> next_window(const std::optional<Window>&);

  /// Random source for decoy donor selection.  Seeded once (from the
  /// `random_seed` parameter) so results are reproducible for the current
  /// single-threaded decoy search.  NOT thread-safe: a shared std::mt19937
  /// must not be advanced concurrently -- give each thread its own generator
  /// before parallelizing the donor loop.  Mutable because the (logically
  /// const) decoy search advances its state.
  std::mt19937::result_type rng_seed_ = 0;  ///< decoy RNG seed (for per-pass reseed)
  mutable std::mt19937 rng_;

  /// [LOCAL-RT] local adaptive RT window state ('local_rt:' parameters).
  bool use_local_rt_ = false;
  bool use_auto_ = false;                          ///< auto-estimate window scales from data
  double host_fwhm_ = 0.0;                          ///< optional host chromatographic FWHM (s; 0 = none)
  std::shared_ptr<const LocalRtModel> local_rt_;   ///< current run-pair model
  double lrt_cap_ = 100.0;                         ///< generous half-window backstop (param)
  double lrt_floor_ = 5.0;                         ///< min half-window (param); widen-ceiling basis
  double lrt_locality_ = 120.0;                    ///< anchor RT neighborhood (param)
  std::size_t lrt_n_anchors_ = 3;                  ///< anchors per side (param)
  double lrt_c_ = 3.0;                             ///< sigma scale for the window (param)
  double lrt_fallback_half_ = 15.0;                ///< half-window when <2 anchors (param)
  double lrt_widen_ = 1.0;                         ///< adaptive-widening multiplier (set per pass)
  mutable std::size_t lrt_supported_ = 0;          ///< queries with a sigma-sized window
  mutable std::size_t lrt_fallback_ = 0;           ///< queries with a tight fallback window
};

} // namespace OpenMS::PipEcho
