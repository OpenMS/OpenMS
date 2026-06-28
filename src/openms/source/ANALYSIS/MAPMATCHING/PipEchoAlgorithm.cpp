// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include "PIPECHO/Impl.h"
#include "PIPECHO/Run.h"
#include "PIPECHO/RunStatistics.h"
#include "PIPECHO/Util.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>

namespace OpenMS
{

namespace
{
/******************************************************************************/
/// Find the smallest and largest m/z value among all maps.
std::pair<double, double> mz_range(const std::vector<FeatureMap>& feature_maps)
{
  double mz_min = std::numeric_limits<double>::max();
  double mz_max = std::numeric_limits<double>::lowest();

  for (auto& map : feature_maps)
  {
    // The input maps are const and their m/z ranges are not guaranteed to be
    // current here (updateRanges() cannot be called on a const map), so we
    // derive the bounds directly from the features. getMinMZ()/getMaxMZ() would
    // throw Exception::InvalidRange on an un-updated (empty) range.
    for (auto& feature : map)
    {
      mz_min = std::min(mz_min, feature.getMZ());
      mz_max = std::max(mz_max, feature.getMZ());
    }
  }

  return std::make_pair(mz_min, mz_max);
}
} // namespace

/******************************************************************************/
PipEchoAlgorithm::PipEchoAlgorithm(): FeatureGroupingAlgorithm()
{
  setName("PipEcho");

  defaults_.setValue("distance_RT:max_difference", 100.0,
                     "Never pair features with a larger RT distance"
                     " (in seconds).");
  defaults_.setMinFloat("distance_RT:max_difference", 0.0);


  defaults_.setValue("distance_MZ:max_difference", 10.0,
                     "Never pair features with larger m/z distance"
                     " (unit defined by 'unit')");
  defaults_.setMinFloat("distance_MZ:max_difference", 0.0);

  defaults_.setValue("distance_MZ:unit", "ppm", "Unit of the 'max_difference' parameter");
  defaults_.setValidStrings("distance_MZ:unit", {"Da", "ppm"});

  defaults_.setValue("fdr", 0.05, "MBR FDR threshold (0.05=5%).");
  defaults_.setMinFloat("fdr", 0.0);
  defaults_.setMaxFloat("fdr", 1.0);

  defaults_.setValue("random_seed", 0,
                     "Seed for the random number generator used to select decoy"
                     " donors. A fixed seed makes results reproducible.");

  defaults_.setValue("min_decoys", 20,
                     "Minimum number of MBR decoy transfers required before any"
                     " transferred feature is kept. Transferred features are also"
                     " dropped whenever the requested 'fdr' cannot be resolved by"
                     " the available decoys (1/decoys > fdr). When transfers are"
                     " dropped, only direct identifications are retained. A 'fdr'"
                     " of 1.0 disables FDR control and keeps all transfers"
                     " regardless of this value.",
                     {"advanced"});
  defaults_.setMinInt("min_decoys", 1);

  defaults_.setValue("max_training_points", 50000,
                     "Upper bound on the number of candidate transfers used to"
                     " TRAIN the transfer-FDR SVM in each cross-validation fold"
                     " (0 = unlimited). Predictions and the FDR/q-values are"
                     " always computed over ALL transfers, so this only bounds"
                     " the SVM-fitting cost and does not change which transfers"
                     " are scored. Runs with very many features per run can"
                     " otherwise make the SVM grid search slow. When the cap is"
                     " hit, a deterministic stratified subsample is used: both"
                     " labelled classes (decoys and positives) are kept whole if"
                     " they fit, otherwise the minority class is kept whole and"
                     " the majority is score-spread sampled to fit the budget.",
                     {"advanced"});
  defaults_.setMinInt("max_training_points", 0);

  // --- Local adaptive RT window (EXPERIMENTAL; default off = the global window) ---
  defaults_.setValue("local_rt:enabled", "false",
                     "Use a LOCAL adaptive retention-time window for MBR candidate"
                     " search instead of the single global RT window"
                     " ('distance_RT:max_difference'). For each donor the expected"
                     " acceptor RT is predicted from nearby peptides identified in BOTH"
                     " runs (local alignment) and the search window is sized from the"
                     " local RT scatter, sharpening the RT feature and removing much of"
                     " the false in-window background. The window is widened adaptively"
                     " if too few decoys are produced to resolve the FDR.",
                     {"advanced"});
  defaults_.setValidStrings("local_rt:enabled", {"true", "false"});

  defaults_.setValue("local_rt:max_window", 100.0,
                     "Maximum half-width (seconds) of the local RT window before"
                     " adaptive widening (backstop cap on the data-driven width).",
                     {"advanced"});
  defaults_.setMinFloat("local_rt:max_window", 0.0);

  defaults_.setValue("local_rt:min_window", 5.0,
                     "Minimum half-width (seconds) of the local RT window.", {"advanced"});
  defaults_.setMinFloat("local_rt:min_window", 0.0);

  defaults_.setValue("local_rt:anchor_window", 120.0,
                     "Only peptides within this RT distance (seconds) of the donor are"
                     " used as local-alignment anchors.", {"advanced"});
  defaults_.setMinFloat("local_rt:anchor_window", 0.0);

  defaults_.setValue("local_rt:anchors", 3,
                     "Maximum number of anchor peptides per side used for the local RT"
                     " prediction.", {"advanced"});
  defaults_.setMinInt("local_rt:anchors", 1);

  defaults_.setValue("local_rt:sigma_scale", 3.0,
                     "Local RT window half-width = this multiple of the local anchor"
                     " RT-shift standard deviation.", {"advanced"});
  defaults_.setMinFloat("local_rt:sigma_scale", 0.0);

  defaults_.setValue("local_rt:fallback_window", 15.0,
                     "Half-width (seconds) used when fewer than two local anchors are"
                     " available.", {"advanced"});
  defaults_.setMinFloat("local_rt:fallback_window", 0.0);

  defaults_.setValue("local_rt:auto", "false",
                     "Estimate the local RT window scales (max_window, min_window,"
                     " anchor_window, fallback_window) from the data instead of using the"
                     " fixed values above, so one configuration adapts across short and"
                     " long gradients. Derived from the RT-shift scatter of shared"
                     " anchors and the anchor density; uses 'median_fwhm' as a physical"
                     " lower guard when provided.", {"advanced"});
  defaults_.setValidStrings("local_rt:auto", {"true", "false"});

  defaults_.setValue("local_rt:median_fwhm", 0.0,
                     "Chromatographic peak FWHM (seconds) used as a physical lower guard"
                     " when 'auto' is enabled (0 = estimate window scales from RT"
                     " residuals and anchor density only). Set by the host tool"
                     " (e.g. ProteomicsLFQ) which has raw-spectra context.", {"advanced"});
  defaults_.setMinFloat("local_rt:median_fwhm", 0.0);

  defaults_.setValue("local_rt:rt_score", "svm_and_mbr",
                     "How the retention-time feature is used on the local RT path."
                     " 'raw': the raw |Δrt| is the SVM predictor (legacy)."
                     " 'svm': a calibrated [0,1] RT-agreement score (FlashLFQ-style"
                     " two-tailed CDF against the data-driven RT prediction-error"
                     " distribution) replaces it as the SVM predictor."
                     " 'svm_and_mbr': as 'svm', and the calibrated score also enters"
                     " the bootstrap geometric mean. Ignored unless 'enabled'.",
                     {"advanced"});
  defaults_.setValidStrings("local_rt:rt_score", {"raw", "svm", "svm_and_mbr"});

  defaultsToParam_();
}

/******************************************************************************/
PipEchoAlgorithm::~PipEchoAlgorithm() = default;

/******************************************************************************/
void PipEchoAlgorithm::group(const std::vector<FeatureMap>& feature_maps, ConsensusMap& consensus_map)
{
  PipEcho::Impl impl(param_, mz_range(feature_maps));
  PipEcho::RunMap runs;

  impl.partition_features(feature_maps, runs);

  // Ion mobility is used as a scoring feature only when EVERY run can build an
  // IM tolerance (>= 2 identified features carrying an IM width). The decision
  // is global -- all runs use IM or none do -- so the geometric-mean MBR score
  // stays on a single comparable scale: a per-run decision would mix 3- and
  // 4-feature scores in the shared SVM `mbr_score` predictor and the fallback
  // q-value ordering, which aggregate acceptors across all runs.
  auto run_supports_im = [](const PipEcho::Run& run) {
    std::size_t n = 0;
    for (const auto& donor : run.donors.storage)
    {
      if (PipEcho::Util::feature_im_width(donor->feature).has_value()) { ++n; }
    }
    return n >= 2;
  };
  const bool enable_im
    = ! runs.empty()
      && std::ranges::all_of(
        runs, [&](const auto& kv) { return run_supports_im(kv.second); });

  // [LOCAL-RT auto] Estimate the global window scales from the data before matching
  // (no-op unless 'local_rt:auto'); sets cap/floor/locality/fallback experiment-wide.
  impl.estimate_auto_params(runs);

  // [LOCAL-RT] Match donors -> acceptors. With the local adaptive window enabled,
  // wrap the O(R^2) matching in an ADAPTIVE-WIDENING loop: if too few decoy
  // transfers are generated to resolve the requested FDR, double the local window
  // scale and re-match, until estimable or the window saturates the global ceiling
  // (then the existing conservative FDR gate drops transfers). Baseline (env unset)
  // runs exactly one pass and is byte-identical.
  const std::size_t target_decoys = impl.target_decoy_count();
  for (double widen = 1.0;;)
  {
    impl.set_widen_factor(widen);

    // Clear matches from any previous (tighter) widening pass.
    for (auto& kv : runs)
    {
      for (auto& a : kv.second.acceptors.storage) { a->target.reset(); a->decoy.reset(); }
    }

    auto logger = ProgressLogger();
    std::size_t progress {};
    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, std::pow(runs.size(), 2), "matching donors and acceptors");
    for (auto& acceptor_run : runs)
    {
      PipEcho::RunStatistics stats(acceptor_run.second, enable_im);
      PipEcho::AcceptorMap& acceptors(acceptor_run.second.acceptors);

      for (auto& donor_run : runs)
      {
        if (donor_run.first != acceptor_run.first)
        {
          PipEcho::DonorMap& donors(donor_run.second.donors);
          // [LOCAL-RT] build this pair's local RT calibration first.
          impl.set_local_rt_model(donors, acceptor_run.second.donors);
          impl.link_donors_and_acceptors(stats, donors, acceptors);
        }
        logger.setProgress(++progress);
      }
    }
    logger.endProgress();

    if (! impl.local_rt_enabled()) { break; }  // baseline: single pass

    std::size_t decoys = 0;
    for (auto& kv : runs)
    {
      for (auto& a : kv.second.acceptors.storage) { if (a->decoy.has_value()) { ++decoys; } }
    }
    if (decoys >= target_decoys)
    {
      OPENMS_LOG_INFO << "[LOCAL-RT] " << decoys << " decoy transfers (>= " << target_decoys
                      << " needed) at widen x" << widen << "; FDR estimable." << std::endl;
      break;
    }
    if (widen >= impl.widen_ceiling())
    {
      OPENMS_LOG_WARN << "[LOCAL-RT] widened the local RT window to the global ceiling but only "
                      << decoys << " decoy transfer(s) (< " << target_decoys
                      << " needed); the conservative FDR gate may drop transfers." << std::endl;
      break;
    }
    widen *= 2.0;
    OPENMS_LOG_INFO << "[LOCAL-RT] only " << decoys << " decoy transfer(s) (< " << target_decoys
                    << " needed); widening local RT window x" << widen << " and re-matching." << std::endl;
  }

  impl.generate_consensus_map(runs, consensus_map);
  postprocess_(feature_maps, consensus_map);
}

} // namespace OpenMS
