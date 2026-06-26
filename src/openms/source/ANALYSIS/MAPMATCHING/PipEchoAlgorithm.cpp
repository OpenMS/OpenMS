// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h>
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
                     " hit, a deterministic stratified subsample is used (all"
                     " decoys and labelled positives are kept; the rest is"
                     " score-spread sampled).",
                     {"advanced"});
  defaults_.setMinInt("max_training_points", 0);

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
        impl.link_donors_and_acceptors(stats, donors, acceptors);
      }

      logger.setProgress(++progress);
    }
  }

  logger.endProgress();
  impl.generate_consensus_map(runs, consensus_map);
  postprocess_(feature_maps, consensus_map);
}

} // namespace OpenMS
