// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h"
#include "OpenMS/CONCEPT/ProgressLogger.h"
#include "OpenMS/ML/CLUSTERING/HashGrid.h"
#include <OpenMS/CONCEPT/LogStream.h>

#include "PIPECHO/Impl.h"
#include "PIPECHO/Run.h"
#include "PIPECHO/RunStatistics.h"

namespace OpenMS
{

/******************************************************************************/
/// Find the largest m/z value among all maps.
std::pair<double, double> mz_range(const std::vector<FeatureMap>& feature_maps)
{
  double mz_min = std::numeric_limits<double>::max();
  double mz_max = std::numeric_limits<double>::lowest();

  for (auto& map : feature_maps)
  {
    // NOTE: map.getMaxMZ() always throws an exception, even if
    // updateRanges was called on the map before calling getMaxMZ.
    // Therefore we need to walk the map manually :(
    for (auto& feature : map)
    {
      mz_min = std::min(mz_min, feature.getMZ());
      mz_max = std::max(mz_max, feature.getMZ());
    }
  }

  return std::make_pair(mz_min, mz_max);
}

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

  auto logger = ProgressLogger();
  std::size_t progress {};

  logger.setLogType(ProgressLogger::CMD);
  logger.startProgress(0, std::pow(runs.size(), 2), "matching donors and acceptors");

  for (auto& acceptor_run : runs)
  {
    PipEcho::RunStatistics stats(acceptor_run.second);
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
