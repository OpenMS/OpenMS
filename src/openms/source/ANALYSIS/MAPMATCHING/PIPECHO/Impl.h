// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
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

#include <random>

namespace OpenMS::PipEcho
{

/******************************************************************************/
// Handy aliases.
using DonorMap = GridWithStorage<Donor>;
using AcceptorMap = GridWithStorage<Acceptor>;
using RunMap = std::map<std::string, Run>;

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
  mutable std::mt19937 rng_;
};

} // namespace OpenMS::PipEcho
