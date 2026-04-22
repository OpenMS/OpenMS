// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "Impl.h"
#include "OpenMS/CONCEPT/Exception.h"
#include "OpenMS/CONCEPT/ProgressLogger.h"
#include "OpenMS/KERNEL/FeatureMap.h"
#include "PeakTypes.h"
#include "Pep.h"
#include "Score.h"
#include "Util.h"

#include <memory>
#include <ranges>
#include <set>

namespace OpenMS::PipEcho
{

/******************************************************************************/
/**
 * Return a key that can be used to link features with the same
 * amino acid sequence and charge state.
 */
std::string feature_sequence_key(const std::string& ident,
                                 BaseFeature::ChargeType charge)
{
  return ident + "_" + std::to_string(charge);
}

/******************************************************************************/
/**
 * Return a key that can be used to link features with the same
 * amino acid sequence and charge state.
 */
std::string feature_sequence_key(const Feature& feature)
{
  auto hit = PipEcho::Util::feature_hit(feature);

  if (hit)
  {
    return feature_sequence_key(hit->getSequence().toString(),
                                feature.getCharge());
  }

  std::string msg("donor feature missing peptide sequence");
  throw(Exception::MissingInformation(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION, msg));
}

/******************************************************************************/
/// Ensure that a feature is ready for analysis.
void prepare_feature(Feature& feature)
{
  // Ensure everything is sorted for later use.
  for (auto& peptide : feature.getPeptideIdentifications())
  {
    peptide.sort();
  }

  feature.sortPeptideIdentifications();
}

/******************************************************************************/
Impl::Impl(const Param& params, const std::pair<double, double>& mz_range):
    mz_max_diff(params),
    mz_grid_center(0.5),
    rt_sec_max_window(params.getValue("distance_RT:max_difference")),
    mbr_fdr(params.getValue("fdr"))
{
  OPENMS_LOG_INFO << "MBR via PIP-ECHO(" << mz_range.first << ", "
                  << mz_range.second << ", " << mz_grid_center << ", "
                  << mz_max_diff.mz_diff(mz_range.second - mz_range.first)
                  << ", " << rt_sec_max_window << ", " << mbr_fdr << ")"
                  << std::endl;
}

/******************************************************************************/
void Impl::partition_features(const std::vector<FeatureMap>& maps, RunMap& runs)
{
  std::size_t num_maps = maps.size();

  auto logger = ProgressLogger();
  logger.setLogType(ProgressLogger::CMD);
  logger.startProgress(0, num_maps, "building maps");

  for (std::size_t i = 0; i < num_maps; ++i)
  {
    logger.setProgress(i);

    const FeatureMap& map(maps[i]);
    std::string file_name(path_from_feature_map(map));
    Run& run = get_run_from_file_name(runs, file_name);

    for (auto& feature : map)
    {
      // We need to strip off the const-ness of the Feature to
      // prepare it for use.  We don't really care about it being
      // const, only that we have a reference/pointer.  But a
      // FeatureMap doesn't let us have a non-const reference
      // iterator.  So, we have to cheat and cast it.
      prepare_feature(const_cast<Feature&>(feature));
      run.insert(feature, i);
    }
  }

  logger.endProgress();
}

/******************************************************************************/
void Impl::link_donors_and_acceptors(const RunStatistics& stats,
                                     const DonorMap& donors,
                                     AcceptorMap& acceptors)
{
  for (auto& donor : donors.storage)
  {
    for (auto window = initial_window(*donor); window.has_value();
         window = next_window(window))
    {
      const auto target = find_acceptor_for(stats, acceptors, *donor, *window);

      if (target.has_value())
      {
        auto acceptor = target->second;
        acceptor->push_back(target->first, donor, DonorType::Target);
      }

      const auto random_donor = find_random_donor(donors, *donor, *window);
      match_t decoy = {}; // No std::optional::and_then in C++ 20 :(

      if (random_donor.has_value())
      {
        decoy = find_acceptor_for(stats, acceptors, *donor, *window,
                                  (*random_donor)->feature.getRT());

        if (decoy.has_value())
        {
          auto acceptor = decoy->second;
          acceptor->push_back(decoy->first, donor, DonorType::Decoy);
        }
      }

      if (target.has_value() || decoy.has_value())
      {
        break; // We can stop searching.
      }
    }
  }
}

/******************************************************************************/
Impl::match_t Impl::find_acceptor_for(const RunStatistics& stats,
                                      const AcceptorMap& acceptors,
                                      const Donor& donor,
                                      const Window& window,
                                      const std::optional<double> rt_override)
{
  // FIXME: What does the paper say about finding acceptor peaks?
  // In FlashLFQ they do some weird envelope cutting and don't use
  // the actual found peak (FindIndividualAcceptorPeak).
  const double rt = rt_override.value_or(donor.feature.getRT());
  const AcceptorMap::grid_center_t center(rt, donor.feature.getMZ());
  const AcceptorMap::grid_index_t index
    = acceptors.grid.cellIndexAtClusterCenter(center);

  return acceptors.nearby<match_t>(
    index, window.grid_neighbors, match_t {},
    [&](match_t best_match, AcceptorMap::value_type acceptor) -> match_t {
      if (acceptor->is_donor_compatible(donor, window))
      {
        Score score(stats.score(donor.feature, acceptor->feature));

        if (! best_match.has_value()
            || best_match->first.mbr_score < score.mbr_score)
        {
          return std::make_pair(score, acceptor);
        }
      }

      return best_match;
    });
}

/******************************************************************************/
std::optional<std::shared_ptr<const Donor>>
Impl::find_random_donor(const DonorMap& donors,
                        const Donor& start,
                        const Window& window) const
{
  /// FIXME: Is this mass calculation good enough?
  const auto mass = [](const Donor& donor) -> double {
    return (donor.feature.getMZ() - Constants::PROTON_MASS_U)
           * donor.feature.getCharge();
  };

  /// FIXME: Is this sequence calculation good enough?
  /// FIXME: What is a "Base Sequence"?
  const auto base_seq = [](const Donor& donor) -> std::string {
    const auto hit(PipEcho::Util::feature_hit(donor.feature));
    assert(hit.has_value());
    return hit->getSequence().toUnmodifiedString();
  };

  const double mass_min_diff = 5.0 * Constants::PROTON_MASS_U;
  const double mass_max_diff = 11.0 * Constants::PROTON_MASS_U;
  const double rt_min_diff = window.rt_tol * 2.0; // FIXME: Is this okay?

  const double start_mass = mass(start);
  const double start_rt = start.feature.getRT();
  const std::string start_seq = base_seq(start);

  std::vector<std::shared_ptr<const Donor>> matching_donors;

  auto explore
    = [&](double distance, std::shared_ptr<const Donor> donor) -> double {
    const double mass_diff = std::fabs(start_mass - mass(*donor));
    const double rt_diff = std::fabs(start_rt - donor->feature.getRT());

    if (rt_diff > rt_min_diff && mass_diff > mass_min_diff
        && mass_diff < mass_max_diff && start_seq != base_seq(*donor))
    {
      matching_donors.push_back(donor);
    }

    return std::max(distance, mass_diff);
  };

  double cur_max_mass_diff = mass_max_diff;
  std::size_t neighbors = window.grid_neighbors;

  const DonorMap::grid_center_t center(start.feature.getRT(),
                                       start.feature.getMZ());
  const DonorMap::grid_index_t index
    = donors.grid.cellIndexAtClusterCenter(center);

  do
  {
    const double distance = donors.nearby<double>(index, neighbors, 0, explore);

    // Adjust search params for next round.
    if (distance < cur_max_mass_diff)
    {
      // We're not exploring far enough.
      ++neighbors;
    }
    else
    {
      // We need to increase our tolerance.
      cur_max_mass_diff *= 10;
    }
  } while (matching_donors.empty() && cur_max_mass_diff < 1e5 && neighbors < 5);

  if (! matching_donors.empty())
  {
    std::srand(std::time(nullptr));
    return matching_donors[std::rand() % matching_donors.size()];
  }

  return {};
}

/******************************************************************************/
void Impl::generate_consensus_map(RunMap& runs, ConsensusMap& consensus_map)
{
  using ident_val_t = std::set<Peak, PeakCmp>;
  using ident_map_t = std::map<std::string, ident_val_t>;

  std::size_t acceptors_seen {}, decoys_seen {};
  ident_map_t ident_map;

  auto insert = [&](const std::string key, const Peak& peak) {
    auto [place, inserted]
      = ident_map.insert(std::make_pair(key, ident_val_t {peak}));
    if (! inserted) place->second.insert(peak);
  };

  for (auto& run : runs)
  {
    // Place each donor into the identity map.
    for (auto& donor : run.second.donors.storage)
    {
      std::string key = feature_sequence_key(donor->feature);
      insert(key, *donor);
    }
  }

  auto view = runs | std::views::transform([](auto& run) {
                return run.second.acceptors.storage;
              })
              | std::views::join;

  std::vector<std::shared_ptr<Acceptor>> all_acceptors;

  //  FIXME
  // all_acceptors.reserve(std::ranges::size(view));
  // all_acceptors.insert(std::ranges::begin(view),
  // std::ranges::end(view));

  for (auto& acpt : view)
  {
    all_acceptors.push_back(acpt);
  }

  PipEcho::Pep pep(all_acceptors);
  const PipEcho::Pep::group_t& acceptors = pep.run(mbr_fdr);

  // Put each acceptor into the correct feature bucket.
  for (auto& acceptor : acceptors)
  {
    std::string key
      = feature_sequence_key(acceptor->donor_ident, acceptor->donor_charge);

    insert(key, acceptor->acceptor);

    ++acceptors_seen;
    if (! acceptor->is_target()) ++decoys_seen;
  }

  OPENMS_LOG_INFO << "PIP-ECHO: added " << acceptors_seen
                  << " acceptors to the consensus map with " << decoys_seen
                  << " decoys seen" << std::endl;

  // We can now turn that IdentMap into a ConsensusMap.
  for (auto& group : ident_map)
  {
    ConsensusFeature consensus_feature;

    for (auto& peak : group.second)
    {
      // Note: This is apparently the correct way to add a feature
      // to a ConsensusFeature.  It does more work, such as copying
      // identifications, that the other insert methods don't do.
      consensus_feature.insert(peak.map_index, peak.feature);
    }

    consensus_feature.computeConsensus();
    consensus_map.push_back(consensus_feature);
  }
}

/******************************************************************************/
std::string Impl::path_from_feature_map(const FeatureMap& map)
{
  StringList paths;
  map.getPrimaryMSRunPath(paths);

  if (paths.empty())
  {
    OPENMS_LOG_WARN << "Warning: feature map has no MS run path!" << std::endl;
    return "unknown";
  }
  else
  {
    if (paths.size() > 1)
    {
      OPENMS_LOG_WARN << "Warning: feature map has more than one MS run path, "
                      << "using first path: " << paths[0] << std::endl;
    }

    return paths[0];
  }
}

/******************************************************************************/
Run& Impl::get_run_from_file_name(RunMap& runs, const std::string& file_name)
{
  auto run_it = runs.find(file_name);

  if (run_it == runs.end())
  {
    const auto [it, status]
      = runs.try_emplace(file_name, rt_sec_max_window, mz_grid_center);

    assert(status);
    run_it = it;
  }

  return run_it->second;
  ;
}

/******************************************************************************/
std::optional<Window> Impl::initial_window(const Donor& donor)
{
  double mz_diff = mz_max_diff.mz_diff(donor.feature.getMZ());

  return Window {rt_sec_max_window, mz_diff, 1};
}

/******************************************************************************/
std::optional<Window> Impl::next_window(const std::optional<Window>& prev)
{
  if (prev && prev->rt_tol <= rt_sec_max_window && prev->grid_neighbors < 3)
  {
    Window next(*prev);
    next.grid_neighbors++;
    return next;
  }
  else
  {
    return {};
  }
}

} // namespace OpenMS::PipEcho
