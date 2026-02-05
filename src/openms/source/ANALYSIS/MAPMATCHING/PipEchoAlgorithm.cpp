// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <memory>
#include <random>
#include <set>

#include "OpenMS/CONCEPT/Exception.h"
#include "OpenMS/config.h"
#include "OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h"
#include "OpenMS/CONCEPT/ProgressLogger.h"
#include "OpenMS/ML/CLUSTERING/HashGrid.h"
#include <OpenMS/CONCEPT/LogStream.h>

#include "PIPECHO/GridWithStorage.h"
#include "PIPECHO/MzDiff.h"
#include "PIPECHO/PeakTypes.h"
#include "PIPECHO/Run.h"
#include "PIPECHO/Util.h"
#include "PIPECHO/Window.h"

namespace OpenMS {

  /****************************************************************************/
  using OpenMS::PipEcho::Acceptor;
  using OpenMS::PipEcho::Donor;
  using OpenMS::PipEcho::GridWithStorage;
  using OpenMS::PipEcho::MzDiff;
  using OpenMS::PipEcho::Peak;
  using OpenMS::PipEcho::PeakCmp;
  using OpenMS::PipEcho::Run;
  using OpenMS::PipEcho::Window;

  /****************************************************************************/
  // Handy aliases.
  using DonorMap = GridWithStorage<Donor>;
  using AcceptorMap = GridWithStorage<Acceptor>;
  using RunMap = std::map<std::string, Run>;

  /****************************************************************************/
  /**
   * Internal representation of the PIP-ECHO algorithm.
   */
  class PipEchoImpl {
  public:
    /// Construct a new implementation object.
    PipEchoImpl(const Param& params, const std::pair<double, double>& mz_range);

    /// Separate donors from acceptors.
    void partition_features(const std::vector<FeatureMap>&,
                            RunMap&);

    /// Match identified features with unidentified features.
    void link_donors_and_acceptors(const DonorMap&, AcceptorMap&);

    /// Search for a matching acceptor.
    Acceptor::match_t
    find_acceptor_for(const AcceptorMap&,
                      const Donor&,
                      const Window&,
                      const std::optional<double> = {});

    /// Find a random Donor that is dissimilar to a given Donor.
    std::optional<const Donor*>
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

  private:
    std::string path_from_feature_map(const FeatureMap&);
    bool is_donor_feature(const Feature&);
    Run& get_run_from_file_name(RunMap&, const std::string&);

    std::optional<Window> initial_window(const Donor& donor);
    std::optional<Window> next_window(const std::optional<Window>&);
  };

  /****************************************************************************/
  /**
   * Return a key that can be used to link features with the same
   * amino acid sequence and charge state.
   */
  std::string feature_sequence_key(const Feature& feature) {
    auto hit = PipEcho::Util::feature_hit(feature);

    if (hit) {
      return hit->getSequence().toString() + "_" + feature.getCharge();
    }

    std::string msg("donor feature missing peptide sequence");

    throw(Exception::MissingInformation(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        msg));
  }

  /****************************************************************************/
  /// Ensure that a feature is ready for analysis.
  void prepare_feature(Feature& feature) {
    // Ensure everything is sorted for later use.
    for (auto& peptide : feature.getPeptideIdentifications()) {
      peptide.sort();
    }

    feature.sortPeptideIdentifications();
  }

  /****************************************************************************/
  /// Find the largest m/z value among all maps.
  std::pair<double,double> mz_range(const std::vector<FeatureMap>& feature_maps) {
    double mz_min = std::numeric_limits<double>::max();
    double mz_max = std::numeric_limits<double>::lowest();

    for (auto &map : feature_maps) {
      // NOTE: map.getMaxMZ() always throws an exception, even if
      // updateRanges was called on the map before calling getMaxMZ.
      // Therefore we need to walk the map manually :(
      for (auto& feature : map) {
        mz_min = std::min(mz_min, feature.getMZ());
        mz_max = std::max(mz_max, feature.getMZ());
      }
    }

    return std::make_pair(mz_min, mz_max);
  }

  /****************************************************************************/
  PipEchoAlgorithm::PipEchoAlgorithm()
  : FeatureGroupingAlgorithm()
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

    defaults_.setValue("distance_MZ:unit", "Da",
                       "Unit of the 'max_difference' parameter");
    defaults_.setValidStrings("distance_MZ:unit", {"Da","ppm"});

    defaultsToParam_();
  }

  /****************************************************************************/
  PipEchoAlgorithm::~PipEchoAlgorithm() = default;

  /****************************************************************************/
  void PipEchoAlgorithm::group(
    const std::vector<FeatureMap>& feature_maps,
    ConsensusMap& consensus_map
  )
  {
    PipEchoImpl impl(param_, mz_range(feature_maps));
    RunMap runs;

    impl.partition_features(feature_maps, runs);

    auto logger = ProgressLogger();
    std::size_t progress{};

    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, std::pow(runs.size(), 2),
                         "matching donors and acceptors");

    for (auto& acceptor_run : runs) {
      AcceptorMap& acceptors = acceptor_run.second.acceptors;

      for (auto& donor_run : runs) {
        if (donor_run.first != acceptor_run.first) {
          DonorMap& donors = donor_run.second.donors;
          impl.link_donors_and_acceptors(donors, acceptors);
        }

        logger.setProgress(++progress);
      }
    }

    logger.endProgress();
    impl.generate_consensus_map(runs, consensus_map);
    postprocess_(feature_maps, consensus_map);
  }

  /****************************************************************************/
  PipEchoImpl::PipEchoImpl(const Param& params, const std::pair<double, double>& mz_range)
  : mz_max_diff(params),
    mz_grid_center(0.5), // FIXME: Can this be calculated?
    //mz_grid_center(mz_range.second * 1e-6 * 10.0), // 10 PPM of max m/z
    rt_sec_max_window(params.getValue("distance_RT:max_difference"))
  {
    OPENMS_LOG_INFO << "PIP-ECHO("
                    << mz_range.first << ", "
                    << mz_range.second << ", "
                    << mz_grid_center << ", "
                    << mz_max_diff.mz_diff(mz_range.second - mz_range.first) << ", "
                    << rt_sec_max_window << ")"
                    << std::endl;
  }

  /****************************************************************************/
  void PipEchoImpl::partition_features(const std::vector<FeatureMap>& maps,
                                       RunMap& runs)
  {
    std::size_t num_maps = maps.size();

    auto logger = ProgressLogger();
    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, num_maps, "building maps");

    for (std::size_t i=0; i<num_maps; ++i) {
      logger.setProgress(i);

      const FeatureMap& map = maps[i];
      std::string file_name = path_from_feature_map(map);
      Run& run = get_run_from_file_name(runs, file_name);

      for (auto& feature : map) {
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

  /****************************************************************************/
  void PipEchoImpl::link_donors_and_acceptors(const DonorMap& donors,
                                              AcceptorMap& acceptors)
  {
    for (auto& donor : donors.storage) {
      for (auto window=initial_window(*donor); window.has_value();
           window = next_window(window))
        {
          const auto target = find_acceptor_for(acceptors, *donor, *window);
          if (target.has_value()) {
            target->second->update_donor(&Acceptor::target, target, *donor);
          }

          const auto random_donor = find_random_donor(donors, *donor, *window);
          Acceptor::match_t decoy = {}; // No std::optional::and_then in C++ 20 :(

          if (random_donor.has_value()) {
            decoy = find_acceptor_for(acceptors,
                                      *donor,
                                      *window,
                                      (*random_donor)->feature.getRT());

            if (decoy.has_value()) {
              decoy->second->update_donor(&Acceptor::decoy, decoy, *donor);
            }
          }

          if (target.has_value() || decoy.has_value()) {
            break; // We can stop searching.
          }
        }
    }
  }

  /****************************************************************************/
  Acceptor::match_t
  PipEchoImpl::find_acceptor_for(const AcceptorMap& acceptors,
                                 const Donor& donor,
                                 const Window& window,
                                 const std::optional<double> rt_override)
  {
    // FIXME: Remove these later.
    std::random_device randev;
    std::mt19937 randgen(randev());

    // FIXME: What does the paper say about finding acceptor peaks?
    // In FlashLFQ they do some weird envelope cutting and don't use
    // the actual found peak (FindIndividualAcceptorPeak).
    const double rt = rt_override.value_or(donor.feature.getRT());
    const AcceptorMap::grid_center_t center(rt, donor.feature.getMZ());
    const AcceptorMap::grid_index_t index = acceptors.grid.cellIndexAtClusterCenter(center);

    return acceptors.nearby(
      index, window.grid_neighbors, Acceptor::match_t{},
      [&](Acceptor::match_t best_match, Acceptor& acceptor) -> Acceptor::match_t {
        if (acceptor.is_donor_compatible(donor, window)) {
          // FIXME: Fake score given here!
          double score = std::generate_canonical<double, 10>(randgen);

          if (!best_match.has_value() || best_match->first < score) {
            return std::make_pair(score, &acceptor);
          }
        }

        return best_match;
      });
  }

  /****************************************************************************/
  std::optional<const Donor*>
  PipEchoImpl::find_random_donor(const DonorMap& donors,
                                 const Donor& start,
                                 const Window& window) const
  {
    /// FIXME: Is this mass calculation good enough?
    const auto mass = [](const Donor& donor) -> double {
      return (donor.feature.getMZ() - Constants::PROTON_MASS_U) *
             donor.feature.getCharge();
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
    const double rt_min_diff = window.rt_tol * 2.0;  // FIXME: Is this okay?

    const double start_mass = mass(start);
    const double start_rt = start.feature.getRT();
    const std::string start_seq = base_seq(start);

    std::vector<const Donor*> matching_donors;

    auto explore = [&](double distance, Donor& donor) -> double {
      const double mass_diff = std::fabs(start_mass - mass(donor));
      const double rt_diff = std::fabs(start_rt - donor.feature.getRT());

      if (rt_diff > rt_min_diff &&
          mass_diff > mass_min_diff &&
          mass_diff < mass_max_diff &&
          start_seq != base_seq(donor))
        {
          matching_donors.push_back(&donor);
        }

      return std::max(distance, mass_diff);
    };

    double cur_max_mass_diff = mass_max_diff;
    std::size_t neighbors = window.grid_neighbors;

    const DonorMap::grid_center_t center(start.feature.getRT(), start.feature.getMZ());
    const DonorMap::grid_index_t index = donors.grid.cellIndexAtClusterCenter(center);

    do {
      const double distance = donors.nearby(index, neighbors, 0, explore);

      // Adjust search params for next round.
      if (distance < cur_max_mass_diff) {
        // We're not exploring far enough.
        ++neighbors;
      } else {
        // We need to increase our tolerance.
        cur_max_mass_diff *= 10;
      }
    } while (matching_donors.empty() &&
             cur_max_mass_diff < 1e5 &&
             neighbors < 5);

    if (!matching_donors.empty()) {
      std::srand(std::time(nullptr));
      return matching_donors[std::rand() % matching_donors.size()];
    }

    return {};
  }

  /****************************************************************************/
  void PipEchoImpl::generate_consensus_map(RunMap& runs, ConsensusMap& consensus_map) {
    using IdentVal = std::set<Peak, PeakCmp>;
    using IdentMap = std::map<std::string, IdentVal>;

    std::size_t acceptors_added{}, decoys_seen{};
    IdentMap ident_map;

    auto insert = [&](const Donor& donor, const Peak& peak) {
      std::string key = feature_sequence_key(donor.feature);
      auto [place, inserted] = ident_map.insert(std::make_pair(key, IdentVal{peak}));
      if (!inserted) place->second.insert(peak);
    };

    for (auto& run : runs) {
      // Place each donor into the identity map.
      for (auto& donor : run.second.donors.storage) {
        insert(*donor, *donor);
      }

      // Now place the acceptors in there too.
      for (auto& acceptor : run.second.acceptors.storage) {
        if (acceptor->target.has_value()) {
          ++acceptors_added;
          insert(*acceptor->target->second, *acceptor);
        }
        if (acceptor->decoy.has_value()) ++decoys_seen;
      }

      // Don't need this anymore.
      //run.second.clear(); FIXME: Ug, this results in a SIGSEV :(
    }

    OPENMS_LOG_INFO
      << "PIP-ECHO: added "
      << acceptors_added << " acceptors to the consensus_map with "
      << decoys_seen << " decoys seen"
      << std::endl;

    // We can now turn that IdentMap into a ConsensusMap.
    for (auto& group : ident_map) {
      ConsensusFeature consensus_feature;

      for (auto& peak : group.second) {
        // Note: This is apparently the correct way to add a feature
        // to a ConsensusFeature.  It does more work, such as copying
        // identifications, that the other insert methods don't do.
        consensus_feature.insert(peak.map_index, peak.feature);
      }

      consensus_feature.computeConsensus();
      consensus_map.push_back(consensus_feature);
    }
  }

  /****************************************************************************/
  std::string PipEchoImpl::path_from_feature_map(const FeatureMap& map) {
    StringList paths;
    map.getPrimaryMSRunPath(paths);

    if (paths.empty()) {
      OPENMS_LOG_WARN << "Warning: feature map has no MS run path!" << std::endl;
      return "unknown";
    } else {
      if (paths.size() > 1) {
        OPENMS_LOG_WARN << "Warning: feature map has more than one MS run path, "
                        << "using first path: " << paths[0] << std::endl;
      }

      return paths[0];
    }
  }

  /****************************************************************************/
  Run& PipEchoImpl::get_run_from_file_name(RunMap& runs, const std::string& file_name) {
      auto run_it = runs.find(file_name);

      if (run_it == runs.end()) {
        const auto [it, status] =
          runs.try_emplace(file_name,
                           file_name,
                           rt_sec_max_window,
                           mz_grid_center);

        assert(status);
        run_it = it;
      }

      return run_it->second;;
  }

  /****************************************************************************/
  std::optional<Window> PipEchoImpl::initial_window(const Donor& donor) {
    double mz_diff = mz_max_diff.mz_diff(donor.feature.getMZ());

    return Window{
      rt_sec_max_window,
      mz_diff,
      1};
  }

  /****************************************************************************/
  std::optional<Window> PipEchoImpl::next_window(const std::optional<Window>& prev) {
    if (prev && prev->rt_tol <= rt_sec_max_window && prev->grid_neighbors < 3) {
      Window next(*prev);
      next.grid_neighbors++;
      return next;
    } else {
      return {};
    }
  }

  /****************************************************************************/
  // double Score::calc_intensity_score(const Feature& donor, const Feature& acceptor) {
  //   double donor_intensity = donor.getIntensity();
  //   double acceptor_intensity = acceptor.getIntensity();

//     // var logIntensity = Math.Log(acceptorIntensity, 2);
  //   // return CalculateScore(_logIntensityDistribution, logIntensity);
  // }
}
