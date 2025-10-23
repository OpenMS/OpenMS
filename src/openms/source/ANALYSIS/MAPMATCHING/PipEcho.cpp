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
#include <OpenMS/ANALYSIS/MAPMATCHING/PipEcho.h>
#include <OpenMS/ML/CLUSTERING/HashGrid.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS {

  /****************************************************************************/
  /**
   * Window used when searching for acceptor features.
   */
  struct Window {
    double rt_tol; // Seconds.
    double mz_tol; // Daltons.
    std::size_t grid_neighbors;
  };

  /****************************************************************************/
  /**
   * Information about a Feature and which FeatureMap and mzML file it
   * came from.  This is similar to a FeatureHandle but with
   * additional information.
   */
  struct Peak {
    const std::string file_name;
    const std::size_t map_index;
    const Feature& feature;

    Peak(const std::string& file_name,
         const std::size_t map_index,
         const Feature& feature)
    : file_name(file_name),
      map_index(map_index),
      feature(feature)
    { };
  };

  /****************************************************************************/
  /**
   * Used for putting Peak objects in ordered containers.
   */
  struct PeakCmp {
    bool operator()(const Peak& lhs, const Peak& rhs) const {
      if (lhs.map_index == rhs.map_index) {
        return lhs.feature.getUniqueId() < rhs.feature.getUniqueId();
      }
      return lhs.map_index < rhs.map_index;
    }
  };

  /****************************************************************************/
  /**
   * A Feature that has been identified.
   */
  struct Donor : Peak {
    Donor(const Peak& peak)
    : Peak(peak.file_name,
           peak.map_index,
           peak.feature)
    { };
  };

  /****************************************************************************/
  /**
   * A Feature that has no identification information and may match a
   * Donor.
   */
  struct Acceptor : Peak {
    /// Type used when searching for a matching Acceptor.
    using match_t = std::optional<std::pair<double, Acceptor*>>;

    /// Type used for tracking targets and decoys.
    using scored_t = std::optional<std::pair<double, const Donor*>>;

    /// Possible target Donor.
    scored_t target;

    /// Constructor.
    Acceptor(const Peak& peak)
    : Peak(peak.file_name,
           peak.map_index,
           peak.feature)
    {};

    /// Check if a Donor matches this Acceptor.
    bool is_donor_compatible(const Donor&, const Window&) const;

    /// Update a target or decoy slot.
    template <typename Member>
    void update_donor(Member, const scored_t&);

    /// Update a target or decoy slot.
    template <typename Member>
    void update_donor(Member, const match_t&, const Donor&);
  };

  /****************************************************************************/
  /**
   * This class works around an issue with the HashGrid type.
   *
   * In HashGrid, the cell contents are copied after insertion so
   * mutation is not possible.  Additionally, due to limitations in
   * its implementation it doesn't support smart pointers as cell
   * values.
   *
   * Therefore we need to use raw pointers, but need those pointers to
   * be stable and stored somewhere.  This type uses a vector of smart
   * pointers to maintain stable storage and a HashGrid for cell-based
   * access to those pointers.
   */
  template <typename T>
  struct GridWithStorage {
    // Handy aliases.
    using grid_t = HashGrid<T*>;
    using grid_center_t = grid_t::ClusterCenter;
    using grid_index_t = grid_t::CellIndex;

    // Access to the underlying storage.
    std::vector<std::unique_ptr<T>> storage;
    grid_t grid;

    /// Construct a new object given the grid center.
    GridWithStorage(const grid_center_t& center)
    : grid(center) { };

    /// Insert a new object into the storage and grid.
    void insert(const Peak& peak) {
      grid_center_t center(peak.feature.getRT(), peak.feature.getMZ());
      auto acceptor = std::make_unique<T>(peak);
      storage.push_back(std::move(acceptor));
      grid.insert(std::make_pair(center, storage.back().get()));
    }

    /// Remove all objects from the storage and grid.
    void clear() {
      grid.clear();
      storage.clear();
    }
  };

  /****************************************************************************/
  // Handy aliases.
  using DonorMap = std::vector<Donor>;
  using AcceptorMap = GridWithStorage<Acceptor>;

  /****************************************************************************/
  /**
   * Internal representation of the PIP-ECHO algorithm.
   */
  class PipEchoImpl {
  public:
    /// Construct a new implementation object.
    PipEchoImpl(const Param& params);

    /// Separate donors from acceptors.
    void partition_features(const std::vector<FeatureMap>&,
                            DonorMap&, AcceptorMap&);

    /// Match identified features with unidentified features.
    void link_donors_and_acceptors(const DonorMap&, AcceptorMap&);

    /// Search for a matching acceptor.
    Acceptor::match_t
    find_acceptor_for(const AcceptorMap&, const Donor&, const Window&);

    /// Fill in the final ConsensusMap.
    void generate_consensus_map(DonorMap&, AcceptorMap&, ConsensusMap&);
  public:
    /// Max allowed m/z difference between donor and acceptor.
    double mz_dal_max_diff;

    /// Max allowed RT difference between donor and acceptor.
    double rt_sec_max_window;

  private:
    std::string path_from_feature_map(const FeatureMap&);
    bool is_donor_feature(const Feature&);

    std::optional<Window> initial_window();
    std::optional<Window> next_window(const std::optional<Window>&);
  };

  /****************************************************************************/
  /**
   * Return the first peptide hit from a feature.
   *
   * This code is here to isolate checking hits because the process
   * might change in a future version of OpenMS.
   */
  std::optional<PeptideHit> feature_hit(const Feature& feature) {
    const PeptideIdentificationList& peps = feature.getPeptideIdentifications();

    // FIXME: Is this the correct way to know if a feature has an
    // ID?  There seems to be two ways to store an ID on a feature,
    // but most code I see uses the "old way".
    if (!peps.empty() && !peps[0].getHits().empty()) {
      return peps[0].getHits()[0];
    }

    return {};
  }

  /****************************************************************************/
  /**
   * Return a key that can be used to link features with the same
   * amino acid sequence and charge state.
   */
  std::string feature_sequence_key(const Feature& feature) {
    auto hit = feature_hit(feature);

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
  PipEcho::PipEcho()
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
  PipEcho::~PipEcho() = default;

  /****************************************************************************/
  /// Execute the PIP-ECHO algorithm.
  void PipEcho::group(
    const std::vector<FeatureMap>& feature_maps,
    ConsensusMap& consensus_map
  )
  {
    PipEchoImpl impl(param_);
    AcceptorMap acceptors({impl.rt_sec_max_window, impl.mz_dal_max_diff});
    DonorMap donors;

    impl.partition_features(feature_maps, donors, acceptors);
    impl.link_donors_and_acceptors(donors, acceptors);
    impl.generate_consensus_map(donors, acceptors, consensus_map);

    postprocess_(feature_maps, consensus_map);
  }

  /****************************************************************************/
  PipEchoImpl::PipEchoImpl(const Param& params)
  : mz_dal_max_diff(params.getValue("distance_MZ:max_difference")),
    rt_sec_max_window(params.getValue("distance_RT:max_difference"))
  {
    if (params.getValue("distance_MZ:unit") == "ppm") {
      // FIXME: implement this!
      throw(Exception::NotImplemented(__FILE__,
                                      __LINE__,
                                      OPENMS_PRETTY_FUNCTION));
    }

    OPENMS_LOG_INFO << "PIP-ECHO("
                    << rt_sec_max_window << ", "
                    << mz_dal_max_diff << ")"
                    << std::endl;
  }

  /****************************************************************************/
  /// Separate donors and acceptors.
  void PipEchoImpl::partition_features(const std::vector<FeatureMap>& maps,
                                       DonorMap& donors,
                                       AcceptorMap& acceptors)
  {
    std::size_t num_maps = maps.size();

    for (std::size_t i=0; i<num_maps; ++i) {
      const FeatureMap& map = maps[i];
      std::string file_name = path_from_feature_map(map);

      for (auto& feature : map) {
        // We need to strip off the const-ness of the Feature so
        // prepare it for use.  We don't really care about it being
        // const, only that we have a reference/pointer.  But a
        // FeatureMap doesn't let us have a non-const reference
        // iterator.  So, we have to cheat and cast it.
        prepare_feature(const_cast<Feature&>(feature));
        Peak peak(file_name, i, feature);

        if (is_donor_feature(feature)) {
          donors.push_back(Donor(peak));
        } else {
          acceptors.insert(peak);
        }
      }
    }
  }

  /****************************************************************************/
  /// Link matching donors and acceptors.
  void PipEchoImpl::link_donors_and_acceptors(const DonorMap& donors,
                                              AcceptorMap& acceptors)
  {
    for (auto& donor : donors) {
      for (auto window=initial_window(); window.has_value();
           window = next_window(window))
        {
          auto target = find_acceptor_for(acceptors, donor, *window);
          if (target.has_value()) {
            target->second->update_donor(&Acceptor::target, target, donor);
          }

          if (target) {
            // FIXME: create decoys and open window only when there
            // are no decoys and no targets.
            break;
          }
        }
    }
  }

  /****************************************************************************/
  /// Find an acceptor for the given donor.
  Acceptor::match_t
  PipEchoImpl::find_acceptor_for(const AcceptorMap& acceptors,
                                 const Donor& donor,
                                 const Window& window)
  {
    Acceptor::match_t best_match;

    // FIXME: Remove these later.
    std::random_device randev;
    std::mt19937 randgen(randev());

    auto check_bucket = [&](AcceptorMap::grid_t::const_grid_iterator bucket) {
      for (auto& acceptor : bucket->second) {
        if (acceptor.second->is_donor_compatible(donor, window)) {
          // FIXME: Fake score given here!
          double score = std::generate_canonical<double, 10>(randgen);

          if (!best_match.has_value() || best_match->first < score) {
            best_match = {score, acceptor.second};
          }
        }
      }
    };

    AcceptorMap::grid_center_t center(donor.feature.getRT(), donor.feature.getMZ());
    const AcceptorMap::grid_index_t index = acceptors.grid.cellIndexAtClusterCenter(center);

    // Check the cell where we expect to find the acceptor.
    if (auto cell = acceptors.grid.grid_find(index);
        cell != acceptors.grid.grid_end())
      {
        check_bucket(cell);
      }

    // Check nearby cells.
    const auto index_x = index.getX();
    const auto index_y = index.getY();
    const auto offset = window.grid_neighbors;

    for (auto i = index_x - offset; i <= index_x + offset; ++i) {
      for (auto j = index_y - offset; j <= index_y + offset; ++j) {
        const auto addr = AcceptorMap::grid_index_t(i, j);

        if (auto cell = acceptors.grid.grid_find(addr);
            cell != acceptors.grid.grid_end())
          {
            check_bucket(cell);
          }
      }
    }

    return best_match;
  }

  /****************************************************************************/
  /// Turn our internal data structures into a ConsensusMap.
  void PipEchoImpl::generate_consensus_map(DonorMap& donors,
                                           AcceptorMap& acceptors,
                                           ConsensusMap& consensus_map)
  {
    using IdentVal = std::set<Peak, PeakCmp>;
    using IdentMap = std::map<std::string, IdentVal>;
    IdentMap ident_map;

    auto insert = [&](const Donor& donor, const Peak& peak) {
      std::string key = feature_sequence_key(donor.feature);
      auto [place, inserted] = ident_map.insert(std::make_pair(key, IdentVal{peak}));
      if (!inserted) place->second.insert(peak);
    };

    // Place each donor into the identity map.
    for (auto& donor : donors) {
      insert(donor, donor);
    }

    std::size_t acceptors_added{};

    // Now place the acceptors in there too.
    for (auto& acceptor : acceptors.storage) {
      if (acceptor->target.has_value()) {
        ++acceptors_added;
        insert(*acceptor->target->second, *acceptor);
      }
    }

    OPENMS_LOG_INFO
      << "PIP-ECHO: added "
      << acceptors_added << " acceptors to the consensus_map"
      << std::endl;

    // Don't need this anymore.
    acceptors.clear();
    donors.clear();

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
  /// Return `true` if the given Donor is compatible with this Acceptor.
  bool Acceptor::is_donor_compatible(const Donor& donor, const Window& window) const {
    return donor.file_name != this->file_name &&
           donor.feature.getCharge() == this->feature.getCharge() &&
           std::fabs(donor.feature.getRT() - this->feature.getRT()) <= window.rt_tol &&
           std::fabs(donor.feature.getMZ() - this->feature.getMZ()) <= window.mz_tol;
  }

  /****************************************************************************/
  /// Update the `target` or `decoy` members if the given Donor has a
  /// higher score.
  template <typename Member>
  void Acceptor::update_donor(Member slot, const Acceptor::scored_t& donor) {
    if (!donor.has_value()) return;

    if (!(this->*slot).has_value() || (this->*slot)->first < donor->first) {
      this->*slot = donor;
    }
  }

  /****************************************************************************/
  /// Helper to convert from a found Acceptor to a scored Donor.
  template <typename Member>
  void Acceptor::update_donor(Member slot,
                              const Acceptor::match_t& match,
                              const Donor& donor)
  {
    if (!match.has_value()) return;
    this->update_donor(slot, std::make_pair(match->first, &donor));
  }

  /****************************************************************************/
  /// Extract a single file path from a FeatureMap.
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
  /// Return `true` if the given Feature is a donor.
  bool PipEchoImpl::is_donor_feature(const Feature& feature) {
    return feature_hit(feature).has_value();
  }

  /****************************************************************************/
  /// Return a starting window for feature searching.
  std::optional<Window> PipEchoImpl::initial_window() {
    return Window{
      rt_sec_max_window,
      mz_dal_max_diff,
      1};
  }

  /****************************************************************************/
  /// Increment the given window.
  std::optional<Window> PipEchoImpl::next_window(const std::optional<Window>& prev) {
    if (prev && prev->rt_tol <= rt_sec_max_window && prev->grid_neighbors < 3) {
      Window next(*prev);
      next.grid_neighbors++;
      return next;
    } else {
      return {};
    }
  }
}
