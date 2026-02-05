// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>

#include "OpenMS/KERNEL/Feature.h"

#include "Window.h"

namespace OpenMS {
namespace PipEcho {

  /****************************************************************************/
  /**
   * Information about a Feature and which FeatureMap and mzML file it
   * came from.  This is similar to a FeatureHandle but with
   * additional information.
   */
  struct Peak {
    const std::size_t map_index;
    const Feature& feature;

    Peak(const std::size_t map_index,
         const Feature& feature)
    : map_index(map_index),
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
    : Peak(peak.map_index,
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

    /// Possible decoy Donor.
    scored_t decoy;

    /// Constructor.
    Acceptor(const Peak& peak)
    : Peak(peak.map_index,
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
  // NOTE: This method must be in the header due to type magic.
  template <typename Member>
  void Acceptor::update_donor(Member slot, const scored_t& donor) {
    if (!donor.has_value()) return;

    if (!(this->*slot).has_value() || (this->*slot)->first < donor->first) {
      this->*slot = donor;
    }
  }

  /****************************************************************************/
  // NOTE: This method must be in the header due to type magic.
  template <typename Member>
  void Acceptor::update_donor(Member slot,
                              const match_t& match,
                              const Donor& donor)
  {
    if (!match.has_value()) return;
    update_donor(slot, std::make_pair(match->first, &donor));
  }
}}
