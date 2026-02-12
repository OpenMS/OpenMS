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

#include "Score.h"
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
    /// A scored donor.
    using scored_t = std::pair<Score, const Donor*>;

    /// Targets.
    std::vector<scored_t> targets;

    /// Decoys.
    std::vector<scored_t> decoys;

    /// Constructor.
    Acceptor(const Peak& peak)
    : Peak(peak.map_index,
           peak.feature)
    {};

    /// Check if a Donor matches this Acceptor.
    bool is_donor_compatible(const Donor&, const Window&) const;

    /// Return the donor with the highest score.  The feature for this
    /// acceptor will be marked as either a decoy or target depending
    /// on how the donor was selected.
    std::optional<const Donor*> fetch_and_mark_best_donor();
  };
}}
