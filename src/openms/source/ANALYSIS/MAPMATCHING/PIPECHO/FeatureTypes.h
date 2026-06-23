// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "OpenMS/KERNEL/Feature.h"
#include "Score.h"
#include "Util.h"
#include "Window.h"

#include <cstddef>

namespace OpenMS::PipEcho
{

/******************************************************************************/
/**
 * Information about a Feature and which FeatureMap and mzML file it
 * came from.  This is similar to a FeatureHandle but with
 * additional information.
 */
struct FeatureRef
{
  const std::size_t map_index;
  const Feature& feature;

  FeatureRef(const std::size_t map_index, const Feature& feature):
      map_index(map_index),
      feature(feature) {};
};

/******************************************************************************/
/**
 * Used for putting FeatureRef objects in ordered containers.
 */
struct FeatureRefCmp
{
  bool operator()(const FeatureRef& lhs, const FeatureRef& rhs) const
  {
    if (lhs.map_index == rhs.map_index)
    {
      return lhs.feature.getUniqueId() < rhs.feature.getUniqueId();
    }
    return lhs.map_index < rhs.map_index;
  }
};

/******************************************************************************/
/**
 * A Feature that has been identified.
 */
struct Donor : FeatureRef
{
  Donor(const FeatureRef& peak): FeatureRef(peak.map_index, peak.feature) {};

  std::string ident() const
  {
    auto hit = Util::feature_hit(feature);

    if (hit) { return hit->getSequence().toString(); }

    std::string msg("donor feature missing peptide sequence");
    throw(Exception::MissingInformation(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, msg));
  };
};

/******************************************************************************/
/**
 * Type to make it clear if a donor feature is considered a target or
 * a decoy.
 */
enum class DonorType
{
  Target,
  Decoy
};

/******************************************************************************/
/**
 * A Feature that has no identification information and may match a
 * Donor.
 */
struct Acceptor : FeatureRef
{
  /// A scored donor.
  using scored_t = std::pair<Score, std::shared_ptr<const Donor>>;

  /// The best scoring target donor.
  std::optional<scored_t> target;

  /// The best scoring decoy donor;
  std::optional<scored_t> decoy;

  /// Constructor.
  Acceptor(const FeatureRef& peak): FeatureRef(peak.map_index, peak.feature) {};

  /// Check if a Donor matches this Acceptor.
  bool is_donor_compatible(const Donor&, const Window&) const;

  /// Add a donor.
  void push_back(const Score&, std::shared_ptr<const Donor>, DonorType);
};

} // namespace OpenMS::PipEcho
