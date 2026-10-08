// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "FeatureTypes.h"
#include "GridWithStorage.h"
#include <OpenMS/KERNEL/Feature.h>

#include <memory>

namespace OpenMS::PipEcho
{

/******************************************************************************/
class Run
{
public:
  /// Construct a new Run object.
  Run(const double rt_window, const double mz_window):
      donors({rt_window, mz_window}),
      acceptors({rt_window, mz_window})
  {
  }

  /// Insert a peak: a donor if it has an identification (see Util::feature_hit()), else an acceptor.
  void insert(const Feature& feature, const std::size_t map_index, std::optional<Util::FeatureHit> hit)
  {
    const bool donor = hit.has_value();
    FeatureRef ref(map_index, feature, std::move(hit));

    if (donor)
    {
      donors.insert(std::make_shared<Donor>(ref));
    }
    else
    {
      acceptors.insert(std::make_shared<Acceptor>(ref));
    }
  }

  /// Release all storage.
  void clear()
  {
    donors.clear();
    acceptors.clear();
  }

  // Allow direct access to the grid type.
  GridWithStorage<Donor> donors;
  GridWithStorage<Acceptor> acceptors;
};


} // namespace OpenMS::PipEcho
