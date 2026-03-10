// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "GridWithStorage.h"
#include "OpenMS/KERNEL/Feature.h"
#include "PeakTypes.h"

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

  /// Insert a donor peak.
  void insert(const Feature& feature, const std::size_t map_index)
  {
    Peak peak(map_index, feature);

    if (is_donor_feature(feature))
    {
      donors.insert(std::make_shared<Donor>(peak));
    }
    else
    {
      acceptors.insert(std::make_shared<Acceptor>(peak));
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

private:
  bool is_donor_feature(const Feature&);
};


} // namespace OpenMS::PipEcho
