// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>

namespace OpenMS::PipEcho
{

/******************************************************************************/
/**
 * Window used when searching for acceptor features.
 */
struct Window
{
  double rt_tol; // Seconds.
  double mz_tol; // Daltons.
  std::size_t grid_neighbors;
};

} // namespace OpenMS::PipEcho
