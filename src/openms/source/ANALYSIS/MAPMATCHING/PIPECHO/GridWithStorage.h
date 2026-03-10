// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h"

#include <memory>

namespace OpenMS::PipEcho
{

/******************************************************************************/
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
template<typename T>
struct GridWithStorage
{
  // Handy aliases.
  using value_type = std::shared_ptr<T>;
  using grid_t = HashGrid<value_type>;
  using grid_center_t = grid_t::ClusterCenter;
  using grid_index_t = grid_t::CellIndex;

  // Access to the underlying storage.
  std::vector<value_type> storage;
  grid_t grid;

  /// Construct a new object given the grid center.
  GridWithStorage(const grid_center_t& center): grid(center) {};

  /// Insert a new object into the storage and grid.
  void insert(value_type);

  /// Explore cells near a starting point.
  template<typename Result>
  Result nearby(const grid_index_t&,
                std::size_t,
                Result,
                std::function<Result(Result, value_type)>) const;

  /// Remove all objects from the storage and grid.
  void clear()
  {
    grid.clear();
    storage.clear();
  }
};

/******************************************************************************/
template<typename T>
void GridWithStorage<T>::insert(value_type ptr)
{
  grid_center_t center(ptr->feature.getRT(), ptr->feature.getMZ());
  storage.push_back(ptr);
  grid.insert(std::make_pair(center, ptr));
}

/******************************************************************************/
template<typename T>
template<typename Result>
Result
GridWithStorage<T>::nearby(const grid_index_t& index,
                           std::size_t neighbors,
                           Result init,
                           std::function<Result(Result, value_type)> op) const
{
  // Check the cell where we expect to find the starting tval.
  if (auto cell = this->grid.grid_find(index); cell != this->grid.grid_end())
  {
    for (auto& tval : cell->second)
    {
      init = op(std::move(init), tval.second);
    }
  }

  if (neighbors < 1) { return init; }

  // Check nearby cells.
  const auto index_x = index.getX();
  const auto index_y = index.getY();

  for (auto i = index_x - neighbors; i <= index_x + neighbors; ++i)
  {
    for (auto j = index_y - neighbors; j <= index_y + neighbors; ++j)
    {
      const auto addr = grid_index_t(i, j);

      if (auto cell = this->grid.grid_find(addr); cell != this->grid.grid_end())
      {
        for (auto& tval : cell->second)
        {
          init = op(std::move(init), tval.second);
        }
      }
    }
  }

  return init;
}

} // namespace OpenMS::PipEcho
