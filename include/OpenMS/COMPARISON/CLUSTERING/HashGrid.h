// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <boost/array.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/Hasher.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HASHGRID_H
#define OPENMS_COMPARISON_CLUSTERING_HASHGRID_H

namespace OpenMS
{
  /**
   * @brief Container for (2-dimensional coordinate, value) pairs.
   *
   * A hash-grid consists of hash-grid cells. The key of each cell is a pair of integers.
   * Each pair is assigned to a cell using a hash function.
   *
   * This container implements most parts of the C++ standard map interface.
   *
   * @tparam Cluster Type to be stored in the hash grid. (e.g. @ref{HierarchicalClustering::Cluster})
   */
  template <typename Cluster>
  class HashGrid
  {
    public:
      /**
       * @brief Coordinate for stored pairs.
       */
      // XXX: Check is there is another type handy in OpenMS already
      typedef DPosition<2, DoubleReal> ClusterCenter;

      /**
       * @brief Index for cells.
       */
      typedef DPosition<2, UInt> CellIndex;

      /**
       * @brief Contents of a cell.
       */
      typedef typename boost::unordered_multimap<ClusterCenter, Cluster> CellContent;

      /**
       * @brief Map of (cell-index, cell-content).
       */
      typedef boost::unordered_map<CellIndex, CellContent> Grid;

      typedef typename CellContent::key_type key_type;
      typedef typename CellContent::mapped_type mapped_type;
      typedef typename CellContent::value_type value_type;

      typedef typename Grid::const_iterator const_grid_iterator;
      typedef typename Grid::iterator grid_iterator;
      typedef typename CellContent::const_iterator const_cell_iterator;
      typedef typename CellContent::iterator cell_iterator;
      typedef typename CellContent::size_type size_type;

    private:
      Grid cells_;
      CellIndex grid_dimension_;

    public:
      /**
       * @brief Dimension of cells.
       */
      const ClusterCenter cell_dimension;
      /**
       * @brief Upper-right corner of key space for cells.
       */
      const CellIndex &grid_dimension;

    public:
      HashGrid(const ClusterCenter &cell_dimension)
        : cell_dimension(cell_dimension), grid_dimension(grid_dimension_)
      {
        // XXX: CellIndex needs constructor
        for (typename CellIndex::iterator it = grid_dimension_.begin(); it != grid_dimension_.end(); ++it) *it = 0;
      }

      /**
       * @brief Inserts a (2-dimensional coordinate, value) pair.
       * @param v Pair to be inserted.
       * @return Iterator that points to the inserted pair.
       */
      cell_iterator insert(const value_type &v)
      {
        const CellIndex cellkey = cellindex_at_clustercenter(v.first);
        CellContent &cell = cells_[cellkey];
        update_grid_dimension(cellkey);
        return cell.insert(v);
      }

      /**
       * @brief Erases elements matching the 2-dimensional coordinate.
       * @param x Key of element to be erased.
       * @return Number of elements erased.
       */
      size_type erase(const key_type &x)
      {   
        const CellIndex cellkey = cellindex_at_clustercenter(x);
        try
        {
          CellContent &cell = cells_.at(cellkey);
          return cell.erase(x);
        }
        catch (std::out_of_range &) { }
        return 0;
      }

      /**
       * @brief Clears the map.
       */
      void clear() { cells_.clear(); }

      /**
       * @brief Returns iterator to first grid cell.
       */
      const_grid_iterator grid_begin() const { return cells_.begin(); }

      /**
       * @brief Returns iterator to on after last grid cell.
       */
      const_grid_iterator grid_end() const { return cells_.end(); }

      /**
       * @brief Returns the grid cell at given index.
       */
      const typename Grid::mapped_type &grid_at(const CellIndex &x) const { return cells_.at(x); }

      /**
       * @warning Currently needed non-const by HierarchicalClustering.
       */
      grid_iterator grid_begin() { return cells_.begin(); }
      /**
       * @warning Currently needed non-const by HierarchicalClustering.
       */
      grid_iterator grid_end() { return cells_.end(); }
      /**
       * @warning Currently needed non-const by HierarchicalClustering.
       */
      typename Grid::mapped_type &grid_at(const CellIndex &x) { return cells_.at(x); }

    private:
      // XXX: Replace with proper operator
      CellIndex cellindex_at_clustercenter(const ClusterCenter &key)
      {
        CellIndex ret;
        typename CellIndex::iterator it = ret.begin();
        typename ClusterCenter::const_iterator lit = key.begin(), rit = cell_dimension.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit / *rit;
        return ret;
      }

      void update_grid_dimension(const CellIndex &d)
      {
        typename CellIndex::const_iterator it_new = d.begin();
        typename CellIndex::iterator it_cur = grid_dimension_.begin();
        for (; it_new != d.end(); ++it_new, ++it_cur)
        {
          if (*it_cur < *it_new)
            *it_cur = *it_new;
        }
      }
  };
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HASHGRID_H */
