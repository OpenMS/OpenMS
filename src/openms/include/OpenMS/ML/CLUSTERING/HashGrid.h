// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <cmath>
#include <iterator>
#include <limits>

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
   * @tparam Cluster Type to be stored in the hash grid. (e.g. HierarchicalClustering::Cluster)
   */
  template<typename Cluster>
  class HashGrid
  {
  public:
    /**
     * @brief Coordinate for stored pairs.
     */
    // XXX: Check is there is another type handy in OpenMS already
    typedef DPosition<2, double> ClusterCenter;

    /**
     * @brief Index for cells.
     */
    typedef DPosition<2, Int64> CellIndex;

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

  private:
    /**
     * @brief Element iterator for the hash grid.
     */
    class Iterator
    {
    private:
      friend class HashGrid;

      typedef typename Grid::iterator grid_iterator;
      typedef typename CellContent::iterator cell_iterator;

      /** @name Typedefs for STL compliance, these replace std::iterator
       */
      //@{
      /// The iterator's category type
      using iterator_category = std::input_iterator_tag;
      /// The iterator's value type
      using value_type = HashGrid::value_type;
      /// The reference type as returned by operator*()
      using reference = value_type&;
      /// The pointer type as returned by operator->()
      using pointer = value_type*;
      /// The difference type
      using difference_type = std::ptrdiff_t;
      //@}


      Grid& grid_;
      grid_iterator grid_it_;
      cell_iterator cell_it_;

      // Search for next non-empty cell
      void searchNextCell_()
      {
        while (cell_it_ == grid_it_->second.end())
        {
          grid_it_++;

          // If we are at the last cell, set cell iterator to something well-known
          if (grid_it_ == grid_.end())
          {
            cell_it_ = cell_iterator();
            return;
          }

          cell_it_ = grid_it_->second.begin();
        }
      }

    public:
      explicit Iterator(Grid& grid) : grid_(grid), grid_it_(grid.end())
      {
      }

      Iterator(Grid& grid, grid_iterator grid_it, cell_iterator cell_it) : grid_(grid), grid_it_(grid_it), cell_it_(cell_it)
      {
        searchNextCell_();
      }

      Iterator& operator++()
      {
        ++cell_it_;
        searchNextCell_();
        return *this;
      }

      const Iterator operator++(int)
      {
        Iterator ret(*this);
        operator++();
        return ret;
      }

      bool operator==(const Iterator& rhs) const
      {
        return grid_it_ == rhs.grid_it_ && cell_it_ == rhs.cell_it_;
      }

      bool operator!=(const Iterator& rhs) const
      {
        return !(*this == rhs);
      }

      value_type& operator*() const
      {
        return *cell_it_;
      }

      value_type* operator->() const
      {
        return &*cell_it_;
      }

      const CellIndex index() const
      {
        return grid_it_->first;
      }
    };

    /**
     * @brief Constant element iterator for the hash grid.
     */
    class ConstIterator
    {
    private:
      friend class HashGrid;

      typedef typename Grid::const_iterator grid_iterator;
      typedef typename CellContent::const_iterator cell_iterator;

      /** @name Typedefs for STL compliance, these replace std::iterator
       */
      //@{
      /// The iterator's category type
      using iterator_category = std::input_iterator_tag;
      /// The iterator's value type
      using value_type = HashGrid::value_type;
      /// The reference type as returned by operator*()
      using reference = value_type&;
      /// The pointer type as returned by operator->()
      using pointer = value_type*;
      /// The difference type
      using difference_type = std::ptrdiff_t;
      //@}

      const Grid& grid_;
      grid_iterator grid_it_;
      cell_iterator cell_it_;

      // Search for next non-empty cell
      void searchNextCell_()
      {
        while (cell_it_ == grid_it_->second.end())
        {
          grid_it_++;

          // If we are at the last cell, set cell iterator to something well-known
          if (grid_it_ == grid_.end())
          {
            cell_it_ = cell_iterator();
            return;
          }

          cell_it_ = grid_it_->second.begin();
        }
      }

    public:
      ConstIterator(const Grid& grid) : grid_(grid), grid_it_(grid.end())
      {
      }

      ConstIterator(const Grid& grid, grid_iterator grid_it, cell_iterator cell_it) : grid_(grid), grid_it_(grid_it), cell_it_(cell_it)
      {
        searchNextCell_();
      }

      ConstIterator(const Iterator& it) : grid_(it.grid_), grid_it_(it.grid_it_), cell_it_(it.cell_it_)
      {
      }

      ConstIterator& operator++()
      {
        ++cell_it_;
        searchNextCell_();
        return *this;
      }

      const ConstIterator operator++(int)
      {
        ConstIterator ret(*this);
        operator++();
        return ret;
      }

      bool operator==(const ConstIterator& rhs) const
      {
        return grid_it_ == rhs.grid_it_ && cell_it_ == rhs.cell_it_;
      }

      bool operator!=(const ConstIterator& rhs) const
      {
        return !(*this == rhs);
      }

      const value_type& operator*() const
      {
        return *cell_it_;
      }

      const value_type* operator->() const
      {
        return &*cell_it_;
      }

      const CellIndex index() const
      {
        return grid_it_->first;
      }
    };

  public:
    typedef ConstIterator const_iterator;
    typedef Iterator iterator;
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
    const CellIndex& grid_dimension;

  public:
    explicit HashGrid(const ClusterCenter& c_dimension) : cells_(), grid_dimension_(), cell_dimension(c_dimension), grid_dimension(grid_dimension_)
    {
    }

    /**
     * @brief Inserts a (2-dimensional coordinate, value) pair.
     * @param v Pair to be inserted.
     * @return Iterator that points to the inserted pair.
     */
    cell_iterator insert(const value_type& v)
    {
      const CellIndex cellkey = cellindexAtClustercenter_(v.first);
      CellContent& cell = cells_[cellkey];
      updateGridDimension_(cellkey);
      return cell.insert(v);
    }

    /**
     * @brief Erases element on given iterator.
     */
    void erase(iterator pos)
    {
      CellContent& cell = pos.grid_it_->second;
      cell.erase(pos.cell_it_);
    }

    /**
     * @brief Erases elements matching the 2-dimensional coordinate.
     * @param key Key of element to be erased.
     * @return Number of elements erased.
     */
    size_type erase(const key_type& key)
    {
      const CellIndex cellkey = cellindexAtClustercenter_(key);
      auto cell = cells_.find(cellkey);
      if (cell == cells_.end())
        return 0;
      return cell->second.erase(key);
    }

    /**
     * @brief Clears the map.
     */
    void clear()
    {
      cells_.clear();
    }

    /**
     * @brief Returns iterator to first element.
     */
    iterator begin()
    {
      grid_iterator grid_it = cells_.begin();
      if (grid_it == cells_.end())
        return end();

      cell_iterator cell_it = grid_it->second.begin();
      return iterator(cells_, grid_it, cell_it);
    }

    /**
     * @brief Returns iterator to first element.
     */
    const_iterator begin() const
    {
      const_grid_iterator grid_it = cells_.begin();
      if (grid_it == cells_.end())
        return end();

      const_cell_iterator cell_it = grid_it->second.begin();
      return const_iterator(cells_, grid_it, cell_it);
    }

    /**
     * @brief Returns iterator to first element.
     */
    iterator end()
    {
      return iterator(cells_);
    }

    /**
     * @brief Returns iterator to first element.
     */
    const_iterator end() const
    {
      return const_iterator(cells_);
    }

    /**
     * @brief Return true if HashGrid is empty.
     */
    bool empty() const
    {
      return size() == 0;
    }

    /**
     * @brief Return number of elements.
     */
    size_type size() const
    {
      size_type ret = 0;

      for (const_grid_iterator it = grid_begin(); it != grid_end(); ++it)
      {
        ret += it->second.size();
      }

      return ret;
    }

    /**
     * @brief Returns iterator to first grid cell.
     */
    const_grid_iterator grid_begin() const
    {
      return cells_.begin();
    }

    /**
     * @brief Returns iterator to on after last grid cell.
     */
    const_grid_iterator grid_end() const
    {
      return cells_.end();
    }

    /**
     * @brief Returns the grid cell at given index.
     */
    const typename Grid::mapped_type& grid_at(const CellIndex& x) const
    {
      return cells_.at(x);
    }

    /**
     * @warning Currently needed non-const by HierarchicalClustering.
     */
    typename Grid::mapped_type& grid_at(const CellIndex& x)
    {
      return cells_.at(x);
    }

    /**
     * @brief Returns the grid cell at given index if present, otherwise the grid_end iterator.
     */
    const_grid_iterator grid_find(const CellIndex& x) const
    {
      return cells_.find(x);
    }

    /**
     * @brief Returns the grid cell at given index if present, otherwise the grid_end iterator.
     */
    grid_iterator grid_find(const CellIndex& x)
    {
      return cells_.find(x);
    }


    /**
     * @warning Currently needed non-const by HierarchicalClustering.
     */
    grid_iterator grid_begin()
    {
      return cells_.begin();
    }
    /**
     * @warning Currently needed non-const by HierarchicalClustering.
     */
    grid_iterator grid_end()
    {
      return cells_.end();
    }

  private:
    // XXX: Replace with proper operator
    CellIndex cellindexAtClustercenter_(const ClusterCenter& key)
    {
      CellIndex ret;
      typename CellIndex::iterator it = ret.begin();
      typename ClusterCenter::const_iterator lit = key.begin(), rit = cell_dimension.begin();
      for (; it != ret.end(); ++it, ++lit, ++rit)
      {
        double t = std::floor(*lit / *rit);
        if (t < std::numeric_limits<Int64>::min() || t > std::numeric_limits<Int64>::max())
          throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        *it = static_cast<Int64>(t);
      }
      return ret;
    }

    void updateGridDimension_(const CellIndex& d)
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

  /** Hash value for OpenMS::DPosition. */
  template<UInt N, typename T>
  std::size_t hash_value(const DPosition<N, T>& b)
  {
    boost::hash<T> hasher;
    std::size_t hash = 0;
    for (typename DPosition<N, T>::const_iterator it = b.begin(); it != b.end(); ++it)
      hash ^= hasher(*it);
    return hash;
  }

} // namespace OpenMS

#endif /* OPENMS_COMPARISON_CLUSTERING_HASHGRID_H */
