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

#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <boost/unordered/unordered_set.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

namespace OpenMS
{
  /**
   * @brief Generic 2-dimensional hierarchical clustering with geometric hashing.
   *
   * The input data is saved into a hash grid. The dimension of the hash cells
   * is also the maximum cluster dimension.
   *
   * The clustering is performed on a 5x5 subsets of the hash grid. Only
   * clusters with all points in the inner 3x3 subset are accepted into the
   * output; all others are discarded. This makes sure that all clusters are
   * maximal and can't get larger with points not visible.
   *
   * This clustering only supports centroid linkage. It uses a priority queue
   * to save minimum distances between two subsets (proto-cluster?). No full
   * distance matrix is required.
   *
   * @tparam PointRef Reference associated with every point. Must have a default constructor.
   */
  template <typename PointRef>
  class HierarchicalClustering
  {
    public:
      /**
       * @brief Coordinate of a point to be clustered.
       *  @attention To be replaced by a %OpenMS coordinate type.
       */
      typedef DPosition<2, DoubleReal> PointCoordinate;

      /**
       *  @brief Bounding box of cluster.
       *  @attention To be replaced by OpenMS bounding box.
       */
      class BoundingBox
        : public std::pair<PointCoordinate, PointCoordinate>
      {
        public:
          BoundingBox(const PointCoordinate &p)
            : std::pair<PointCoordinate, PointCoordinate>(std::make_pair(p, p))
          { }

          BoundingBox(const BoundingBox &b)
            : std::pair<PointCoordinate, PointCoordinate>(b)
          { }

          PointCoordinate size() const
          {
            return this->second - this->first;
          }

          /** @brief Intersection of bounding box. */
          BoundingBox &operator|=(const BoundingBox &rhs)
          {
            typename PointCoordinate::iterator lit;
            typename PointCoordinate::const_iterator rit;

            // Calculate lower bound
            lit = this->first.begin(); rit = rhs.first.begin();
            for (; lit != this->first.end(); ++lit, ++rit) *lit = std::min(*lit, *rit);

            // Calculate upper bound
            lit = this->second.begin(); rit = rhs.second.begin();
            for (; lit != this->second.end(); ++lit, ++rit) *lit = std::max(*lit, *rit);

            return *this;
          }

          /** @brief Intersection of bounding box. */
          BoundingBox operator|(const BoundingBox &rhs) const
          {
            BoundingBox ret(*this);
            ret |= rhs;
            return ret;
          }

          operator PointCoordinate() const
          {
            // (first + second) / 2
            return coordScalarDiv_(this->first + this->second, 2);
          }
      };

      /**
       * @brief Set of points.
       * Describes a cluster on the grid. A point consists of a PointCoordinate and a PointRef.
       */
      class Cluster
        : public boost::unordered_multimap<PointCoordinate, PointRef>
      {
        public:
          BoundingBox bbox;

          Cluster(const BoundingBox &bbox)
            : bbox(bbox)
          { }
      };

      /**
       * @brief The hash grid data type.
       */
      typedef HashGrid<Cluster> Grid;

      /**
       * @brief The hash grid.
       *
       * It contains clusters.
       */
      Grid grid;

    protected:
      /** @brief Tree node used for clustering. */
      class TreeNode
      {
        public:
          const PointCoordinate coord;
          const BoundingBox bbox;
          TreeNode *left, *right;
          UInt points;
          const bool center;
          const PointRef ref;

          TreeNode(const PointCoordinate &coord, const PointRef &ref, bool center)
            : coord(coord), bbox(coord), left(0), right(0), points(1), center(center), ref(ref)
          { }

          TreeNode(const PointCoordinate &coord, const BoundingBox &bbox, TreeNode *left, TreeNode *right)
            : coord(coord), bbox(bbox),
              left(left), right(right),
              points(left->points + right->points),
              center(left->center && right->center),
              ref(PointRef())
          { }
      };

      typedef std::map<typename Grid::CellIndex, std::pair<typename Grid::CellContent *, bool> > ClusterCells;
      typedef boost::unordered_set<TreeNode *> ClusterTrees;

      /** @brief Wrapper class for two trees and the corresponding distance. */
      class TreeDistance
      {
        public:
          DoubleReal distance;
          TreeNode *left, *right;

          TreeDistance(const DoubleReal &distance, TreeNode *left, TreeNode *right)
            : distance(distance), left(left), right(right)
          { }

          bool operator>(const TreeDistance &rhs) const
          {
            return (distance > rhs.distance);
          }
      };

      /** @brief Priority queue queue used to find minimum distances. */
      typedef std::priority_queue<TreeDistance, std::vector<TreeDistance>, std::greater<TreeDistance> > TreeDistanceQueue;

    public:
      /**
       * @brief Constructor
       * @param cluster_dimension Max dimension of cluster
       */
      HierarchicalClustering(const PointCoordinate &cluster_dimension)
        : grid(cluster_dimension)
      { }

      /**
       * @brief Insert new PointCoordinate into grid.
       * @param d PointCoordinate to insert.
       * @param ref Associated caller specified info.
       * @return iterator to inserted cluster.
       */
      typename Grid::cell_iterator insertPoint(const PointCoordinate &d, const PointRef &ref)
      {
        typename Grid::cell_iterator it = insertCluster_(d);
        it->second.insert(std::make_pair(d, ref));
        return it;
      }

      /**
       * @brief Perform clustering of all existing points.
       */
      void cluster()
      {
        // Collect coordinates of all active cells
        std::vector<typename Grid::CellIndex> cells;
        for (typename Grid::const_grid_iterator it = grid.grid_begin(); it != grid.grid_end(); ++it)
          cells.push_back(it->first);
        // Cluster each available cell
        for (typename std::vector<typename Grid::CellIndex>::const_iterator it = cells.begin(); it != cells.end(); ++it)
          clusterIndex_(*it);
      }

    protected:
      /**
       * @brief Insert new Cluster into grid.
       * @param p Point to insert.
       * @return iterator to inserted cluster.
       */
      template <class P>
      typename Grid::cell_iterator insertCluster_(const P &p)
      {
        return grid.insert(std::make_pair(p, Cluster(p)));
      }

      /**
       * @brief Perform clustering at given cell index.
       * @param p Cell index.
       */
      void clusterIndex_(const typename Grid::CellIndex &p);

      /**
       * @brief Collect all cells used to cluster at given cell index.
       *
       * This function collects all cells in a 5x5 array.
       *
       * @param cur Cell index.
       * @param cells List of cells to be used.
       */
      void gridCells5x5_(typename Grid::CellIndex cur, ClusterCells &cells);

      /**
       * @brief Collect one cell.
       * @param cur Cell index.
       * @param cells List of cells.
       * @param center Is the given cell in the center.
       * @param ignore_missing Defines if non-existent errors should be ignored.
       */
      void gridCell_(const typename Grid::CellIndex &cur, ClusterCells &cells, bool center = false, bool ignore_missing = true)
      {
        try
        {
          cells.insert(std::make_pair(cur, std::make_pair(&grid.grid_at(cur), center)));
        }
        catch (std::out_of_range &)
        {
          if (!ignore_missing) throw;
        }
      }

      /**
       * @brief Add a new tree to the set of trees and distance queue
       */
      void addTreeDistance_(TreeNode *tree, ClusterTrees &trees, TreeDistanceQueue &dists)
      {
        // Infinity: no valid distance
        DoubleReal dist_min = std::numeric_limits<DoubleReal>::infinity();
        typename ClusterTrees::const_iterator dist_it = trees.end();

        // Generate minimal distance to existing trees
        for (typename ClusterTrees::const_iterator it = trees.begin(); it != trees.end(); ++it)
        {
          if (tree == *it) continue;
          DoubleReal dist = treeDistance_(tree, *it);
          if (dist < dist_min)
          {
            dist_min = dist;
            dist_it = it;
          }
        }

        // Insert distance if valid one found.
        if (dist_it != trees.end()) dists.push(TreeDistance(dist_min, tree, *dist_it));

        // Insert tree.
        trees.insert(tree);
      }

      /**
       * @brief Returns distance of two tree nodes
       * Returns the euclidean distance of the coordinates of the two trees.
       * It checks the size of the bounding box and returns INFINITY if it gets
       * to large.
       */
      DoubleReal treeDistance_(TreeNode *left, TreeNode *right)
      {
        const BoundingBox bbox = left->bbox | right->bbox;
        if (coordElemGreater_(bbox.size(), grid.cell_dimension))
        {
          return std::numeric_limits<DoubleReal>::infinity();
        }

        const PointCoordinate left_scaled = coordElemDiv_(left->coord, grid.cell_dimension);
        const PointCoordinate right_scaled = coordElemDiv_(right->coord, grid.cell_dimension);
        return coordDist_(left_scaled, right_scaled);
      }

      /**
       * @brief Recursively add the points of a finished cluster into the hash grid.
       * All points are saved in the leafs of the tree.
       * @param tree The tree
       * @param cluster The cluster
       */
      void tree2Cluster_(const TreeNode *tree, Cluster &cluster)
      {
        if (tree->left && tree->right)
        {
          tree2Cluster_(tree->left, cluster);
          tree2Cluster_(tree->right, cluster);
        }
        else
        {
          cluster.insert(std::make_pair(tree->bbox.first, tree->ref));
        }
        delete tree->left;
        delete tree->right;
      }

      /**
       * @brief Recursively add the points of an unfinished cluster back to the grid.
       * All points are saved in the leafs of the tree.
       * @param tree The tree
       */
      void tree2Points_(const TreeNode *tree)
      {
        if (tree->left && tree->right)
        {
          tree2Points_(tree->left);
          tree2Points_(tree->right);
        }
        else
        {
          insertPoint(tree->bbox.first, tree->ref);
        }
        delete tree->left;
        delete tree->right;
      }

      static PointCoordinate coordScalarDiv_(const PointCoordinate &lhs, const DoubleReal &rhs)
      {
        PointCoordinate ret;
        typename PointCoordinate::iterator it = ret.begin();
        typename PointCoordinate::const_iterator lit = lhs.begin();
        for (; it != ret.end(); ++it, ++lit) *it = *lit / rhs;
        return ret;
      }

      static PointCoordinate coordElemDiv_(const PointCoordinate &lhs, const PointCoordinate &rhs)
      {
        PointCoordinate ret;
        typename PointCoordinate::iterator it = ret.begin();
        typename PointCoordinate::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit / *rit;
        return ret;
      }

      static bool coordElemGreater_(const PointCoordinate &lhs, const PointCoordinate &rhs)
      {
        typename PointCoordinate::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; lit != lhs.end(); ++lit, ++rit) 
        {
          if (*lit > *rit) return true;
        }
        return false;
      }

      static DoubleReal coordDist_(const PointCoordinate &lhs, const PointCoordinate &rhs)
      {
        DoubleReal ret = 0;
        PointCoordinate p = lhs - rhs;
        typename PointCoordinate::const_iterator it = p.begin();
        for (; it != p.end(); ++it) ret += std::pow(*it, 2.);
        return std::sqrt(ret);
      }
  };

  template <typename I>
  void HierarchicalClustering<I>::clusterIndex_(const typename Grid::CellIndex &cur)
  {
    ClusterCells cells;
    ClusterTrees trees;
    TreeDistanceQueue dists;

    // Collect all cells we need
    try
    {
      gridCells5x5_(cur, cells);
    }
    catch (std::out_of_range &)
    { return; }

    // Collect and remove existing points from cells
    for (typename ClusterCells::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it)
    {
      typename Grid::CellContent &cell_cur = *cell_it->second.first;
      const bool &cell_center = cell_it->second.second;

      // Iterate per cluster
      typename Grid::cell_iterator cluster_tmp_it = cell_cur.begin();
      while (cluster_tmp_it != cell_cur.end())
      {
        typename Grid::cell_iterator cluster_it = cluster_tmp_it;
        ++cluster_tmp_it;

        // Check if it is not yet a cluster, aka have only one point
        if (cluster_it->second.size() == 1)
        {
          // Add each point to hash grid
          for (typename Cluster::const_iterator point_it = cluster_it->second.begin(); point_it != cluster_it->second.end(); ++point_it)
          {
            const PointCoordinate &coord = point_it->first;
            TreeNode *tree(new TreeNode(coord, point_it->second, cell_center));
            addTreeDistance_(tree, trees, dists);
          }

          // Remove point from hash grid cell
          cell_cur.erase(cluster_it);
        }
      }
    }

    // Try to join two subsets with minimum distance
    while (!dists.empty())
    {
      const typename TreeDistanceQueue::value_type cur_dist = dists.top();
      TreeNode *tree_left(cur_dist.left), *tree_right(cur_dist.right);
      dists.pop();
 
      // Check if both trees are not yet used with a smaller distance
      Size count_left = trees.count(tree_left), count_right = trees.count(tree_right);
      if (count_left && count_right)
      {
        trees.erase(tree_left);
        trees.erase(tree_right);

        const BoundingBox bbox = tree_left->bbox | tree_right->bbox;

        // Arithmethic mean: (left * left.points + right * right.points) / (left.points + right.points)
        const PointCoordinate &left = tree_left->coord, &right = tree_right->coord;
        const UInt &left_points = tree_left->points, &right_points = tree_right->points;
        const PointCoordinate coord = coordScalarDiv_(left * left_points + right * right_points, left_points + right_points);

        TreeNode *tree(new TreeNode(coord, bbox, tree_left, tree_right));

        addTreeDistance_(tree, trees, dists);
      }
      // Re-add a distance for the tree not yet used.
      // Otherwise this subset is lost even if it is not yet maximal.
      else if (count_left) addTreeDistance_(tree_left, trees, dists);
      else if (count_right) addTreeDistance_(tree_right, trees, dists);
    }

    // Add data back to grid
    for (typename ClusterTrees::iterator tree_it = trees.begin(); tree_it != trees.end(); ++tree_it)
    {
      // We got a finished tree with all points in the center, add cluster at centroid
      if ((**tree_it).center)
      {
        Cluster &cluster = insertCluster_((**tree_it).bbox)->second;
        tree2Cluster_(*tree_it, cluster);
      }
      // We got a finished tree but not all points in the center, readd as single points
      else
      {
        tree2Points_(*tree_it);
      }
      delete *tree_it;
    }
  }

  template <typename I>
  void HierarchicalClustering<I>::gridCells5x5_(typename Grid::CellIndex base, ClusterCells &cells)
  {
    // (0, 0)
    gridCell_(base, cells, true, false);

    typename Grid::CellIndex cur = base;
    cur[0] -= 2;
    // (-2, -2)
    cur[1] -= 2; gridCell_(cur, cells);
    // (-2, -1)
    cur[1] += 1; gridCell_(cur, cells);
    // (-2, 0)
    cur[1] += 1; gridCell_(cur, cells);
    // (-2, 1)
    cur[1] += 1; gridCell_(cur, cells);
    // (-2, 2)
    cur[1] += 1; gridCell_(cur, cells);

    cur = base; cur[0] -= 1;
    // (-1, -2)
    cur[1] -= 2; gridCell_(cur, cells);
    // (-1, -1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (-1, 0)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (-1, 1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (-1, 2)
    cur[1] += 1; gridCell_(cur, cells);

    cur = base;
    // (0, -2)
    cur[1] -= 2; gridCell_(cur, cells);
    // (0, -1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (0, 0)
    cur[1] += 1;
    // (0, 1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (0, 2)
    cur[1] += 1; gridCell_(cur, cells);

    cur = base; cur[0] += 1;
    // (1, -2)
    cur[1] -= 2; gridCell_(cur, cells);
    // (1, -1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (1, 0)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (1, 1)
    cur[1] += 1; gridCell_(cur, cells, true);
    // (1, 2)
    cur[1] += 1; gridCell_(cur, cells);

    cur = base; cur[0] += 2;
    // (2, -2)
    cur[1] -= 2; gridCell_(cur, cells);
    // (2, -1)
    cur[1] += 1; gridCell_(cur, cells);
    // (2, 0)
    cur[1] += 1; gridCell_(cur, cells);
    // (2, 1)
    cur[1] += 1; gridCell_(cur, cells);
    // (2, 2)
    cur[1] += 1; gridCell_(cur, cells);
  }
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
