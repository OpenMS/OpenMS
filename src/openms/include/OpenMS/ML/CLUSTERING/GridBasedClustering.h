// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Johannes Veit $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ML/CLUSTERING/ClusteringGrid.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>

#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>

namespace OpenMS
{

  /**
   * @brief basic data structure for distances between clusters
   */
  class OPENMS_DLLAPI MinimumDistance
  {
public:

    /**
    * @brief constructor
    */
    MinimumDistance(const int& cluster_index, const int& nearest_neighbour_index, const double& distance);

    /**
    * @brief returns cluster index
    */
    int getClusterIndex() const;

    /**
    * @brief returns index of nearest cluster
    */
    int getNearestNeighbourIndex() const;

    /**
     * @brief operators for comparisons
     * (for multiset)
     */
    bool operator<(const MinimumDistance& other) const;
    bool operator>(const MinimumDistance& other) const;
    bool operator==(const MinimumDistance& other) const;

private:

    /// hide default constructor
    MinimumDistance();

    /**
    * @brief index in the cluster list
    */
    int cluster_index_;

    /**
    * @brief index of the nearest neighbour of the above cluster
    */
    int nearest_neighbour_index_;

    /**
    * @brief distance between cluster and its nearest neighbour
    */
    double distance_;

  };

  /**
  * @brief 2D hierarchical clustering implementation
  * optimized for large data sets containing many small clusters
  * i.e. dimensions of clusters << dimension of entire dataset
  *
  * The clustering problem therefore simplifies to a number of
  * local clustering problems. Each of the local problems can be
  * solved on a couple of adjacent cells on a larger grid. The grid
  * spacing is determined by the expected typical cluster size
  * in this region.
  *
  * Each data point can have two additional properties A and B.
  * In each cluster all properties A need to be the same,
  * all properties B different.
  */
  template <typename Metric>
  class GridBasedClustering :
    public ProgressLogger
  {
public:
    /**
    * @brief cluster centre, cluster bounding box, grid index
    */
    typedef GridBasedCluster::Point Point; // DPosition<2>
    typedef GridBasedCluster::Rectangle Rectangle; // DBoundingBox<2>
    typedef ClusteringGrid::CellIndex CellIndex; // std::pair<int,int>
    typedef std::multiset<MinimumDistance>::const_iterator MultisetIterator;
    typedef std::unordered_multimap<int, MultisetIterator>::const_iterator NNIterator;

    /**
     * @brief initialises all data structures
     * 
     * @param metric Metric for measuring the distance between points in the 2D plane
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     */
    GridBasedClustering(Metric metric, const std::vector<double>& data_x,
                        const std::vector<double>& data_y, const std::vector<int>& properties_A,
                        const std::vector<int>& properties_B, std::vector<double> grid_spacing_x,
                        std::vector<double> grid_spacing_y) :
      metric_(metric),
      grid_(grid_spacing_x, grid_spacing_y)
    {
      init_(data_x, data_y, properties_A, properties_B);
    }

    /**
     * @brief initialises all data structures
     *
     * @param metric Metric for measuring the distance between points in the 2D plane
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     */
    GridBasedClustering(Metric metric, const std::vector<double>& data_x,
                        const std::vector<double>& data_y, std::vector<double> grid_spacing_x,
                        std::vector<double> grid_spacing_y) :
      metric_(metric),
      grid_(grid_spacing_x, grid_spacing_y)
    {
      // set properties A and B to -1, i.e. ignore properties when clustering
      std::vector<int> properties_A(data_x.size(), -1);
      std::vector<int> properties_B(data_x.size(), -1);
      init_(data_x, data_y, properties_A, properties_B);
    }

    /**
     * @brief performs the hierarchical clustering
     * (merges clusters until their dimension exceeds that of  cell)
     */
    void cluster()
    {
      // progress logger
      // NOTE: for some reason, gcc7 chokes if we remove the OpenMS::String
      // below, so lets just not change it.
      Size clusters_start = clusters_.size();
      startProgress(0, clusters_start, OpenMS::String("clustering"));

      MinimumDistance zero_distance(-1, -1, 0);

      // combine clusters until all have been moved to the final list
      while (!clusters_.empty())
      {
        setProgress(clusters_start - clusters_.size());

        MultisetIterator smallest_distance_it = distances_.lower_bound(zero_distance);

        int cluster_index1 = smallest_distance_it->getClusterIndex();
        int cluster_index2 = smallest_distance_it->getNearestNeighbourIndex();

        eraseMinDistance_(smallest_distance_it);

        // update cluster list
        std::map<int, GridBasedCluster>::iterator cluster1_it = clusters_.find(cluster_index1);
        std::map<int, GridBasedCluster>::iterator cluster2_it = clusters_.find(cluster_index2);
        const GridBasedCluster& cluster1 = cluster1_it->second;
        const GridBasedCluster& cluster2 = cluster2_it->second;
        const std::vector<int>& points1 = cluster1.getPoints();
        const std::vector<int>& points2 = cluster2.getPoints();
        std::vector<int> new_points;
        new_points.reserve(points1.size() + points2.size());
        new_points.insert(new_points.end(), points1.begin(), points1.end());
        new_points.insert(new_points.end(), points2.begin(), points2.end());

        double new_x = (cluster1.getCentre().getX() * points1.size() + cluster2.getCentre().getX() * points2.size()) / (points1.size() + points2.size());
        double new_y = (cluster1.getCentre().getY() * points1.size() + cluster2.getCentre().getY() * points2.size()) / (points1.size() + points2.size());

        // update grid
        CellIndex cell_for_cluster1 = grid_.getIndex(cluster1.getCentre());
        CellIndex cell_for_cluster2 = grid_.getIndex(cluster2.getCentre());
        CellIndex cell_for_new_cluster = grid_.getIndex(DPosition<2>(new_x, new_y));
        grid_.removeCluster(cell_for_cluster1, cluster_index1);
        grid_.removeCluster(cell_for_cluster2, cluster_index2);
        grid_.addCluster(cell_for_new_cluster, cluster_index1);

        // merge clusters
        const Rectangle& box1 = cluster1.getBoundingBox();
        const Rectangle& box2 = cluster2.getBoundingBox();
        Rectangle new_box(box1);
        new_box.enlarge(box2.minPosition());
        new_box.enlarge(box2.maxPosition());

        // Properties A of both clusters should by now be the same. The merge veto has been checked
        // when a new entry to the minimum distance list was added, @see findNearestNeighbour_.
        if (cluster1.getPropertyA() != cluster2.getPropertyA())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Property A of both clusters not the same. ", "A");
        }
        int new_A = cluster1.getPropertyA();

        const std::vector<int>& B1 = cluster1.getPropertiesB();
        const std::vector<int>& B2 = cluster2.getPropertiesB();
        std::vector<int> new_B;
        new_B.reserve(B1.size() + B2.size());
        new_B.insert(new_B.end(), B1.begin(), B1.end());
        new_B.insert(new_B.end(), B2.begin(), B2.end());

        GridBasedCluster new_cluster(DPosition<2>(new_x, new_y), new_box, new_points, new_A, new_B);

        clusters_.erase(cluster1_it);
        clusters_.erase(cluster2_it);
        clusters_.insert(std::make_pair(cluster_index1, new_cluster));

        std::set<int> clusters_to_be_updated;
        clusters_to_be_updated.insert(cluster_index1);

        // erase distance object of cluster with cluster_index2 without updating (does not exist anymore!)
        // (the one with cluster_index1 has already been erased at the top of the while loop)
        eraseMinDistance_(distance_it_for_cluster_idx_[cluster_index2]);

        // find out which clusters need to be updated
        std::pair<NNIterator, NNIterator> nn_range = reverse_nns_.equal_range(cluster_index1);
        for (NNIterator nn_it = nn_range.first; nn_it != nn_range.second;)
        {
          clusters_to_be_updated.insert(nn_it->second->getClusterIndex());
          eraseMinDistance_((nn_it++)->second);
        }
        nn_range = reverse_nns_.equal_range(cluster_index2);
        for (NNIterator nn_it = nn_range.first; nn_it != nn_range.second;)
        {
          clusters_to_be_updated.insert(nn_it->second->getClusterIndex());
          eraseMinDistance_((nn_it++)->second);
        }

        // update clusters
        for (std::set<int>::const_iterator cluster_index = clusters_to_be_updated.begin(); cluster_index != clusters_to_be_updated.end(); ++cluster_index)
        {
          std::map<int, GridBasedCluster>::iterator c_it = clusters_.find(*cluster_index);
          const GridBasedCluster& c = c_it->second;
          if (findNearestNeighbour_(c, *cluster_index))
          {
            grid_.removeCluster(grid_.getIndex(c.getCentre()), *cluster_index);          // remove from grid
            clusters_.erase(c_it);          // remove from cluster list
          }
        }
      }

      endProgress();
    }

    /**
     * @brief extends clusters in y-direction if possible
     * (merges clusters further in y-direction, i.e. clusters can now span multiple cells)
     */
    void extendClustersY()
    {

      // construct new grid (grid only in x-direction, single cell in y-direction)
      std::vector<double> grid_spacing_x = grid_.getGridSpacingX();
      std::vector<double> grid_spacing_y = grid_.getGridSpacingY();
      std::vector<double> grid_spacing_y_new;
      grid_spacing_y_new.push_back(grid_spacing_y.front());
      grid_spacing_y_new.push_back(grid_spacing_y.back());
      ClusteringGrid grid_x_only(grid_spacing_x, grid_spacing_y_new);

      // register final clusters on the new grid
      for (std::map<int, GridBasedCluster>::const_iterator it = clusters_final_.begin(); it != clusters_final_.end(); ++it)
      {
        int cluster_index = it->first;
        GridBasedCluster cluster = it->second;
        grid_x_only.addCluster(grid_x_only.getIndex(cluster.getCentre()), cluster_index);
      }


      // scan through x on the grid
      for (unsigned cell = 0; cell < grid_spacing_x.size(); ++cell)
      {
        CellIndex grid_index(cell, 1);
        if (grid_x_only.isNonEmptyCell(grid_index))
        {
          std::list<int> cluster_indices = grid_x_only.getClusters(grid_index);          // indices of clusters in this x-range
          if (cluster_indices.size() > 1)
          {
            // First, order the clusters in ascending y.
            std::list<GridBasedCluster> cluster_list;            // list to order clusters in y
            std::map<GridBasedCluster, int> index_list;           // allows us to keep track of cluster indices after sorting
            for (std::list<int>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
            {
              cluster_list.push_back(clusters_final_.find(*it)->second);
              index_list.insert(std::make_pair(clusters_final_.find(*it)->second, *it));
            }
            cluster_list.sort();

            // Now check if two adjacent clusters c1 and c2 can be merged.
            std::list<GridBasedCluster>::const_iterator c1 = cluster_list.begin();
            std::list<GridBasedCluster>::const_iterator c2 = cluster_list.begin();
            ++c1;
            while (c1 != cluster_list.end())
            {
              double centre1x = (*c1).getCentre().getX();
              double centre1y = (*c1).getCentre().getY();
              double centre2x = (*c2).getCentre().getX();
              double centre2y = (*c2).getCentre().getY();

              double box1x_min = (*c1).getBoundingBox().minX();
              double box1x_max = (*c1).getBoundingBox().maxX();
              double box1y_min = (*c1).getBoundingBox().minY();
              double box1y_max = (*c1).getBoundingBox().maxY();
              double box2x_min = (*c2).getBoundingBox().minX();
              double box2x_max = (*c2).getBoundingBox().maxX();
              double box2y_min = (*c2).getBoundingBox().minY();
              double box2y_max = (*c2).getBoundingBox().maxY();

              //double y_range1 = box1y_max - box1y_min;
              //double y_range2 = box2y_max - box2y_min;
              //double y_gap = box1y_min - box2y_max;

              // Is there an overlap of the two clusters in x?
              bool overlap = (box1x_min <= box2x_max && box1x_min >= box2x_min) || (box1x_max >= box2x_min && box1x_max <= box2x_max);

              // Is the x-centre of one cluster in the x-range of the other?
              //bool centre_in_range1 = (box2x_min <= centre1x && centre1x <= box2x_max);
              //bool centre_in_range2 = (box1x_min <= centre2x && centre2x <= box1x_max);

              // Is the y-gap between the two clusters smaller than 1/s of their average y-range?
              //double s = 6;    // scaling factor
              //bool clusters_close = (y_gap * s <= (y_range1 - y_range2)/2);

              // Shall we merge the two adjacent clusters?
              //if ((centre_in_range1 || centre_in_range2) && clusters_close)
              if (overlap)
              {
                std::vector<int> points1 = (*c1).getPoints();
                std::vector<int> points2 = (*c2).getPoints();
                std::vector<int> new_points;
                new_points.reserve(points1.size() + points2.size());
                new_points.insert(new_points.end(), points1.begin(), points1.end());
                new_points.insert(new_points.end(), points2.begin(), points2.end());

                double new_x = (centre1x * points1.size() + centre2x * points2.size()) / (points1.size() + points2.size());
                double new_y = (centre1y * points1.size() + centre2y * points2.size()) / (points1.size() + points2.size());

                double min_x = std::min(box1x_min, box2x_min);
                double min_y = std::min(box1y_min, box2y_min);
                double max_x = std::max(box1x_max, box2x_max);
                double max_y = std::max(box1y_max, box2y_max);

                Point new_centre(new_x, new_y);
                Point position_min(min_x, min_y);
                Point position_max(max_x, max_y);
                Rectangle new_bounding_box(position_min, position_max);

                GridBasedCluster new_cluster(new_centre, new_bounding_box, new_points);

                // update final cluster list
                clusters_final_.erase(clusters_final_.find(index_list.find(*c1)->second));
                clusters_final_.erase(clusters_final_.find(index_list.find(*c2)->second));
                clusters_final_.insert(std::make_pair(index_list.find(*c1)->second, new_cluster));

                // update grid
                CellIndex cell_for_cluster1 = grid_x_only.getIndex((*c1).getCentre());
                CellIndex cell_for_cluster2 = grid_x_only.getIndex((*c2).getCentre());
                CellIndex cell_for_new_cluster = grid_x_only.getIndex(new_centre);

                grid_x_only.removeCluster(cell_for_cluster1, index_list.find(*c1)->second);
                grid_x_only.removeCluster(cell_for_cluster2, index_list.find(*c2)->second);
                grid_x_only.addCluster(cell_for_new_cluster, index_list.find(*c1)->second);
              }
              ++c1;
              ++c2;
            }
          }
        }
      }

    }

    /**
     * @brief removes clusters with bounding box dimension in y-direction below certain threshold
     * @param threshold_y    minimal dimension of the cluster bounding box
     */
    void removeSmallClustersY(double threshold_y)
    {
      std::map<int, GridBasedCluster>::iterator it = clusters_final_.begin();
      while (it != clusters_final_.end())
      {
        Rectangle box = it->second.getBoundingBox();
        if (box.maxY() - box.minY() < threshold_y)
        {
          clusters_final_.erase(it++);
        }
        else
        {
          ++it;
        }
      }
    }

    /**
     * @brief returns final results (mapping of cluster indices to clusters)
     */
    std::map<int, GridBasedCluster> getResults() const
    {
      return clusters_final_;
    }

private:

    /**
     * @brief metric for measuring the distance between points in the 2D plane
     */
    Metric metric_;

    /**
    * @brief grid on which the position of the clusters are registered
    * used in cluster method
    */
    ClusteringGrid grid_;

    /**
    * @brief list of clusters
    * maps cluster indices to clusters
    */
    std::map<int, GridBasedCluster> clusters_;

    /**
    * @brief list of final clusters
    * i.e. clusters that are no longer merged
    */
    std::map<int, GridBasedCluster> clusters_final_;

    /**
    * @brief list of minimum distances
    * stores the smallest of the distances in the head
    */
    std::multiset<MinimumDistance> distances_;

    /**
     * @brief reverse nearest neighbor lookup table
     * for finding out which clusters need to be updated faster
     */
    std::unordered_multimap<int, std::multiset<MinimumDistance>::const_iterator> reverse_nns_;

    /**
     * @brief cluster index to distance iterator lookup table
     * for finding out which clusters need to be updated faster
     */
    std::unordered_map<int, std::multiset<MinimumDistance>::const_iterator> distance_it_for_cluster_idx_;

    /**
     * @brief initialises all data structures
     *
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     */
    void init_(const std::vector<double>& data_x, const std::vector<double>& data_y,
               const std::vector<int>& properties_A, const std::vector<int>& properties_B)
    {
      // fill the grid with points to be clustered (initially each cluster contains a single point)
      for (unsigned i = 0; i < data_x.size(); ++i)
      {
        Point position(data_x[i], data_y[i]);
        Rectangle box(position, position);

        std::vector<int> pi;        // point indices
        pi.push_back(i);
        std::vector<int> pb;        // properties B
        pb.push_back(properties_B[i]);

        // add to cluster list
        GridBasedCluster cluster(position, box, pi, properties_A[i], pb);
        clusters_.insert(std::make_pair(i, cluster));

        // register on grid
        grid_.addCluster(grid_.getIndex(position), i);
      }

      // fill list of minimum distances
      std::map<int, GridBasedCluster>::iterator iterator = clusters_.begin();
      while (iterator != clusters_.end())
      {
        int cluster_index = iterator->first;
        const GridBasedCluster& cluster = iterator->second;

        if (findNearestNeighbour_(cluster, cluster_index))
        {
          // remove from grid
          grid_.removeCluster(grid_.getIndex(cluster.getCentre()), cluster_index);
          // remove from cluster list
          clusters_.erase(iterator++);
        }
        else
        {
          ++iterator;
        }
      }
    }

    /**
    * @brief checks if two clusters can be merged
    * Each point in a cluster can (optionally) have two properties A and B.
    * Properties A need to be the same, properties B need to differ in each cluster.
    * Methods checks if this is violated in the merged cluster.
    *
    * @param c1    cluster 1
    * @param c2    cluster 2
    *
    * @return    veto for merging clusters
    * true -> clusters can be merged
    * false -> clusters cannot be merged
    */
    bool mergeVeto_(const GridBasedCluster& c1, const GridBasedCluster& c2) const
    {
      int A1 = c1.getPropertyA();
      int A2 = c2.getPropertyA();

      // check if properties A of both clusters is set or not (not set := -1)
      if (A1 == -1 || A2 == -1)
      {
        return false;
      }

      // Will the merged cluster have the same properties A?
      if (A1 != A2) return true;

      std::vector<int> B1 = c1.getPropertiesB();
      std::vector<int> B2 = c2.getPropertiesB();

      // check if properties B of both clusters is set or not (not set := -1)
      if (std::find(B1.begin(), B1.end(), -1) != B1.end() || std::find(B2.begin(), B2.end(), -1) != B2.end())
      {
        return false;
      }

      // Will the merged cluster have different properties B?
      // (Hence the intersection of properties B of cluster 1 and cluster 2 should be empty.)
      std::vector<int> B_intersection;
      sort(B1.begin(), B1.end());
      sort(B2.begin(), B2.end());
      set_intersection(B1.begin(), B1.end(), B2.begin(), B2.end(), back_inserter(B_intersection));

      return !B_intersection.empty();
    }

    /**
     * @brief determines the nearest neighbour for each cluster
     *
     * @note If no nearest neighbour exists, the cluster is removed from the list.
     * The deletion is done outside of the method, see return boolean.
     * @note If two clusters cannot be merged (merge veto), they are no
     * viable nearest neighbours.
     *
     * @param cluster    cluster for which the nearest neighbour should be found
     * @param cluster_index    index of cluster
     *
     * @return Should the cluster be removed from the cluster list?
     */
    bool findNearestNeighbour_(const GridBasedCluster& cluster, int cluster_index)
    {
      const Point& centre = cluster.getCentre();
      const CellIndex& cell_index = grid_.getIndex(centre);
      double min_dist = 0;
      int nearest_neighbour = -1;

      // search in the grid cell and its 8 neighbouring cells for the nearest neighbouring cluster
      for (int i = -1; i <= 1; ++i)
      {
        for (int j = -1; j <= 1; ++j)
        {
          CellIndex cell_index2(cell_index);
          cell_index2.first += i;
          cell_index2.second += j;
          if (grid_.isNonEmptyCell(cell_index2))
          {
            const std::list<int>& cluster_indices = grid_.getClusters(cell_index2);
            for (std::list<int>::const_iterator cluster_index2 = cluster_indices.begin(); cluster_index2 != cluster_indices.end(); ++cluster_index2)
            {
              if (*cluster_index2 != cluster_index)
              {
                const GridBasedCluster& cluster2 = clusters_.find(*cluster_index2)->second;
                const Point& centre2 = cluster2.getCentre();
                double distance = metric_(centre, centre2);

                if (distance < min_dist || nearest_neighbour == -1)
                {
                  bool veto = mergeVeto_(cluster, cluster2); // If clusters cannot be merged anyhow, they are no nearest neighbours.
                  if (!veto)
                  {
                      min_dist = distance;
                      nearest_neighbour = *cluster_index2;
                  }
                }
              }
            }
          }
        }
      }

      if (nearest_neighbour == -1)
      {
        // no other cluster nearby, hence move the cluster to the final results
        clusters_final_.insert(std::make_pair(cluster_index, clusters_.find(cluster_index)->second));
        return true;
      }

      // add to the list of minimal distances
      std::multiset<MinimumDistance>::const_iterator it = distances_.insert(MinimumDistance(cluster_index, nearest_neighbour, min_dist));
      // add to reverse nearest neighbor lookup table
      reverse_nns_.insert(std::make_pair(nearest_neighbour, it));
      // add to cluster index -> distance lookup table
      distance_it_for_cluster_idx_[cluster_index] = it;

      return false;
    }

    /**
     * @brief remove minimum distance object and its related data
     *
     * Remove the distance object behind @p it from distances_ and remove all
     * corresponding data from the auxiliary data structures reverse_nns_ and
     * distance_it_for_cluster_idx_.
     *
     * @param it    Iterator of distance to be removed from distances_
     */
    void eraseMinDistance_(const std::multiset<MinimumDistance>::const_iterator it)
    {
      // remove corresponding entries from nearest neighbor lookup table
      std::pair<NNIterator, NNIterator> nn_range = reverse_nns_.equal_range(it->getNearestNeighbourIndex());
      for (NNIterator nn_it = nn_range.first; nn_it != nn_range.second; ++nn_it)
      {
        if (nn_it->second == it)
        {
          reverse_nns_.erase(nn_it);
          break;
        }
      }

      // remove corresponding entry from cluster index -> distance lookup table
      distance_it_for_cluster_idx_.erase(it->getClusterIndex());

      // remove from distances_
      distances_.erase(it);
    }
  };
}
