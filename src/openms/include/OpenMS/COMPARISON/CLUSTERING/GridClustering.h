// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>

//#include <OpenMS/KERNEL/StandardTypes.h>
//#include <OpenMS/DATASTRUCTURES/DRange.h>
//#include <OpenMS/COMPARISON/CLUSTERING/MultiplexGrid.h>
//#include <OpenMS/COMPARISON/CLUSTERING/MultiplexCluster.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_GRIDCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_GRIDCLUSTERING_H

using std::vector;

namespace OpenMS
{
/**
* @brief 2D hierarchical clustering implementation
* optimized for large datasets containing many small clusters
* i.e. dimensions of clusters << dimension of entire dataset
* 
* The clustering problem therefore simplifies to a number of
* local clustering problems. Each of the local problems can be
* solved on a couple of adjacent cells on a larger grid. The grid
* spacing is determined by the expected typical cluster size
* in this region.
* 
* Each datapoint can have two additional properties A and B.
* In each cluster all properties A need to be the same,
* all properties B different.
*/
template <typename Metric>
class OPENMS_DLLAPI GridClustering
{
    public:
    /**
    * @brief cluster centre, cluster bounding box, grid index
    */
    typedef MultiplexCluster::Point Point;    // DPosition<2>
    typedef MultiplexCluster::Rectangle Rectangle;    // DBoundingBox<2>
    typedef MultiplexGrid::CellIndex CellIndex;    // std::pair<int,int>
    
    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     */
    GridClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y)
    : metric_(metric), grid_(grid_spacing_x,grid_spacing_y)
    {
    }
    
    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     */
    GridClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y)
    : metric_(metric), grid_(grid_spacing_x,grid_spacing_y)
    {
    }
    
    /**
     * @brief basic data structure for distances between clusters
     */
    class OPENMS_DLLAPI MinimumDistance
    {
        public:
        /**
        * @brief constructor
        */
        MinimumDistance(const int &cluster_index, const int &nearest_neighbour_index, const double &distance)
        :cluster_index_(cluster_index), nearest_neighbour_index_(nearest_neighbour_index), distance_(distance)
        {
        }
    
        /**
        * @brief returns cluster index
        */
        int getClusterIndex() const
        {
            return cluster_index_;
        }
        
        /**
        * @brief returns index of nearest cluster
        */
        int getNearestNeighbourIndex() const
        {
            return  nearest_neighbour_index_;
        }
    
        /**
         * @brief operators for comparisons
         * (for priority queue)
         */
        bool operator<(MinimumDistance other) const
        {
            return distance_ < other.distance_;
        }
        bool operator>(MinimumDistance other) const
        {
            return distance_ > other.distance_;
        }
        bool operator==(MinimumDistance other) const
        {
            return distance_ == other.distance_;
        }

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
    
    private:
    /**
     * @brief metric for measuring the distance between points in the 2D plane
     */
    Metric metric_;
    
    /**
    * @brief grid on which the position of the clusters are registered
    * used in cluster method
    */
    MultiplexGrid grid_;

    /**
    * @brief list of clusters
    * maps cluster indices to clusters
    */
    std::map<int, MultiplexCluster> clusters_;

    /**
    * @brief list of final clusters
    * i.e. clusters that are no longer merged
    */
    std::map<int, MultiplexCluster> clusters_final_;

    /**
    * @brief list of minimum distances
    * stores the smallest of the distances in the head
    */
    //std::multiset<MinimumDistance, std::less<vector<MinimumDistance>::value_type> > distances_;
    
    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     */
    void init(const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B)
    {
        // fill the grid with points to be clustered (initially each cluster contains a single point)
        for (unsigned i = 0; i < data_x.size(); ++i)
        {
            Point position(data_x[i],data_y[i]);
            Rectangle box(position,position);
            
            std::vector<int> pi;    // point indicies
            pi.push_back(i);
            std::vector<int> pb;    // properties B
            pb.push_back(properties_B[i]);
            
            // add to cluster list
            MultiplexCluster cluster(position, box, pi, properties_A[i], pb);
            clusters_.insert(std::make_pair(i,cluster));
            
            // register on grid
            grid_.addCluster(grid_.getIndex(position), i);
        }
        
        // fill list of minimum distances
        std::map<int, MultiplexCluster>::iterator iterator = clusters_.begin();
        while (iterator != clusters_.end())
        {
            int cluster_index = iterator->first;
            MultiplexCluster cluster = iterator->second;
            
            if (findNearestNeighbour(cluster, cluster_index))
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
    bool mergeVeto(MultiplexCluster c1, MultiplexCluster c2) const
    {
        int A1 = c1.getPropertyA();
        int A2 = c2.getPropertyA();
        std::vector<int> B1 = c1.getPropertiesB();
        std::vector<int> B2 = c2.getPropertiesB();
        
        // check if any of the properties A and B is not set i.e. =-1
        if (A1 == -1 || A2 == -1 || std::find(B1.begin(), B1.end(), -1) != B1.end() || std::find(B2.begin(), B2.end(), -1) != B2.end())
        {
            return false;
        }
        
        // Will the merged cluster have the same properties A?
        bool vetoA = !(A1 == A2);
        
        // Will the merged cluster have different properties B?
        std::vector<int> B;
        sort(B1.begin(), B1.end());
        sort(B2.begin(), B2.end());
        set_intersection(B1.begin(),B1.end(),B2.begin(),B2.end(),back_inserter(B));
        bool vetoB = !B.empty();
        
        return (vetoA || vetoB);
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
     * @param Should the cluster be removed from the cluster list? 
     */
    bool findNearestNeighbour(MultiplexCluster cluster, int cluster_index)
    {
        Point centre = cluster.getCentre();
        CellIndex cell_index = grid_.getIndex(centre);
        
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
                    std::list<int> cluster_indices = grid_.getClusters(cell_index2);
                    for (std::list<int>::const_iterator cluster_index2 = cluster_indices.begin(); cluster_index2 != cluster_indices.end(); ++cluster_index2) {
                        if (*cluster_index2 != cluster_index)
                        {
                            MultiplexCluster cluster2 = clusters_.find(*cluster_index2)->second;
                            Point centre2 = cluster2.getCentre();
                            double distance = metric_(centre, centre2);
                            bool veto = mergeVeto(cluster, cluster2);    // If clusters cannot be merged anyhow, they are no nearest neighbours.
                            if (!veto && (distance < min_dist || nearest_neighbour == -1))
                            {
                                min_dist = distance;
                                nearest_neighbour = *cluster_index2;
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
        //distances_.insert(MinimumDistance(cluster_index, nearest_neighbour, min_dist));
        return false;
    }

};

}
#endif /* OPENMS_COMPARISON_CLUSTERING_GRIDCLUSTERING_H */
