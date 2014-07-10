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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/COMPARISON/CLUSTERING/HashGrid2.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_LOCALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_LOCALCLUSTERING_H

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
class OPENMS_DLLAPI LocalClustering
{
    public:
    /**
    * @brief cluster centre, cluster bounding box, hash grid index
    */
    typedef Cluster::Point Point;    // DPosition<2>
    typedef Cluster::Rectangle Rectangle;    // DBoundingBox<2>
    typedef HashGrid2::CellIndex CellIndex;    // std::pair<int,int>

    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     * @param scaling_y    scaling in y-direction
     */
    LocalClustering(const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y, double scaling_y);

    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param grid_spacing_x    grid spacing in x-direction
     * @param grid_spacing_y    grid spacing in y-direction
     * @param scaling_y    scaling in y-direction
     */
    LocalClustering(const std::vector<double> &data_x, const std::vector<double> &data_y, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y, double scaling_y);

    /**
     * @brief performs the hierarchical clustering
     * (merges clusters until their dimension exceeds that of  cell)
     */
    void cluster();
    
    /**
     * @brief extends clusters in y-direction if possible
     * (merges clusters further in y-direction, i.e. clusters can now span multiple cells)
     */
    void extendClustersY();
    
    /**
     * @brief removes clusters with bounding boxe dimension in y-direction below certain threshold
     * @param threshold_y    minimal dimension of the cluster bounding box
     */
    void removeSmallClustersY(double threshold_y);

    /**
     * @brief returns final results (mapping of cluster indices to clusters)
     */
    std::map<int, Cluster> getResults() const;

    /**
    * @brief basic data structure for distances between clusters
    */
    class OPENMS_DLLAPI MinimumDistance
    {
        public:
        /**
        * @brief constructor
        */
        MinimumDistance(const int &cluster_index, const int &nearest_neighbour_index, const double &distance);
    
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
         * (for priority queue)
         */
        bool operator<(MinimumDistance other) const;
        bool operator>(MinimumDistance other) const;
        bool operator==(MinimumDistance other) const;
    
        private:        
        /// hide default Ctor
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
    * @brief hash grid on which the position of the clusters are registered
    * used in cluster method
    */
    HashGrid2 grid_;

    /**
    * @brief list of clusters
    * maps cluster indices to clusters
    */
    std::map<int, Cluster> clusters_;

    /**
    * @brief list of final clusters
    * i.e. clusters that are no longer merged
    */
    std::map<int, Cluster> clusters_final_;

    /**
    * @brief list of minimum distances
    * stores the smallest of the distances in the head
    */
    std::multiset<MinimumDistance, std::less<vector<MinimumDistance>::value_type> > distances_;
    
    /**
    * @brief scaling factor in y-direction
    * After scaling, typical clusters should be symmetric,
    * i.e. about the same dimension in x and y direction.
    */
    double scaling_y_;
 
    /**
     * @brief initialises all data structures
     * 
     * @param data_x    x-coordinates of points to be clustered
     * @param data_y    y-coordinates of points to be clustered
     * @param properties_A    property A of points (same in each cluster)
     * @param properties_B    property B of points (different in each cluster)
     */
    void init(const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B);

    /**
    * @brief set of indices referencing the points in the cluster
    */
    std::vector<int> minusOnes(int l);

    /**
    * @brief returns the distance between centres of two clusters
    * 
    * @param c1    centre of cluster 1
    * @param c2    centre of cluster 2
    * @param scaling    y-coordinates are multiplied by this factor
    * (might be necessary to ensure that typical clusters are about as high as they are wide.)
    * 
    * @return    distance between cluster centres
    */
    double clusterDistance(Point c1, Point c2, double scaling) const;
    
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
    bool mergeVeto(Cluster c1, Cluster c2) const;
    
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
    bool findNearestNeighbour(Cluster cluster, int cluster_index);

};

}

#endif /* OPENMS_COMPARISON_CLUSTERING_LOCALCLUSTERING_H */
