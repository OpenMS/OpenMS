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
     */
    GridClustering(Metric metric)
    : metric_(metric)
    {
    }
    
    private:
    /**
     * @brief metric for measuring the distance between points to be clustered 
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
    * @brief scaling factor in y-direction
    * After scaling, typical clusters should be symmetric,
    * i.e. about the same dimension in x and y direction.
    */
    double scaling_y_;
    
};

}
#endif /* OPENMS_COMPARISON_CLUSTERING_GRIDCLUSTERING_H */
