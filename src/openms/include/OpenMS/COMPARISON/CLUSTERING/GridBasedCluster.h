// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_GRIDBASEDCLUSTER_H
#define OPENMS_COMPARISON_CLUSTERING_GRIDBASEDCLUSTER_H

namespace OpenMS
{
/**
* @brief basic data structure for clustering
*/
class OPENMS_DLLAPI GridBasedCluster
{
    public:
    /**
     * centre of a cluster
     */
    typedef DPosition<2> Point;

    /**
     * bounding box of a cluster
     */
    typedef DBoundingBox<2> Rectangle;

    /**
     * @brief initialises all data structures
     */
    GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices, const int &property_A, const std::vector<int> &properties_B);

    /**
     * @brief initialises all data structures
     */
    GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices);

    /**
     * @brief returns cluster centre
     */
    const Point& getCentre() const;

    /**
     * @brief returns bounding box
     */
    const Rectangle& getBoundingBox() const;

    /**
     * @brief returns indices of points in cluster
     */
    const std::vector<int>& getPoints() const;

    /**
     * @brief returns property A
     */
    int getPropertyA() const;

    /**
     * @brief returns properties B of all points
     */
    const std::vector<int>& getPropertiesB() const;
    
    /**
     * @brief operators for comparisons
     */
    bool operator<(const GridBasedCluster& other) const;
    bool operator>(const GridBasedCluster& other) const;
    bool operator==(const GridBasedCluster& other) const;
    
    private:
    /**
    * @brief centre of the cluster
    */
    Point centre_;

    /**
    * @brief bounding box of the cluster
    * i.e. (min,max) in x and y direction
    */
    Rectangle bounding_box_;

    /**
    * @brief set of indices referencing the points in the cluster
    */
    std::vector<int> point_indices_;

    /**
    * @brief properties A and B
    * Each point in a cluster can (optionally) possess two properties A and B.
    * For two points to be in the same cluster, they need to have the same
    * property A, e.g. the same charge. For two points to be in the same cluster,
    * they need to have different properties B, e.g. originate from two
    * different maps. -1 means properties are not set.
    */
    int property_A_;
    std::vector<int> properties_B_;
    
};

}

#endif /* OPENMS_COMPARISON_CLUSTERING_GRIDBASEDCLUSTER_H */
