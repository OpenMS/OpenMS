// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <vector>

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
