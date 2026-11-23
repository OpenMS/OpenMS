// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/ClusteringGrid.h>
#include <OpenMS/CONCEPT/HashUtils.h>

#include <functional>
#include <sstream>

namespace OpenMS
{
    
ClusteringGrid::ClusteringGrid(const std::vector<double> &grid_spacing_x, const std::vector<double> &grid_spacing_y)
:grid_spacing_x_(grid_spacing_x), grid_spacing_y_(grid_spacing_y), range_x_(grid_spacing_x.front(),grid_spacing_x.back()), range_y_(grid_spacing_y.front(),grid_spacing_y.back())
{
}

std::vector<double> ClusteringGrid::getGridSpacingX() const
{
    return grid_spacing_x_;
}

std::vector<double> ClusteringGrid::getGridSpacingY() const
{
    return grid_spacing_y_;
}

size_t ClusteringGrid::CellIndexHash::operator()(const CellIndex& cell_index) const noexcept
{
    size_t seed = hash_int(cell_index.first);
    hash_combine(seed, hash_int(cell_index.second));
    return seed;
}

void ClusteringGrid::addCluster(const CellIndex &cell_index, const int &cluster_index)
{
    // Add the new cluster index to the list of clusters of the hash grid cell (created if it does not yet exist).
    cells_[cell_index].push_back(cluster_index);
}

void ClusteringGrid::removeCluster(const CellIndex &cell_index, const int &cluster_index)
{
    auto it = cells_.find(cell_index);
    if (it != cells_.end())
    {
        it->second.remove(cluster_index);
        if (it->second.empty())
        {
            cells_.erase(it);
        }
    }
}

void ClusteringGrid::removeAllClusters()
{
    cells_.clear();
}

const std::list<int>& ClusteringGrid::getClusters(const CellIndex &cell_index) const
{
    static const std::list<int> no_clusters;
    auto it = cells_.find(cell_index);
    if (it == cells_.end())
    {
        return no_clusters; // empty/never-populated cell -> no clusters
    }
    return it->second;
}

ClusteringGrid::CellIndex ClusteringGrid::getIndex(const Point &position) const
{
    if (position.getX() < range_x_.first || position.getX() > range_x_.second || position.getY() < range_y_.first || position.getY() > range_y_.second)
    {
        std::stringstream stream;
        stream << "This position (x,y)=(" << position.getX() << "," << position.getY() << ") is outside the range of the grid. (" << range_x_.first << " <= x <= " << range_x_.second << ", " << range_y_.first << " <= y <= " << range_y_.second << ")";
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
    }
    
    int i = std::lower_bound(grid_spacing_x_.begin(), grid_spacing_x_.end(), position.getX(), std::less_equal< double >()) - grid_spacing_x_.begin();
    int j = std::lower_bound(grid_spacing_y_.begin(), grid_spacing_y_.end(), position.getY(), std::less_equal< double >()) - grid_spacing_y_.begin();
        
    return ClusteringGrid::CellIndex (i,j);
}

bool ClusteringGrid::isNonEmptyCell(const CellIndex &cell_index) const
{
    return cells_.contains(cell_index);
}

int ClusteringGrid::getCellCount() const
{
    return cells_.size();
}

}
