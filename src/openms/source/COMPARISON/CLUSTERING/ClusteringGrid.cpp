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
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusteringGrid.h>

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

void ClusteringGrid::addCluster(const CellIndex &cell_index, const int &cluster_index)
{
    if (cells_.find(cell_index) == cells_.end())
    {
        // If hash grid cell does not yet exist, create a new one.
        std::list<int> clusters;
        clusters.push_back(cluster_index);
        cells_.insert(std::make_pair(cell_index, clusters));
    }
    else
    {
        // If hash grid cell already exists, add the new cluster index to the existing list of clusters.
        cells_.find(cell_index)->second.push_back(cluster_index);
    }
}

void ClusteringGrid::removeCluster(const CellIndex &cell_index, const int &cluster_index)
{
    if (cells_.find(cell_index) != cells_.end())
    {
        cells_.find(cell_index)->second.remove(cluster_index);
        if (cells_.find(cell_index)->second.empty())
        {
            cells_.erase(cell_index);
        }
    }
}

void ClusteringGrid::removeAllClusters()
{
    cells_.clear();
}

std::list<int> ClusteringGrid::getClusters(const CellIndex &cell_index) const
{
    return cells_.find(cell_index)->second;
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
    return cells_.find(cell_index) != cells_.end();
}

int ClusteringGrid::getCellCount() const
{
    return cells_.size();
}

}
