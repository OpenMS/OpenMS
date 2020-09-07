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
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <map>
#include <vector>
#include <list>

namespace OpenMS
{
/**
* @brief data structure to store 2D data to be clustered
* e.g. (m/z, retention time) coordinates from multiplex filtering
* 
* @see LocalClustering
*/
class OPENMS_DLLAPI ClusteringGrid
{
    public:
    /**
     * coordinates of a grid cell
     */
    typedef std::pair<int,int> CellIndex;
    
    /**
     * coordinates in x-y-plane
     */
    typedef DPosition<2> Point;
    
    /**
     * @brief constructor taking two vectors
     * @param grid_spacing_x    grid spacing in x direction
     * @param grid_spacing_y    grid spacing in y direction
     *
     * @note Vectors are assumed to be sorted.
     */
    ClusteringGrid(const std::vector<double> &grid_spacing_x, const std::vector<double> &grid_spacing_y);
    
    /**
    * @brief returns grid spacing in x direction
    */
    std::vector<double> getGridSpacingX() const;
    
    /**
    * @brief returns grid spacing in y direction
    */
    std::vector<double> getGridSpacingY() const;
    
    /**
    * @brief adds a cluster to this grid cell
    * 
    * @param cell_index    cell index (i,j) on the grid
    * @param cluster_index    index of the cluster in the cluster list
    */
    void addCluster(const CellIndex &cell_index, const int &cluster_index);
    
    /**
    * @brief removes a cluster from this grid cell
    * and removes the cell if no other cluster left
    * 
    * @param cell_index    cell index (i,j) on the grid
    * @param cluster_index    index of the cluster in the cluster list
    */
    void removeCluster(const CellIndex &cell_index, const int &cluster_index);

    /**
    * @brief removes all clusters from this grid (and hence all cells)
    */
    void removeAllClusters();

    /**
    * @brief returns clusters in this grid cell
    * 
    * @param cell_index    cell index (i,j) on the grid
    * @return list of cluster indices (from the list of clusters) which are centred in this cell
    */
    std::list<int> getClusters(const CellIndex &cell_index) const;

    /**
    * @brief returns grid cell index (i,j) for the positions (x,y)
    * 
    * @param position    coordinates (x,y) on the grid
    * @return cell index (i,j) of the cell in which (x,y) lies
    */
    CellIndex getIndex(const Point &position) const;

    /**
    * @brief checks if there are clusters at this cell index
    * 
    * @param cell_index    cell index (i,j) on the grid
    * @return true if there are clusters in this cell
    * 
    * @throw Exception::IllegalArgument if the coordinates (x,y) lie outside the grid.
    * @throw Exception::InvalidValue if one of the two indices is negative.
    */
    bool isNonEmptyCell(const CellIndex &cell_index) const;

    /**
    * @brief returns number of grid cells occupied by one or more clusters
    * 
    * @return number of non-empty cells
    */
    int getCellCount() const;

    private:
    /**
    * @brief spacing of the grid in x and y direction
    */
    const std::vector<double> grid_spacing_x_;
    const std::vector<double> grid_spacing_y_;

    /**
    * @brief [min, max] of the grid in x and y direction
    */
    std::pair <double,double> range_x_;  
    std::pair <double,double> range_y_;  

    /**
    * @brief grid cell index mapped to a list of clusters in it
    */
    std::map<CellIndex, std::list<int> > cells_;

};

}

