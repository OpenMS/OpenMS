// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusteringGrid.h>

#include <limits>
#include <iostream>

using namespace OpenMS;

START_TEST(ClusteringGrid, "$Id$")

std::vector<double> grid_spacing_x;
std::vector<double> grid_spacing_y;
for (double i = 0; i <=10; ++i)
{
    grid_spacing_x.push_back(i);
    grid_spacing_y.push_back(i);
}

ClusteringGrid* nullPointer = nullptr;
ClusteringGrid* ptr;

START_SECTION(ClusteringGrid(const std::vector<double> &grid_spacing_x, const std::vector<double> &grid_spacing_y))
    ClusteringGrid grid(grid_spacing_x, grid_spacing_y);
    TEST_EQUAL(grid.getGridSpacingX()[3], 3);
    ptr = new ClusteringGrid(grid_spacing_x, grid_spacing_y);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

ClusteringGrid grid(grid_spacing_x, grid_spacing_y);
ClusteringGrid::CellIndex index1 = std::make_pair(2,3);
ClusteringGrid::CellIndex index2 = std::make_pair(5,4);
ClusteringGrid::CellIndex index3 = std::make_pair(7,7);
ClusteringGrid::Point point(6.6,7.7);

START_SECTION(std::vector<double> getGridSpacingX() const)
    TEST_EQUAL(grid.getGridSpacingX()[3], 3);
    TEST_EQUAL(grid.getGridSpacingX()[10], 10);
END_SECTION

START_SECTION(std::vector<double> getGridSpacingY() const)
    TEST_EQUAL(grid.getGridSpacingY()[3], 3);
    TEST_EQUAL(grid.getGridSpacingY()[10], 10);
END_SECTION

START_SECTION(void addCluster(const CellIndex &cell_index, const int &cluster_index))
    grid.addCluster(index1,1);
    grid.addCluster(index2,2);
    TEST_EQUAL(grid.getCellCount(), 2);
END_SECTION

START_SECTION(void removeCluster(const CellIndex &cell_index, const int &cluster_index))
    grid.addCluster(index1,1);
    grid.addCluster(index2,2);
    grid.removeCluster(index2,2);
    TEST_EQUAL(grid.getCellCount(), 1);
END_SECTION

START_SECTION(void removeAllClusters())
    grid.addCluster(index1,1);
    grid.addCluster(index2,2);
    grid.removeAllClusters();
    TEST_EQUAL(grid.getCellCount(), 0);
END_SECTION

START_SECTION(std::list<int> getClusters(const CellIndex &cell_index) const)
    grid.addCluster(index1,1);
    grid.addCluster(index2,2);
    TEST_EQUAL(grid.getClusters(index1).front(), 1);
END_SECTION

START_SECTION(CellIndex getIndex(const Point &position) const)
    TEST_EQUAL(grid.getIndex(point).first, 7);
    TEST_EQUAL(grid.getIndex(point).second, 8);
END_SECTION

START_SECTION(bool isNonEmptyCell(const CellIndex &cell_index) const)
    grid.addCluster(index1,1);
    TEST_EQUAL(grid.isNonEmptyCell(index1), true);
    TEST_EQUAL(grid.isNonEmptyCell(index3), false);
END_SECTION

START_SECTION(int getCellCount() const)
    grid.addCluster(index1,1);
    grid.addCluster(index2,2);
    TEST_EQUAL(grid.getCellCount(), 2);
END_SECTION


END_TEST
