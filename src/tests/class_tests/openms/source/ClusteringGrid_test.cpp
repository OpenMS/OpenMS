// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ML/CLUSTERING/ClusteringGrid.h>

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
