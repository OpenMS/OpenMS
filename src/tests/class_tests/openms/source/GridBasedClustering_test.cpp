// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ML/CLUSTERING/GridBasedClustering.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>

using namespace OpenMS;

START_TEST(GridBasedClustering, "$Id$")

MultiplexClustering::MultiplexDistance metric(1);

std::vector<double> grid_spacing_x;
std::vector<double> grid_spacing_y;
for (double i = 0; i <= 10; ++i)
{
    grid_spacing_x.push_back(i);
    grid_spacing_y.push_back(i);
}

std::vector<double> data_x;
std::vector<double> data_y;
std::vector<int> properties_A;
std::vector<int> properties_B;
for (int i = 0; i < 1000; ++i)
{
    data_x.push_back(5*(sin(static_cast<double>(i))+1));
    data_y.push_back(5*(sin(static_cast<double>(i+18))+1));
    properties_A.push_back(1);    // Should be the same within each cluster.
    properties_B.push_back(i);    // Should be different within each cluster.
}

GridBasedClustering<MultiplexClustering::MultiplexDistance>* nullPointer = nullptr;
GridBasedClustering<MultiplexClustering::MultiplexDistance>* ptr;

START_SECTION(GridBasedClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y))
    GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, properties_A, properties_B, grid_spacing_x, grid_spacing_y);
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
    ptr = new GridBasedClustering<MultiplexClustering::MultiplexDistance>(metric, data_x, data_y, properties_A, properties_B, grid_spacing_x, grid_spacing_y);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(GridBasedClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y))
    GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
    ptr = new GridBasedClustering<MultiplexClustering::MultiplexDistance>(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);

START_SECTION(void cluster())
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
END_SECTION

START_SECTION(std::map<int Cluster> getResults() const)
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
END_SECTION

START_SECTION(void extendClustersY())
    clustering.cluster();
    clustering.extendClustersY();
    TEST_EQUAL(clustering.getResults().size(), 11);
END_SECTION

START_SECTION(void removeSmallClustersY(double threshold_y))
    clustering.cluster();
    clustering.removeSmallClustersY(0.9);
    TEST_EQUAL(clustering.getResults().size(), 8);
END_SECTION

END_TEST
