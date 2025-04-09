// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;

START_TEST(GridBasedCluster, "$Id$")

GridBasedCluster::Point position(4.5, 5.5);
GridBasedCluster::Rectangle box(position,position);
std::vector<int> points;
points.push_back(1);
points.push_back(6);
points.push_back(2);
int propA = 1;
std::vector<int> propB;
propB.push_back(1);
propB.push_back(2);
propB.push_back(3);

GridBasedCluster* nullPointer = nullptr;
GridBasedCluster* ptr;

START_SECTION(GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices, const int &property_A, const std::vector<int> &properties_B))
    GridBasedCluster cluster(position, box, points, propA, propB);
    TEST_EQUAL(cluster.getCentre().getX(), 4.5);
    ptr = new GridBasedCluster(position, box, points, propA, propB);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices))
    GridBasedCluster cluster(position, box, points);
    TEST_EQUAL(cluster.getCentre().getX(), 4.5);
    ptr = new GridBasedCluster(position, box, points);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

GridBasedCluster cluster(position, box, points, propA, propB);

START_SECTION(Point getCentre() const)
    TEST_EQUAL(cluster.getCentre().getX(), 4.5);
    TEST_EQUAL(cluster.getCentre().getY(), 5.5);
END_SECTION

START_SECTION(Rectangle getBoundingBox() const)
    TEST_EQUAL(cluster.getBoundingBox().minX(), 4.5);
    TEST_EQUAL(cluster.getBoundingBox().maxY(), 5.5);
END_SECTION

START_SECTION(std::vector<int> getPoints() const)
    TEST_EQUAL(cluster.getPoints()[0], 1);
    TEST_EQUAL(cluster.getPoints()[2], 2);
END_SECTION

START_SECTION(int getPropertyA() const)
    TEST_EQUAL(cluster.getPropertyA(), 1);
END_SECTION

START_SECTION(std::vector<int> getPropertiesB() const)
    TEST_EQUAL(cluster.getPropertiesB()[0], 1);
    TEST_EQUAL(cluster.getPropertiesB()[2], 3);
END_SECTION

GridBasedCluster::Point position1(4.5, 5.5);
GridBasedCluster::Point position2(4.5, 6.5);
GridBasedCluster cluster1(position1, box, points, propA, propB);
GridBasedCluster cluster2(position2, box, points, propA, propB);

START_SECTION(bool operator<(GridBasedCluster other) const)
    TEST_EQUAL(cluster1 < cluster2, true);
END_SECTION

START_SECTION(bool operator>(GridBasedCluster other) const)
    TEST_EQUAL(cluster2 > cluster1, true);
END_SECTION

START_SECTION(bool operator==(GridBasedCluster other) const)
    TEST_TRUE(cluster1 == cluster1);
END_SECTION

END_TEST
