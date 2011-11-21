// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>

using namespace OpenMS;

typedef OpenMS::HierarchicalClustering<UInt> Test;
typedef Test::PointCoordinate Coordinate;
const Coordinate cluster_dimension(1, 1);

START_TEST(HierarchicalClustering, "$Id$")

START_SECTION(HierarchicalClustering(const PointCoordinate &cluster_dimension))
  Test t(cluster_dimension);
  TEST_EQUAL(t.grid.cell_dimension, cluster_dimension);
END_SECTION

START_SECTION(Grid::cell_iterator insertPoint(const PointCoordinate &d, const PointRef &ref))
  Test t(cluster_dimension);
  const Coordinate coord(0, 0);
  t.insertPoint(coord, 0);
  t.insertPoint(coord, 1);
  TEST_EQUAL(t.grid.grid_begin()->second.size(), 2);
END_SECTION

START_SECTION(void cluster())
  Test t(cluster_dimension);
  const Coordinate coord(0, 0);
  t.insertPoint(coord, 0);
  t.insertPoint(coord, 1);
  t.cluster();
  TEST_EQUAL(t.grid.grid_begin()->second.size(), 1);
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] BoundingBox(const PointCoordinate &p))
  const Coordinate coord(1, 1);
  const Test::BoundingBox b(coord);
  TEST_EQUAL(b.first, coord);
  TEST_EQUAL(b.second, coord);
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] BoundingBox(const BoundingBox &b))
  const Coordinate coord(1, 1);
  const Test::BoundingBox b1(coord);
  const Test::BoundingBox b2(b1);
  TEST_EQUAL(b2.first, coord);
  TEST_EQUAL(b2.second, coord);
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] PointCoordinate size() const)
  const Coordinate coord1(1, 1);
  const Coordinate coord2(2, 2);
  Test::BoundingBox b(coord1);
  TEST_EQUAL(b.size(), Coordinate(0, 0));
  b |= coord2;
  TEST_EQUAL(b.size(), Coordinate(1, 1));
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] BoundingBox& operator|=(const BoundingBox &rhs))
  const Coordinate coord1(1, 1);
  const Coordinate coord2(2, 2);
  Test::BoundingBox b1(coord1);
  const Test::BoundingBox b2(coord2);
  b1 |= b2;
  TEST_EQUAL(b1.first, coord1);
  TEST_EQUAL(b1.second, coord2);
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] BoundingBox operator|(const BoundingBox &rhs) const)
  const Coordinate coord1(1, 1);
  const Coordinate coord2(2, 2);
  const Test::BoundingBox b1(coord1);
  const Test::BoundingBox b2(coord2);
  const Test::BoundingBox b3 = b1 | b2;
  TEST_EQUAL(b3.first, coord1);
  TEST_EQUAL(b3.second, coord2);
END_SECTION

START_SECTION([HierarchicalClustering::BoundingBox] operator PointCoordinate() const)
  const Coordinate coord1(1, 1);
  const Coordinate coord2(2, 2);
  Test::BoundingBox b(coord1);
  b |= coord2;
  const Coordinate c = b;
  TEST_EQUAL(c, Coordinate(1.5, 1.5));
END_SECTION

START_SECTION([HierarchicalClustering::Cluster] Cluster(const BoundingBox &bbox))
  const Coordinate coord(1, 1);
  const Test::BoundingBox b(coord);
  const Test::Cluster c(b);
  TEST_EQUAL(c.bbox.first, coord);
END_SECTION

END_TEST
