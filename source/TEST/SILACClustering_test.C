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
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>

#include <algorithm>
#include <set>

using namespace OpenMS;

typedef SILACClustering Test;
typedef Test::PointCoordinate Coordinate;
typedef std::map<Test::Grid::CellIndex, Test::Grid::CellContent> Cells;

const Coordinate cluster_dimension(1, 1);

Cells initClustering(DoubleReal rt_min, DoubleReal rt_max_spacing)
{
  Test t(cluster_dimension, rt_min, rt_max_spacing);

  t.insertPoint(Coordinate(0, 0), 0);
  t.insertPoint(Coordinate(0, .25), 0);
  t.insertPoint(Coordinate(0, .5), 0);
  t.insertPoint(Coordinate(.25, 0), 0);
  t.insertPoint(Coordinate(.25, .25), 0);
  t.insertPoint(Coordinate(.25, .5), 0);
  t.insertPoint(Coordinate(.5, 0), 0);
  t.insertPoint(Coordinate(.5, .25), 0);
  t.insertPoint(Coordinate(.5, .5), 0);

  t.insertPoint(Coordinate(1.5, 0), 0);
  t.insertPoint(Coordinate(1.5, .25), 0);
  t.insertPoint(Coordinate(1.5, .5), 0);
  t.insertPoint(Coordinate(1.75, 0), 0);
  t.insertPoint(Coordinate(1.75, .25), 0);
  t.insertPoint(Coordinate(1.75, .5), 0);

  t.cluster();

  Cells c;
  std::copy(t.grid.grid_begin(), t.grid.grid_end(), std::inserter(c, c.begin()));
  return c;
}

START_TEST(SILACClustering, "$Id$")

START_SECTION(SILACClustering(const PointCoordinate &cluster_dimension, DoubleReal rt_min, DoubleReal rt_max_spacing))
  Test t(cluster_dimension, 1, 2);
  TEST_EQUAL(t.grid.cell_dimension, cluster_dimension);
  TEST_EQUAL(t.rt_min, 1);
  TEST_EQUAL(t.rt_max_spacing, 2);
END_SECTION

START_SECTION(void cluster())
// Test main clustering
{
  Cells c(initClustering(0, 0));

  TEST_EQUAL(c.size(), 2);
  Cells::const_iterator it = c.begin();
  TEST_EQUAL(it->second.size(), 1);
  it++;
  TEST_EQUAL(it->second.size(), 1);
}
// Test rt_min
{
  Cells c(initClustering(.5, 0));

  // The 2 depends on the implementation
  TEST_EQUAL(c.size(), 2);
  Cells::const_iterator it = c.begin();
  TEST_EQUAL(it->second.size(), 1);
  it++;
  TEST_EQUAL(it->second.size(), 0);
}
// Test rt_max_spacing
{
  Cells c(initClustering(0, 1));

  // The 1 depends on the implementation
  TEST_EQUAL(c.size(), 1);
  Cells::const_iterator it = c.begin();
  TEST_EQUAL(it->second.size(), 1);
}
END_SECTION

END_TEST
