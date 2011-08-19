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

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>

using namespace OpenMS;

class Value
{
};

typedef OpenMS::HashGrid<Value> TestGrid;
const TestGrid::ClusterCenter cell_dimension(1, 1);

START_TEST(HashGrid, "$Id$")

START_SECTION(HashGrid(const ClusterCenter &cell_dimension))
{
  TestGrid t(cell_dimension);
  TEST_EQUAL(t.grid_dimension[0], 0);
  TEST_EQUAL(t.grid_dimension[1], 0);
}
END_SECTION

START_SECTION(cell_iterator insert(const value_type &v))
{
  TestGrid t(cell_dimension);
  const TestGrid::ClusterCenter key1(1, 2);
  TestGrid::cell_iterator it = t.insert(std::make_pair(key1, TestGrid::mapped_type()));
  TEST_EQUAL(t.grid_dimension[0], key1[0]);
  TEST_EQUAL(t.grid_dimension[1], key1[1]);
  TEST_EQUAL(it->first[0], key1[0]);
  TEST_EQUAL(it->first[1], key1[1]);
  const TestGrid::ClusterCenter key2(2, 3);
  it = t.insert(std::make_pair(key2, TestGrid::mapped_type()));
  TEST_EQUAL(t.grid_dimension[0], key2[0]);
  TEST_EQUAL(t.grid_dimension[1], key2[1]);
  TEST_EQUAL(it->first[0], key2[0]);
  TEST_EQUAL(it->first[1], key2[1]);
}
END_SECTION

START_SECTION(size_type erase(const key_type& key))
{
  TestGrid t(cell_dimension);
  const TestGrid::ClusterCenter key(1, 2);
  t.insert(std::make_pair(key, TestGrid::mapped_type()));
  TEST_EQUAL(t.erase(key), 1);
}
END_SECTION

START_SECTION(void clear())
{
  TestGrid t(cell_dimension);
  t.clear();
}
END_SECTION

START_SECTION(const_grid_iterator grid_begin() const)
{
  TestGrid t(cell_dimension);
  const TestGrid ct(t);
  const TestGrid::ClusterCenter key(1, 2);
  t.insert(std::make_pair(key, TestGrid::mapped_type()));
}
END_SECTION

START_SECTION(const_grid_iterator grid_end() const)
{
  TestGrid t(cell_dimension);
  const TestGrid ct(t);
  // Does not work, no output defined for iterators
  //TEST_EQUAL(ct.grid_begin(), ct.grid_end());
}
END_SECTION

START_SECTION(const typename Grid::mapped_type &grid_at(const CellIndex &x) const)
{
  TestGrid t(cell_dimension);
  const TestGrid ct(t);
  const TestGrid::CellIndex i(0, 0);
  TEST_EXCEPTION(std::out_of_range, ct.grid_at(i));
}
END_SECTION

END_TEST

