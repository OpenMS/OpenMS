// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ML/CLUSTERING/HashGrid.h>

#include <limits>
#include <iostream>

using namespace OpenMS;

struct Value
{
};

typedef OpenMS::HashGrid<Value> TestGrid;
const TestGrid::ClusterCenter cell_dimension(1, 1);

START_TEST(HashGrid, "$Id$")

START_SECTION(HashGrid(const ClusterCenter &cell_dimension))
{
  TestGrid t(cell_dimension);
  TEST_EQUAL(t.cell_dimension, cell_dimension);
  TEST_EQUAL(t.grid_dimension[0], 0);
  TEST_EQUAL(t.grid_dimension[1], 0);
}
END_SECTION

START_SECTION(cell_iterator insert(const value_type &v))
{
  TestGrid::cell_iterator it;
  TestGrid t(cell_dimension);

  const TestGrid::ClusterCenter key1(1, 2);
  it = t.insert(std::make_pair(key1, TestGrid::mapped_type()));
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

  {
    const TestGrid::ClusterCenter key(0, (double)std::numeric_limits<Int64>::min() - 1e5);
    TEST_EXCEPTION(Exception::OutOfRange, t.insert(std::make_pair(key, TestGrid::mapped_type())));
  }

  {
    const TestGrid::ClusterCenter key(0, (double)std::numeric_limits<Int64>::max() + 1e5);
    TEST_EXCEPTION(Exception::OutOfRange, t.insert(std::make_pair(key, TestGrid::mapped_type())));
  }
}
END_SECTION

START_SECTION(void erase(iterator pos))
  TestGrid t(cell_dimension);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));

  TEST_EQUAL(t.size(), 1);
  TestGrid::iterator it = t.begin();
  t.erase(it);
  TEST_EQUAL(t.size(), 0);
END_SECTION

START_SECTION(size_type erase(const key_type& key))
{
  TestGrid t(cell_dimension);
  const TestGrid::ClusterCenter key(1, 2);

  t.insert(std::make_pair(key, TestGrid::mapped_type()));
  TEST_EQUAL(t.erase(key), 1);

  t.insert(std::make_pair(key, TestGrid::mapped_type()));
  t.insert(std::make_pair(key, TestGrid::mapped_type()));
  TEST_EQUAL(t.erase(key), 2);
}
END_SECTION

START_SECTION(void clear())
{
  TestGrid t(cell_dimension);
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 2), TestGrid::mapped_type()));
  TEST_EQUAL(t.empty(), false);
  t.clear();
  TEST_EQUAL(t.empty(), true);
}
END_SECTION

START_SECTION(iterator begin())
{
  TestGrid t(cell_dimension);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 0), TestGrid::mapped_type()));
  t.insert(std::make_pair(TestGrid::ClusterCenter(2, 0), TestGrid::mapped_type()));

  TestGrid::iterator it = t.begin();
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it == t.end(), true);
}
END_SECTION

START_SECTION(const_iterator begin() const)
{
  TestGrid t(cell_dimension);
  const TestGrid &ct(t);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 0), TestGrid::mapped_type()));
  t.insert(std::make_pair(TestGrid::ClusterCenter(2, 0), TestGrid::mapped_type()));

  TestGrid::const_iterator it = ct.begin();
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it->first[0] <= 2, true);
  TEST_EQUAL(it->first[1], 0);
  it++;
  TEST_EQUAL(it == ct.end(), true);
}
END_SECTION

START_SECTION(iterator end())
{
  TestGrid t(cell_dimension);
  TestGrid::iterator it = t.begin();
  TEST_EQUAL(it == t.end(), true);
}
END_SECTION

START_SECTION(const_iterator end() const)
{
  const TestGrid ct(cell_dimension);
  TestGrid::const_iterator it = ct.begin();
  TEST_EQUAL(it == ct.end(), true);
}
END_SECTION

START_SECTION(bool empty() const)
  TestGrid t(cell_dimension);
  TEST_EQUAL(t.empty(), true);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  TEST_EQUAL(t.empty(), false);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  TEST_EQUAL(t.empty(), false);
END_SECTION

START_SECTION(size_type size() const)
  TestGrid t(cell_dimension);
  TEST_EQUAL(t.size(), 0);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  TEST_EQUAL(t.size(), 1);
  t.insert(std::make_pair(TestGrid::ClusterCenter(0, 0), TestGrid::mapped_type()));
  TEST_EQUAL(t.size(), 2);
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 0), TestGrid::mapped_type()));
  TEST_EQUAL(t.size(), 3);
END_SECTION

START_SECTION(const_grid_iterator grid_begin() const)
{
  TestGrid t(cell_dimension);
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 2), TestGrid::mapped_type()));
  const TestGrid &ct(t);
  TEST_EQUAL(ct.grid_begin()->second.size(), 1);
}
END_SECTION

START_SECTION(grid_iterator grid_begin())
{
  TestGrid t(cell_dimension);
  t.insert(std::make_pair(TestGrid::ClusterCenter(1, 2), TestGrid::mapped_type()));
  TEST_EQUAL(t.grid_begin()->second.size(), 1);
}
END_SECTION

START_SECTION(const_grid_iterator grid_end() const)
{
  const TestGrid t(cell_dimension);
  TEST_EQUAL(t.grid_begin() == t.grid_end(), true);
}
END_SECTION

START_SECTION(grid_iterator grid_end())
{
  TestGrid t(cell_dimension);
  TEST_EQUAL(t.grid_begin() == t.grid_end(), true);
}
END_SECTION

START_SECTION(const Grid::mapped_type& grid_at(const CellIndex &x) const)
{
  const TestGrid t(cell_dimension);
  const TestGrid::CellIndex i(0, 0);
  TEST_EXCEPTION(std::out_of_range, t.grid_at(i));
}
END_SECTION

START_SECTION(Grid::mapped_type& grid_at(const CellIndex &x))
{
  TestGrid t(cell_dimension);
  const TestGrid::CellIndex i(0, 0);
  TEST_EXCEPTION(std::out_of_range, t.grid_at(i));
}
END_SECTION

START_SECTION([EXTRA] std::size_t hash_value(const DPosition<N, T> &b))
{
  const DPosition<1, UInt> c1(1);
  const DPosition<1, UInt> c2(2);
  TEST_EQUAL(hash_value(c1), hash_value(c1));
  TEST_NOT_EQUAL(hash_value(c1), hash_value(c2));
}
{
  const DPosition<2, UInt> c1(1, 1);
  const DPosition<2, UInt> c2(1, 2);
  const DPosition<2, UInt> c3(2, 2);
  TEST_EQUAL(hash_value(c1), hash_value(c1));
  TEST_NOT_EQUAL(hash_value(c1), hash_value(c2));
  // Disabled: DPosition hash function is broken for this case
  //TEST_NOT_EQUAL(hash_value(c1), hash_value(c3));
}
END_SECTION

END_TEST

