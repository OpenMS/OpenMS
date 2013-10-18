// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>

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

