// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
