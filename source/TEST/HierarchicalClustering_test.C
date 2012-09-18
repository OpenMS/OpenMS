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
