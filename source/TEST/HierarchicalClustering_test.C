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
{
  Test t(cluster_dimension);
}
END_SECTION

START_SECTION(typename Grid::cell_iterator insertPoint(const PointCoordinate &d, const PointRef &ref))
{
  Test t(cluster_dimension);
  const Coordinate coord1(0, 0);
  t.insertPoint(coord1, 0);
  const Coordinate coord2(1, 1);
  t.insertPoint(coord2, 0);
}
END_SECTION

START_SECTION(void cluster())
{
  Test t(cluster_dimension);
  t.cluster();
}
END_SECTION

END_TEST
