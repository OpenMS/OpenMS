// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HashClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HashClustering* ptr = 0;
START_SECTION(HashClustering())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(~HashClustering())
{
	delete ptr;
}
END_SECTION

START_SECTION((HashClustering(std::vector< DataPoint > &data, DoubleReal rt_threshold, DoubleReal mz_threshold, ClusteringMethod &method_)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void performClustering()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void getSubtrees(std::vector< std::vector< SILACTreeNode > > &subtrees)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void createClusters(std::vector< std::vector< DataPoint * > > &clusters)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((std::vector<std::vector<Real> > getSilhouetteValues()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([HashClustering::InsufficientInput] InsufficientInput(const char *file, int line, const char *function, const char *message="not enough data points to cluster anything")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([HashClustering::InsufficientInput] virtual ~InsufficientInput()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

