// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterFunctor, "$Id: ClusterFunctor_test.C 5256 2009-05-13 06:25:08Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(ClusterFunctor())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(~ClusterFunctor())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ClusterFunctor(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ClusterFunctor& operator=(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void operator()(DistanceMatrix< Real > &original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold=1) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void registerChildren()))
{
  ClusterFunctor* cfp = Factory<ClusterFunctor>::create("AverageLinkage");
	TEST_NOT_EQUAL( dynamic_cast<AverageLinkage*>(cfp) , 0)
  cfp = Factory<ClusterFunctor>::create("SingleLinkage");
	TEST_NOT_EQUAL( dynamic_cast<SingleLinkage*>(cfp) , 0)
  cfp = Factory<ClusterFunctor>::create("CompleteLinkage");
	TEST_NOT_EQUAL( dynamic_cast<CompleteLinkage*>(cfp) , 0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



