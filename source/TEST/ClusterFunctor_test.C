// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterFunctor, "$Id$")

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

//interface class is not testable

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

START_SECTION((virtual void cluster(const DistanceMatrix< double > &original_distance, DistanceMatrix< double > &actual_distance, vector< vector< UInt > > &clusters, const String filepath="", const double threshold=1) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void registerChildren()))
{
  ClusterFunctor* cfp = Factory<ClusterFunctor>::create("AverageLinkage");
	TEST_EQUAL( cfp->getName() , "AverageLinkage")
  cfp = Factory<ClusterFunctor>::create("SingleLinkage");
	TEST_EQUAL( cfp->getName() , "SingleLinkage")
  cfp = Factory<ClusterFunctor>::create("CompleteLinkage");
	TEST_EQUAL( cfp->getName() , "CompleteLinkage")
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(ClusterFunctor::getProductName(),"ClusterFunctor")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



