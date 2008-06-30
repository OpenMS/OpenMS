// -*- Mode: C++; tab-width: 2; -*-
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

CHECK(ClusterFunctor())
{
  NOT_TESTABLE
}	
RESULT

CHECK(~ClusterFunctor())
{
  NOT_TESTABLE
}
RESULT

//interface class is not testable

CHECK((ClusterFunctor(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
RESULT

CHECK((ClusterFunctor& operator=(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual void cluster(const DistanceMatrix< double > &originalDist, DistanceMatrix< double > &actualDist, vector< vector< UInt > > &clusters, const String filepath="", const double threshold=1) const =0))
{
  NOT_TESTABLE
}
RESULT

CHECK((static void registerChildren()))
{
  ClusterFunctor* cfp = Factory<ClusterFunctor>::create("AverageLinkage");
	TEST_EQUAL( cfp->getName() , "AverageLinkage")
  cfp = Factory<ClusterFunctor>::create("SingleLinkage");
	TEST_EQUAL( cfp->getName() , "SingleLinkage")
  cfp = Factory<ClusterFunctor>::create("CompleteLinkage");
	TEST_EQUAL( cfp->getName() , "CompleteLinkage")
}
RESULT

CHECK((static const String getProductName()))
{
  TEST_EQUAL(ClusterFunctor::getProductName(),"ClusterFunctor")
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



