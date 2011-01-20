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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/EuclideanSimilarity.h>
#include <cmath>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EuclideanSimilarity, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EuclideanSimilarity* ptr = 0;
START_SECTION(EuclideanSimilarity())
{
	ptr = new EuclideanSimilarity();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~EuclideanSimilarity())
{
	delete ptr;
}
END_SECTION

START_SECTION((EuclideanSimilarity(const EuclideanSimilarity &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((EuclideanSimilarity& operator=(const EuclideanSimilarity &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Real operator()(const std::pair< Real, Real > &a, const std::pair< Real, Real > &b) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2.0f,2.0f),make_pair(4.f,4.f)), 1-sqrt(8.0));
			TEST_REAL_SIMILAR(es(make_pair(9.f,0.1f),make_pair(2.8f,2.f)), 1-sqrt(42.05));
			TEST_REAL_SIMILAR(es(make_pair(12.f,0.0f),make_pair(2.f,0.0f)), 1-sqrt(100.0));
			es.setScale(sqrt(233.28f));
}
END_SECTION

START_SECTION((Real operator()(const std::pair< Real, Real > &c) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(9.0f,0.1f)), 1-0);
			TEST_REAL_SIMILAR(es(make_pair(2.8f,2.0f)), 1-0);
}
END_SECTION

START_SECTION((void setScale(Real x)))
{
			EuclideanSimilarity es;
			es.setScale(10);
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2.0f,2.0f),make_pair(4.f,4.f)), 1-(sqrt(8.0)/10));
			TEST_REAL_SIMILAR(es(make_pair(9.0f,0.1f),make_pair(2.8f,2.f)), 1-(sqrt(42.05)/10));
			TEST_REAL_SIMILAR(es(make_pair(12.0f,0.0f),make_pair(2.f,0.0f)), 1-(sqrt(100.0)/10));
			es.setScale(sqrt(233.28f));
			TEST_REAL_SIMILAR(es(make_pair(0.1f,0.1f),make_pair(10.9f,10.9f)), 1-1);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



