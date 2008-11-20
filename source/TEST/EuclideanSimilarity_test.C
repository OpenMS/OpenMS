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
// $Maintainer: $
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

START_SECTION((double operator()(const std::pair< double, double > &a, const std::pair< double, double > &b) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2,2),make_pair(4,4)), 1-sqrt(8));
			TEST_REAL_SIMILAR(es(make_pair(9,0.1),make_pair(2.8,2)), 1-sqrt(42.05));
			TEST_REAL_SIMILAR(es(make_pair(12,0),make_pair(2,0)), 1-sqrt(100));
			es.setScale(sqrt(233.28));
}
END_SECTION

START_SECTION((double operator()(const std::pair< double, double > &c) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(9,0.1)), 1-0);
			TEST_REAL_SIMILAR(es(make_pair(2.8,2)), 1-0);
}
END_SECTION

START_SECTION((void setScale(UInt x)))
{
			EuclideanSimilarity es;
			es.setScale(10);
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2,2),make_pair(4,4)), 1-(sqrt(8)/10));
			TEST_REAL_SIMILAR(es(make_pair(9,0.1),make_pair(2.8,2)), 1-(sqrt(42.05)/10));
			TEST_REAL_SIMILAR(es(make_pair(12,0),make_pair(2,0)), 1-(sqrt(100)/10));
			es.setScale(sqrt(233.28));
			TEST_REAL_SIMILAR(es(make_pair(0.1,0.1),make_pair(10.9,10.9)), 1-1);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



