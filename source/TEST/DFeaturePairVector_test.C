// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/KERNEL/DFeature.h>

#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DFeaturePairVector<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DFeaturePairVector<2>* pl_ptr = 0;
CHECK(DFeaturePairVector())
	pl_ptr = new DFeaturePairVector<2>();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK(~DFeaturePairVector())
	delete pl_ptr;
RESULT

CHECK( DFeaturePairVector& operator = (const DFeaturePairVector& rhs) )
	
	DFeaturePairVector<2> pvec;
		
	DFeature<2> f1;
	f1.getPosition()[0] = 2.0;
	f1.getPosition()[1] = 4.0;
	
	DFeature<2> f2;
	f2.getPosition()[0] = 3.0;
	f2.getPosition()[1] = 6.0;
	
	DFeaturePair<2> pair1;
	pair1.setFirst(f1);
	pair1.setSecond(f2);
	pvec.push_back(pair1);
	
	DFeature<2> f3;
	f3.getPosition()[0] = 4.0;
	f3.getPosition()[1] = 8.0;
	
	DFeature<2> f4;
	f4.getPosition()[0] = 5.0;
	f4.getPosition()[1] = 10.0;
	
	DFeaturePair<2> pair2;
	pair2.setFirst(f3);
	pair2.setSecond(f4);
	pvec.push_back(pair2);
	
	DFeaturePairVector<2> pvec_copy = pvec;
	
	TEST_EQUAL(pvec_copy.size(),2);
	
	DFeaturePairVector<2>::const_iterator cit = pvec_copy.begin();
	TEST_EQUAL(cit->getFirst().getPosition()[0],2.0);
	TEST_EQUAL(cit->getFirst().getPosition()[1],4.0);		
	
	cit++;
	TEST_EQUAL(cit->getSecond().getPosition()[0],5.0);
	TEST_EQUAL(cit->getSecond().getPosition()[1],10.0);
	
RESULT

CHECK(DFeaturePairVector(const DFeaturePairVector& vec))
	
	DFeaturePairVector<2> pvec;
		
	DFeature<2> f1;
	f1.getPosition()[0] = 2.0;
	f1.getPosition()[1] = 4.0;
	
	DFeature<2> f2;
	f2.getPosition()[0] = 3.0;
	f2.getPosition()[1] = 6.0;
	
	DFeaturePair<2> pair1;
	pair1.setFirst(f1);
	pair1.setSecond(f2);
	pvec.push_back(pair1);
	
	DFeature<2> f3;
	f3.getPosition()[0] = 4.0;
	f3.getPosition()[1] = 8.0;
	
	DFeature<2> f4;
	f4.getPosition()[0] = 5.0;
	f4.getPosition()[1] = 10.0;
	
	DFeaturePair<2> pair2;
	pair2.setFirst(f3);
	pair2.setSecond(f4);
	pvec.push_back(pair2);
	
	DFeaturePairVector<2> pvec_copy(pvec);
	
	TEST_EQUAL(pvec_copy.size(),2);
	
	DFeaturePairVector<2>::const_iterator cit = pvec_copy.begin();
	TEST_EQUAL(cit->getFirst().getPosition()[0],2.0);
	TEST_EQUAL(cit->getFirst().getPosition()[1],4.0);		
	
	cit++;
	TEST_EQUAL(cit->getSecond().getPosition()[0],5.0);
	TEST_EQUAL(cit->getSecond().getPosition()[1],10.0);
	
RESULT


CHECK( bool operator == (const DFeaturePairVector& rhs) const )
		
	DFeature<2> f1;
	f1.getPosition()[0] = 2.0;
	f1.getPosition()[1] = 4.0;
	
	DFeature<2> f2;
	f2.getPosition()[0] = 3.0;
	f2.getPosition()[1] = 6.0;
	
	DFeaturePair<2> pair1;
	pair1.setFirst(f1);
	pair1.setSecond(f2);
		
	DFeature<2> f3;
	f3.getPosition()[0] = 4.0;
	f3.getPosition()[1] = 8.0;
	
	DFeature<2> f4;
	f4.getPosition()[0] = 5.0;
	f4.getPosition()[1] = 10.0;
	
	DFeaturePair<2> pair2;
	pair2.setFirst(f3);
	pair2.setSecond(f4);
		
	DFeaturePairVector<2> pvec1;
	pvec1.push_back(pair1);
	pvec1.push_back(pair2);
	
	DFeaturePairVector<2> pvec2;
	pvec2.push_back(pair1);
	pvec2.push_back(pair2);
	
	TEST_EQUAL(pvec1 == pvec2, true)

RESULT

CHECK( bool operator != (const DFeaturePairVector& rhs) const )
	DFeature<2> f1;
	f1.getPosition()[0] = 2.0;
	f1.getPosition()[1] = 4.0;
	
	DFeature<2> f2;
	f2.getPosition()[0] = 3.0;
	f2.getPosition()[1] = 6.0;
	
	DFeaturePair<2> pair1;
	pair1.setFirst(f1);
	pair1.setSecond(f2);
		
	DFeature<2> f3;
	f3.getPosition()[0] = 4.0;
	f3.getPosition()[1] = 8.0;
	
	DFeature<2> f4;
	f4.getPosition()[0] = 5.0;
	f4.getPosition()[1] = 10.0;
	
	DFeaturePair<2> pair2;
	pair2.setFirst(f3);
	pair2.setSecond(f4);
		
	DFeaturePairVector<2> pvec1;
	pvec1.push_back(pair1);
	pvec1.push_back(pair2);
	
	DFeaturePairVector<2> pvec2;
	pair1.setQuality(1.0);
	pvec2.push_back(pair1);
	pvec2.push_back(pair2);
	
	TEST_EQUAL(pvec1 != pvec2, true)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
