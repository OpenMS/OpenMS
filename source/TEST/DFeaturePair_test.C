// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DFeaturePair_test.C,v 1.1 2006/04/10 16:11:13 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/KERNEL/DFeature.h>

///////////////////////////

START_TEST(DFeaturePair<D>, "$Id: DFeaturePair_test.C,v 1.1 2006/04/10 16:11:13 ole_st Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DFeaturePair<10>* d10_ptr = 0;
CHECK(DFeaturePair())
	d10_ptr = new DFeaturePair<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~DFeaturePair())
	delete d10_ptr;
RESULT

CHECK(==)

	DFeaturePair<2> p1;
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	DFeature<2> f2;
	f2.getPosition()[0] = 3.0;
	f2.getPosition()[1] = 4.0;
	
	p1.setFirst(f1);
	p1.setSecond(f2);
	p1.setQuality(5.0);
		
	DFeaturePair<2> p2;
	DFeature<2> f3;
	f3.getPosition()[0] = 1.0;
	f3.getPosition()[1] = 2.0;
	DFeature<2> f4;
	f4.getPosition()[0] = 3.0;
	f4.getPosition()[1] = 4.0;
	
	p2.setFirst(f3);
	p2.setSecond(f4);
	p2.setQuality(5.0);
	
	TEST_EQUAL(p1==p2,true);

RESULT

CHECK(!=)

	DFeaturePair<2> p1;
	DFeature<2> f1;
	f1.getPosition()[0] = 2.0;
	f1.getPosition()[1] = 2.0;
	DFeature<2> f2;
	f2.getPosition()[0] = 2.0;
	f2.getPosition()[1] = 2.0;
	
	p1.setFirst(f1);
	p1.setSecond(f2);
	
	DFeaturePair<2> p2;
	DFeature<2> f3;
	f3.getPosition()[0] = 1.0;
	f3.getPosition()[1] = 1.0;
	DFeature<2> f4;
	f4.getPosition()[0] = 1.0;
	f4.getPosition()[1] = 1.0;
	
	p2.setFirst(f3);
	p2.setSecond(f4);
	
	TEST_EQUAL(p1!=p2,true);

RESULT

CHECK(const QualityType& getOverallQuality() const)
	const DFeature<10> p;
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT

CHECK(QualityType& getQuality())
	DFeaturePair<3> p;
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
	p.getQuality() = 123.456;
	TEST_REAL_EQUAL(p.getQuality(), 123.456)
	p.getQuality() = -0.12345;
	TEST_REAL_EQUAL(p.getQuality(), -0.12345)
	p.getQuality() = 0.0;
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
RESULT

CHECK(setQuality(QualityType))
	DFeaturePair<3> p;
	p.setQuality(123.456);
	TEST_REAL_EQUAL(p.getQuality(), 123.456)
	p.setQuality(-0.12345);
	TEST_REAL_EQUAL(p.getQuality(), -0.12345)
	p.setQuality(0.0);
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
RESULT


CHECK(setFirst())
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setFirst(f1);

	DFeature<2> f2;
	f2 = p.getFirst();

	TEST_EQUAL(f1,f2)

RESULT

CHECK(setSecond())
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setSecond(f1);

	DFeature<2> f2;
	f2 = p.getSecond();

	TEST_EQUAL(f1,f2)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
