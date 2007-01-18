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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/KERNEL/DFeature.h>

///////////////////////////

START_TEST(DFeaturePair<D>, "$Id$")

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

CHECK( DFeaturePair(const DFeaturePair& fp) )
	DFeaturePair<2> p1;
	p1.setQuality(5.0);
	
	DFeaturePair<2> p2(p1);
	
	TEST_REAL_EQUAL( p1.getQuality(),p2.getQuality() )

RESULT

CHECK( DFeaturePair(FeatureType const & first, FeatureType const & second, QualityType const & quality = QualityType(0)) )
	DFeature<2> f1;
	DFeature<2> f2;
	
	DFeaturePair<2> pair(f1,f2);
	
	TEST_EQUAL( f1, pair.getFirst() )
	TEST_EQUAL( f2, pair.getSecond() )
RESULT

CHECK( DFeaturePair& operator = (const DFeaturePair& rhs) )
	DFeaturePair<2> p1;
	p1.setQuality( 5.0 );
	
	DFeaturePair<2> p2 = p1;
	
	TEST_REAL_EQUAL( p1.getQuality(),p2.getQuality() )
RESULT

CHECK(bool operator == (const DFeaturePair& rhs) const)

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

CHECK(bool operator != (const DFeaturePair& rhs) const)

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

CHECK(void setQuality(const QualityType& ql))
	DFeaturePair<3> p;
	p.setQuality(123.456);
	TEST_REAL_EQUAL(p.getQuality(), 123.456)
	p.setQuality(-0.12345);
	TEST_REAL_EQUAL(p.getQuality(), -0.12345)
	p.setQuality(0.0);
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
RESULT


CHECK(FeatureType& getFirst())
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setFirst(f1);

	DFeature<2> f2;
	f2 = p.getFirst();

	TEST_EQUAL(f1,f2)

RESULT

CHECK(FeatureType& getSecond())
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setSecond(f1);

	DFeature<2> f2;
	f2 = p.getSecond();

	TEST_EQUAL(f1,f2)

RESULT

CHECK(const FeatureType& getFirst() const)
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setFirst(f1);

	const DFeature<2> f2 = p.getFirst();
	TEST_EQUAL(f1,f2)

RESULT

CHECK(const FeatureType& getSecond() const)
	DFeaturePair<2> p;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 1.0;
	f1.getPosition()[1] = 2.0;
	p.setSecond(f1);

	const DFeature<2> f2 = p.getSecond();
	TEST_EQUAL(f1,f2)
	
RESULT

CHECK( const QualityType& getQuality() const )
	DFeaturePair<2> p;
	p.setQuality(3.0);
	const DFeaturePair<2>::QualityType q = p.getQuality();
	
	TEST_REAL_EQUAL( q,p.getQuality() )

RESULT

CHECK( void setFirst(const FeatureType& frt) )
	DFeaturePair<2> p;
	const DFeature<2> f;
	p.setFirst(f);
	
	TEST_EQUAL( f, p.getFirst() )
RESULT

CHECK( void setQuality(const QualityType& ql) )
	DFeaturePair<2> p;
	const DFeaturePair<2>::QualityType q = 10.0;
	p.setQuality(q);
	
	TEST_EQUAL( q, p.getQuality() )	
RESULT

CHECK( void setSecond(const FeatureType& sec) )
	DFeaturePair<2> p;
	const DFeature<2> f;
	p.setSecond(f);
	
	TEST_EQUAL( f, p.getSecond() )
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
