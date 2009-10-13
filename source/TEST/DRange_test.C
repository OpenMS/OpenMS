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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DRange.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;

START_TEST(DRange<D>, "$id$")

/////////////////////////////////////////////////////////////

//do not modify these points, they are used in many tests
DPosition<2> p1,p2,p3,one,two;
p1[0]=-1.0f;
p1[1]=-2.0f;
p2[0]=3.0f;
p2[1]=4.0f;	
p3[0]=-10.0f;
p3[1]=20.0f;
one[0]=1;
one[1]=1;
two[0]=2;
two[1]=2;

//do not modify these points, they are used in many tests

std::cout.precision(writtenDigits<>(DoubleReal()));
std::cerr.precision(writtenDigits<>(DoubleReal()));

DRange<2>* ptr = 0;
START_SECTION(DRange())
	ptr = new DRange<2>;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~DRange())
	delete ptr;
END_SECTION

START_SECTION(DRange(const PositionType& lower, const PositionType& upper))
	DRange<2> r(p1,p2);
	TEST_REAL_SIMILAR(r.min()[0],-1.0f);
	TEST_REAL_SIMILAR(r.min()[1],-2.0f);
	TEST_REAL_SIMILAR(r.max()[0],3.0f);
	TEST_REAL_SIMILAR(r.max()[1],4.0f);
END_SECTION

//do not modify this range, it is used in many tests
DRange<2> r(p1,p2);
//do not modify this range, it is used in many tests

START_SECTION(DRange(const DRange& range))
	DRange<2> r2(r);
	TEST_REAL_SIMILAR(r2.min()[0],-1.0f);
	TEST_REAL_SIMILAR(r2.min()[1],-2.0f);
	TEST_REAL_SIMILAR(r2.max()[0],3.0f);
	TEST_REAL_SIMILAR(r2.max()[1],4.0f);
END_SECTION

START_SECTION(DRange(const Base& range))
	Internal::DIntervalBase<2> ib(r);
	DRange<2> r2(ib);
	TEST_REAL_SIMILAR(r2.min()[0],-1.0f);
	TEST_REAL_SIMILAR(r2.min()[1],-2.0f);
	TEST_REAL_SIMILAR(r2.max()[0],3.0f);
	TEST_REAL_SIMILAR(r2.max()[1],4.0f);
END_SECTION

START_SECTION(DRange& operator=(const Base& rhs))
	Internal::DIntervalBase<2> ib(r);
	DRange<2> r2;
	r2 = ib;
	TEST_REAL_SIMILAR(r2.min()[0],-1.0f);
	TEST_REAL_SIMILAR(r2.min()[1],-2.0f);
	TEST_REAL_SIMILAR(r2.max()[0],3.0f);
	TEST_REAL_SIMILAR(r2.max()[1],4.0f);
END_SECTION

START_SECTION(DRange& operator=(const DRange& rhs))
	DRange<2> r2;
	r2 = r;
	TEST_REAL_SIMILAR(r2.min()[0],-1.0f);
	TEST_REAL_SIMILAR(r2.min()[1],-2.0f);
	TEST_REAL_SIMILAR(r2.max()[0],3.0f);
	TEST_REAL_SIMILAR(r2.max()[1],4.0f);
END_SECTION

START_SECTION(DRange(CoordinateType minx, CoordinateType miny, CoordinateType maxx, CoordinateType maxy))
	DRange<2> r2(1.0f,2.0f,3.0f,4.0f);
	TEST_REAL_SIMILAR(r2.min()[0],1.0f);
	TEST_REAL_SIMILAR(r2.min()[1],2.0f);
	TEST_REAL_SIMILAR(r2.max()[0],3.0f);
	TEST_REAL_SIMILAR(r2.max()[1],4.0f);
END_SECTION

START_SECTION(bool operator == (const DRange& rhs) const )
	DRange<2> r2(r);
	TEST_EQUAL(r==r2,true);
	r2.setMinX(0.0f);
	TEST_EQUAL(r==r2,false);
	r2.setMinX(r.min()[0]);
	TEST_EQUAL(r==r2,true);
	r2.setMaxY(0.0f);
	TEST_EQUAL(r==r2,false);
	r2.setMaxY(r.max()[1]);
	TEST_EQUAL(r==r2,true);
END_SECTION

START_SECTION(bool operator == (const Base& rhs) const )
	Internal::DIntervalBase<2> r2(r);
	TEST_EQUAL(r==r2,true);
	r2.setMinX(0.0f);
	TEST_EQUAL(r==r2,false);
	r2.setMinX(r.min()[0]);
	TEST_EQUAL(r==r2,true);
	r2.setMaxY(0.0f);
	TEST_EQUAL(r==r2,false);
	r2.setMaxY(r.max()[1]);
	TEST_EQUAL(r==r2,true);
END_SECTION

START_SECTION(bool encloses(const PositionType& position) const)
	DRange<2> r2(p1,p2);
	DPosition<2> p;
	p[0]=0.0f;
	p[1]=0.0f;
	TEST_EQUAL(r2.encloses(p),true);
	p[0]=-3.0f;
	p[1]=-3.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=-3.0f;
	p[1]=0.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=0.0f;
	p[1]=-3.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=-3.0f;
	p[1]=5.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=0.0f;
	p[1]=5.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=5.0f;
	p[1]=5.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=5.0f;
	p[1]=0.0f;
	TEST_EQUAL(r2.encloses(p),false);
	p[0]=5.0f;
	p[1]=-3.0f;
	TEST_EQUAL(r2.encloses(p),false);
END_SECTION

START_SECTION(DRangeIntersection intersects(const DRange& range) const)
	DRange<2> r2(p1,p2);
	DRange<2> r3(r2);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Inside)
	r3.setMaxX(10.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)
	r3.setMax(r2.max()+one);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)
	r3.setMin(r2.max()+one);
	r3.setMax(r2.max()+two);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMin(r2.min());
	r3.setMinX(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)		
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(1.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)		
	r3.setMinX(-10.0f);
	r3.setMinY(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)	
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)	
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)	
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)	
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Disjoint)
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)	
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)	
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.intersects(r3),DRange<2>::Intersects)		
END_SECTION

START_SECTION(bool isIntersected(const DRange& range) const)
	DRange<2> r2(p1,p2);
	DRange<2> r3(r2);
	TEST_EQUAL(r2.isIntersected(r3),true)
	r3.setMaxX(10.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)
	r3.setMax(r2.max()+one);
	TEST_EQUAL(r2.isIntersected(r3),true)
	r3.setMin(r2.max()+one);
	r3.setMax(r2.max()+two);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMin(r2.min());
	r3.setMinX(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(-10.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(-9.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)		
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(1.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)		
	r3.setMinX(-10.0f);
	r3.setMinY(10.0f);
	r3.setMax(r3.min()+one);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMinX(-10.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(-9.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)	
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)	
	r3.setMinX(9.0f);
	r3.setMinY(0.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(10.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)	
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)	
	r3.setMinX(9.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(10.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.isIntersected(r3),false)
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)	
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)	
	r3.setMinX(-5.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(0.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(0.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)	
	r3.setMinX(0.0f);
	r3.setMinY(-5.0f);
	r3.setMaxX(5.0f);
	r3.setMaxY(5.0f);
	TEST_EQUAL(r2.isIntersected(r3),true)		
END_SECTION

START_SECTION(DRange united(const DRange<D>& other_range) const)
	DRange<2> r2(p1,p2);
	DRange<2> r3(r2);
	TEST_EQUAL(r2 == r2.united(r3), true)
	TEST_EQUAL(r3 == r2.united(r3), true)
	TEST_EQUAL(r2 == r3.united(r2), true)
	TEST_EQUAL(r3 == r3.united(r2), true)
	r3.setMin(r2.max()+one);
	r3.setMax(r2.max()+two);
	DRange<2> r4;
	r4.setMin(r2.min());
	r4.setMax(r3.max());
	TEST_EQUAL(r2.united(r3) == r4, true);
	TEST_EQUAL(r3.united(r2) == r4, true);
END_SECTION


START_SECTION(bool encloses(CoordinateType x, CoordinateType y) const)
	DRange<2> r2(p1,p2);
	TEST_EQUAL(r2.encloses(0.0f,0.0f),true);
	TEST_EQUAL(r2.encloses(-3.0f,-3.0f),false);
	TEST_EQUAL(r2.encloses(-3.0f,0.0f),false);
	TEST_EQUAL(r2.encloses(0.0f,-3.0f),false);
	TEST_EQUAL(r2.encloses(-3.0f,5.0f),false);
	TEST_EQUAL(r2.encloses(0.0f,5.0f),false);
	TEST_EQUAL(r2.encloses(5.0f,5.0f),false);
	TEST_EQUAL(r2.encloses(5.0f,0.0f),false);
	TEST_EQUAL(r2.encloses(5.0f,-3.0f),false);
END_SECTION

START_SECTION(bool isEmpty() const)
	DRange<2> tmp;
	TEST_EQUAL(tmp.isEmpty(), true);
	tmp = DRange<2>::zero;
	TEST_EQUAL(tmp.isEmpty(), true);
	tmp = DRange<2>(p1,p2);		
	TEST_EQUAL(tmp.isEmpty(), false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
