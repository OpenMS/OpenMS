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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/KERNEL/Feature.h>

///////////////////////////

START_TEST(ElementPair<Feature>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

ElementPair<Feature>* ptr = 0;
CHECK((ElementPair()))
	ptr = new ElementPair<Feature>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT


CHECK((virtual ~ElementPair()))
	delete ptr;
RESULT


CHECK((ElementPair(const ElementPair& fp)))
	ElementPair<Feature> p1;
	p1.setQuality(5.0);
	
	ElementPair<Feature> p2(p1);
	
	TEST_REAL_EQUAL( p1.getQuality(),p2.getQuality() )
RESULT


CHECK((ElementPair(const ElementType& first, const ElementType& second, DoubleReal quality=0.0)))
	Feature f1;
	Feature f2;
	
	ElementPair<Feature> pair(f1,f2);
	
	TEST_EQUAL( f1, pair.getFirst() )
	TEST_EQUAL( f2, pair.getSecond() )
RESULT


CHECK((ElementPair& operator = (const ElementPair& rhs)))
	ElementPair<Feature> p1;
	p1.setQuality( 5.0 );
	
	ElementPair<Feature> p2 = p1;
	
	TEST_REAL_EQUAL( p1.getQuality(),p2.getQuality() )
RESULT


CHECK((bool operator == (const ElementPair& rhs) const))
	ElementPair<Feature> p1;
	Feature f1;
  f1.setRT(1.0);
  f1.setMZ(2.0);
	Feature f2;
  f2.setRT(3.0);
  f2.setMZ(4.0);
	
	p1.setFirst(f1);
	p1.setSecond(f2);
	p1.setQuality(5.0);
  
	ElementPair<Feature> p2;
	Feature f3;
  f3.setRT(1.0);
  f3.setMZ(2.0);
  Feature f4;
  f4.setRT(3.0);
  f4.setMZ(4.0);
	
	p2.setFirst(f3);
	p2.setSecond(f4);
	p2.setQuality(5.0);
  
  TEST_EQUAL(p1==p2,true);
RESULT


CHECK((bool operator != (const ElementPair& rhs) const))
	ElementPair<Feature> p1;
	Feature f1;
  f1.setRT(2.0);
  f1.setMZ(2.0);
  Feature f2;
  f2.setRT(2.0);
  f2.setMZ(2.0);
	
	p1.setFirst(f1);
	p1.setSecond(f2);
	
	ElementPair<Feature> p2;
	Feature f3;
  f1.setRT(1.0);
  f1.setMZ(1.0);
  Feature f4;
  f2.setRT(1.0);
  f2.setMZ(1.0);
	
	p2.setFirst(f3);
	p2.setSecond(f4);
	
	TEST_EQUAL(p1!=p2,true);
RESULT


CHECK((DoubleReal getQuality() const))
	ElementPair<Feature> p;
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
RESULT


CHECK((void setQuality(DoubleReal quality)))
	ElementPair<Feature> p;
	p.setQuality(123.456);
	TEST_REAL_EQUAL(p.getQuality(), 123.456)
	p.setQuality(-0.12345);
	TEST_REAL_EQUAL(p.getQuality(), -0.12345)
	p.setQuality(0.0);
	TEST_REAL_EQUAL(p.getQuality(), 0.0)
RESULT


CHECK((const ElementType& getFirst() const))
	ElementPair<Feature> p;
	
	Feature f1;
  f1.setRT(1.0);
  f1.setMZ(2.0);
  p.setFirst(f1);

	const Feature f2 = p.getFirst();
	TEST_EQUAL(f1,f2)

RESULT


CHECK((void setFirst(const ElementType &frt)))
	ElementPair<Feature> p;
	const Feature f;
	p.setFirst(f);
	
	TEST_EQUAL( f, p.getFirst() )
RESULT


CHECK((const ElementType& getSecond() const))
	ElementPair<Feature> p;
	
	Feature f1;
  f1.setRT(1.0);
  f1.setMZ(2.0);
  p.setSecond(f1);

	const Feature f2 = p.getSecond();
	TEST_EQUAL(f1,f2)
RESULT


CHECK((void setSecond(const ElementType &sec)))
	ElementPair<Feature> p;
	const Feature f;
	p.setSecond(f);
	
	TEST_EQUAL( f, p.getSecond() )
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
