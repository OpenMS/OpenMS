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
#include <OpenMS/ANALYSIS/MAPMATCHING/IndexTuple.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap< Feature > ContainerType;
typedef ContainerType::value_type ElementType;
typedef Feature::PositionType PositionType;

START_TEST(IndexTuple, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexTuple* ptr = 0;
CHECK(IndexTuple())
	ptr = new IndexTuple();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~IndexTuple())
	delete ptr;
RESULT

CHECK(IndexTuple& operator = (const IndexTuple& source))
  ElementType e;
  IndexTuple it(1,2,e);

  IndexTuple it_copy;
  it_copy = it;
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
RESULT

CHECK(IndexTuple(const IndexTuple& source))
  ElementType e;
  IndexTuple it(1,2,e);

  IndexTuple it_copy(it);
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
RESULT

CHECK((IndexTuple(UInt map_index, UInt element_index, const ElementType& element)))
  ElementType e;
  IndexTuple it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
  TEST_EQUAL(it.getMapIndex() == 1, true)  
  TEST_EQUAL(it.getPosition() == e.getPosition(), true)  
RESULT

CHECK(bool operator != (const IndexTuple& i) const)
  ElementType e;
  IndexTuple it1(1,2,e);
  IndexTuple it2(2,2,e);
  
  TEST_EQUAL(it1 != it2, true)  
RESULT

CHECK(bool operator == (const IndexTuple& i) const)
  ElementType e;
  IndexTuple it1(2,2,e);
  IndexTuple it2(2,2,e);
  
  TEST_EQUAL(it1 == it2, true)  
RESULT

CHECK(const PositionType& getPosition() const)
  ElementType e;
  IndexTuple it(1,2,e);

  TEST_EQUAL(it.getPosition() == e.getPosition(), true)  
RESULT

CHECK(UInt getElementIndex() const)
  ElementType e;
  IndexTuple it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK(UInt getMapIndex() const)
  ElementType e;
  IndexTuple it(1,2,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)  
RESULT

CHECK(void setElementIndex(UInt e))
  IndexTuple it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK(void setMapIndex(UInt c))
  IndexTuple it;
  it.setMapIndex(2);

  TEST_EQUAL(it.getMapIndex() == 2, true)  
RESULT

CHECK(void setPosition(const PositionType& p))
  ElementType e;
  IndexTuple it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getPosition() == e.getPosition(), true)  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



