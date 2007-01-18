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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/IndexTuple.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DFeatureMap<2, DFeature<2, KernelTraits> > ContainerType;
typedef ContainerType::value_type ElementType;
typedef ElementType::TraitsType TraitsType;
typedef DPosition<2, TraitsType> PositionType;

START_TEST(IndexTuple, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexTuple<>* ptr = 0;
CHECK(IndexTuple())
	ptr = new IndexTuple<>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~IndexTuple())
	delete ptr;
RESULT

CHECK(IndexTuple& operator = (const IndexTuple& source))
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  IndexTuple<ContainerType> it_copy;
  it_copy = it;
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(&(it.getElement()) == &(it_copy.getElement()), true)  
  TEST_EQUAL(it.getTransformedPosition() == it_copy.getTransformedPosition(), true)  
RESULT

CHECK(IndexTuple(const IndexTuple& source))
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  IndexTuple<ContainerType> it_copy(it);
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(&(it.getElement()) == &(it_copy.getElement()), true)  
  TEST_EQUAL(it.getTransformedPosition() == it_copy.getTransformedPosition(), true)  
RESULT

CHECK((IndexTuple(const UnsignedInt& map_index, const UnsignedInt& element_index, const ElementType& element)))
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
  TEST_EQUAL(it.getMapIndex() == 1, true)  
  TEST_EQUAL(&(it.getElement()) == &e, true)  
  TEST_EQUAL(it.getTransformedPosition() == e.getPosition(), true)  
RESULT

CHECK(bool operator != (const IndexTuple& i) const)
  ElementType e;
  IndexTuple<ContainerType> it1(1,2,e);
  IndexTuple<ContainerType> it2(2,2,e);
  
  TEST_EQUAL(it1 != it2, true)  
RESULT

CHECK(bool operator == (const IndexTuple& i) const)
  ElementType e;
  IndexTuple<ContainerType> it1(2,2,e);
  IndexTuple<ContainerType> it2(2,2,e);
  
  TEST_EQUAL(it1 == it2, true)  
RESULT

CHECK(const ElementType& getElement() const)
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  TEST_EQUAL(&(it.getElement()) == &e, true)  
RESULT

CHECK(const PositionType& getTransformedPosition() const)
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  TEST_EQUAL(it.getTransformedPosition() == e.getPosition(), true)  
RESULT

CHECK(const UnsignedInt& getElementIndex() const)
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK(const UnsignedInt& getMapIndex() const)
  ElementType e;
  IndexTuple<ContainerType> it(1,2,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)  
RESULT

CHECK(void setElement(const ElementType& e))
  ElementType e;
  IndexTuple<ContainerType> it;
  it.setElement(e);

  TEST_EQUAL(&(it.getElement()) == &e, true)   
RESULT

CHECK(void setElementIndex(const UnsignedInt& e))
  IndexTuple<ContainerType> it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK(void setMapIndex(const UnsignedInt& c))
  IndexTuple<ContainerType> it;
  it.setMapIndex(2);

  TEST_EQUAL(it.getMapIndex() == 2, true)  
RESULT

CHECK(void setTransformedPosition(const PositionType& p))
  ElementType e;
  IndexTuple<ContainerType> it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getTransformedPosition() == e.getPosition(), true)  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



