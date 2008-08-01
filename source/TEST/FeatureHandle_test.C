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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/FeatureMap.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap< Feature > ContainerType;
typedef ContainerType::value_type ElementType;
typedef Feature::PositionType PositionType;

START_TEST(FeatureHandle, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureHandle* ptr = 0;
CHECK((FeatureHandle()))
	ptr = new FeatureHandle();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~FeatureHandle()))
	delete ptr;
RESULT

CHECK((FeatureHandle& operator=(const FeatureHandle &rhs)))
  ElementType e;
  FeatureHandle it(1,2,e);

  FeatureHandle it_copy;
  it_copy = it;
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
RESULT

CHECK((FeatureHandle(const FeatureHandle &rhs)))
  ElementType e;
  FeatureHandle it(1,2,e);

  FeatureHandle it_copy(it);
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
RESULT

CHECK((FeatureHandle(UInt map_index, UInt element_index, const Peak2D &point)))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
  TEST_EQUAL(it.getMapIndex() == 1, true)  
  TEST_EQUAL(it.getPosition() == e.getPosition(), true)  
RESULT

CHECK((virtual bool operator!=(const FeatureHandle &i) const))
  ElementType e;
  FeatureHandle it1(1,2,e);
  FeatureHandle it2(2,2,e);
  
  TEST_EQUAL(it1 != it2, true)  
RESULT

CHECK((virtual bool operator==(const FeatureHandle &i) const))
  ElementType e;
  FeatureHandle it1(2,2,e);
  FeatureHandle it2(2,2,e);
  
  TEST_EQUAL(it1 == it2, true)  
RESULT

CHECK((UInt getElementIndex() const))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK((UInt getMapIndex() const))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)  
RESULT

CHECK((void setElementIndex(UInt e)))
  FeatureHandle it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
RESULT

CHECK((void setMapIndex(UInt i)))
  FeatureHandle it;
  it.setMapIndex(2);

  TEST_EQUAL(it.getMapIndex() == 2, true)  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



