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
START_SECTION((FeatureHandle()))
	ptr = new FeatureHandle();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureHandle()))
	delete ptr;
END_SECTION

START_SECTION((FeatureHandle& operator=(const FeatureHandle &rhs)))
  ElementType e;
  FeatureHandle it(1,2,e);

  FeatureHandle it_copy;
  it_copy = it;
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
END_SECTION

START_SECTION((FeatureHandle(const FeatureHandle &rhs)))
  ElementType e;
  FeatureHandle it(1,2,e);

  FeatureHandle it_copy(it);
  
  TEST_EQUAL(it.getElementIndex() == it_copy.getElementIndex(), true)  
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)  
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)  
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)  
END_SECTION

START_SECTION((FeatureHandle(UInt map_index, UInt element_index, const Peak2D &point)))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
  TEST_EQUAL(it.getMapIndex() == 1, true)  
  TEST_EQUAL(it.getPosition() == e.getPosition(), true)  
END_SECTION

START_SECTION((virtual bool operator!=(const FeatureHandle &i) const))
  ElementType e;
  FeatureHandle it1(1,2,e);
  FeatureHandle it2(2,2,e);
  
  TEST_EQUAL(it1 != it2, true)  
END_SECTION

START_SECTION((virtual bool operator==(const FeatureHandle &i) const))
  ElementType e;
  FeatureHandle it1(2,2,e);
  FeatureHandle it2(2,2,e);
  
  TEST_EQUAL(it1 == it2, true)  
END_SECTION

START_SECTION((UInt getElementIndex() const))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
END_SECTION

START_SECTION((UInt getMapIndex() const))
  ElementType e;
  FeatureHandle it(1,2,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)  
END_SECTION

START_SECTION((void setElementIndex(UInt e)))
  FeatureHandle it;
  it.setElementIndex(2);

  TEST_EQUAL(it.getElementIndex() == 2, true)  
END_SECTION

START_SECTION((void setMapIndex(UInt i)))
  FeatureHandle it;
  it.setMapIndex(2);

  TEST_EQUAL(it.getMapIndex() == 2, true)  
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



