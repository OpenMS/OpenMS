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
#include <OpenMS/ANALYSIS/MAPMATCHING/Group.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap< Feature > ContainerType;
typedef ContainerType::value_type ElementType;
typedef Feature::PositionType PositionType;

START_TEST(Group, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Group<ContainerType>* ptr = 0;
CHECK(Group())
	ptr = new Group<ContainerType>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Group())
	delete ptr;
RESULT

CHECK(Group& operator= (const Group& source))
  ElementType e;
  IndexTuple it(1,2,e.getIntensity(),e.getPosition());
  Group<ContainerType> group;
  group.insert(it);

  Group<ContainerType> group_copy;
  group_copy = group;
  
  TEST_EQUAL(*(group.begin()) == *(group_copy.begin()), true)  
RESULT

CHECK(Group(const Group& source))
  ElementType e;
  IndexTuple it(1,2,e.getIntensity(),e.getPosition());
  Group<ContainerType> group;
  group.insert(it);

  Group<ContainerType> group_copy(group);
  
  TEST_EQUAL(*(group.begin()) == *(group_copy.begin()), true)  
RESULT

CHECK(bool isEmpty())
  ElementType e;
  IndexTuple it(1,2,e.getIntensity(),e.getPosition());
  Group<ContainerType> group;
  group.insert(it);

  TEST_EQUAL(group.isEmpty(), false) 
RESULT

CHECK(bool operator != (const Group& group) const)
  ElementType e1;
  IndexTuple it1(1,2,e1.getIntensity(),e1.getPosition());
  Group<ContainerType> group1;
  group1.insert(it1);
  
  ElementType e2;
  IndexTuple it2(1,2,e2.getIntensity(),e2.getPosition());
  Group<ContainerType> group2;
  group2.insert(it2);

  TEST_EQUAL(group1 != group2, false) 
RESULT

CHECK(bool operator == (const Group& group) const)
  ElementType e1;
  IndexTuple it1(1,2,e1.getIntensity(),e1.getPosition());
  Group<ContainerType> group1;
  group1.insert(it1);
  
  IndexTuple it2(1,2,e1.getIntensity(),e1.getPosition());
  Group<ContainerType> group2;
  group2.insert(it2);

  TEST_EQUAL(group1 == group2, true) 
RESULT

CHECK(unsigned int count() const)
  ElementType e;
  IndexTuple it(1,2,e.getIntensity(),e.getPosition());
  Group<ContainerType> group;
  group.insert(it);

  TEST_EQUAL(group.count() == 1, true) 
RESULT

CHECK((std::pair< typename Base::iterator, bool > insert(const Element& elem) throw(Exception::InvalidValue)))
  ElementType e;
  IndexTuple it1(1,2,e.getIntensity(),e.getPosition());
  Group<ContainerType> group;
  group.insert(it1);

  IndexTuple it2(1,3,e.getIntensity(),e.getPosition());
  
  TEST_EXCEPTION(Exception::InvalidValue,group.insert(it2));
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



