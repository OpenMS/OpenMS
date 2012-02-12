// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Map.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Map, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Map<int, int>* map_ptr = 0;
Map<int, int>* map_nullPointer = 0;
START_SECTION((Map()))
	map_ptr = new Map<int, int>;
  TEST_NOT_EQUAL(map_ptr, map_nullPointer)
END_SECTION

START_SECTION((~Map()))
	delete map_ptr;
END_SECTION

START_SECTION((T& operator [] (const Key& key)))
	Map<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	TEST_EQUAL(hm.size(), 6)
	TEST_EQUAL(hm[0], 1)
	TEST_EQUAL(hm[1], 2)
	TEST_EQUAL(hm[2], 4)
	TEST_EQUAL(hm[3], 8)
	TEST_EQUAL(hm[4], 16)
	TEST_EQUAL(hm[5], 32)
END_SECTION

START_SECTION((const T & operator[](const Key &key) const ))
	Map<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	const Map<int, int>& const_map = const_cast<const Map<int, int>&>(hm);
	TEST_EQUAL(const_map.size(), 6)
	TEST_EQUAL(const_map[0], 1)
	TEST_EQUAL(const_map[1], 2)
	TEST_EQUAL(const_map[2], 4)
	TEST_EQUAL(const_map[3], 8)
	TEST_EQUAL(const_map[4], 16)
	TEST_EQUAL(const_map[5], 32)
	typedef Map<int,int> MyMap; // otherwise next line wont work
	TEST_EXCEPTION(MyMap::IllegalKey, const_map[6])
END_SECTION

START_SECTION((bool has(const Key& key) const))
	Map<int, int> hm;
	hm.insert(Map<int, int>::ValueType(0, 0));
	hm.insert(Map<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm.has(0), true)
	TEST_EQUAL(hm.has(1), true)
	TEST_EQUAL(hm.has(2), false)
END_SECTION

START_SECTION(([Map::IllegalKey] IllegalKey(const char *file, int line, const char *function)))
	// already tested in const T & operator[](const Key &key) const
	NOT_TESTABLE
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
