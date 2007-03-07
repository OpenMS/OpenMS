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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummySeeder.h>

using namespace OpenMS;

START_TEST(DummySeeder, "$Id: DummySeeder_test.C 1586 2007-03-01 17:59:10Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
	
using namespace OpenMS;
using namespace std;

// default ctor
DummySeeder* ptr = 0;
CHECK(DummySeeder())
	ptr = new DummySeeder();
  TEST_EQUAL(ptr->getName(), "DummySeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DummySeeder())
	delete ptr;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
