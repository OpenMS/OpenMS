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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/ANDIFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;

///////////////////////////

START_TEST(ANDIFile, "$Id$")

/////////////////////////////////////////////////////////////

ANDIFile* ptr = 0;
CHECK((ANDIFile()))
	ptr = new ANDIFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ANDIFile()))
	delete ptr;
RESULT

CHECK((template<typename MapType> void load(const String& filename, MapType& map) throw(Exception::FileNotFound, Exception::ParseError)))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	ANDIFile andi;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , andi.load("dummy/dummy.cdf",e) )

	// real test
	andi.load("data/ANDIFile_test.cdf",e);
  //---------------------------------------------------------------------------
  // 60 : (120,100) 
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100) 
	//---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 3)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[0].getRetentionTime(), 60)
	TEST_REAL_EQUAL(e[1].getRetentionTime(), 120)
	TEST_REAL_EQUAL(e[2].getRetentionTime(), 180)
	TEST_REAL_EQUAL(e[0].getContainer().size(), 1)
	TEST_REAL_EQUAL(e[1].getContainer().size(), 3)
	TEST_REAL_EQUAL(e[2].getContainer().size(), 5)

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getIntensity(), 100)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
