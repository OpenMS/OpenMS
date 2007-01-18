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

#include <OpenMS/VISUAL/MappingInfo.h>
#include <OpenMS/FORMAT/Param.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MappingInfo* ptr = 0;
CHECK(MappingInfo())
	ptr = new MappingInfo();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK([EXTRA] ~MappingInfo())
	delete ptr;
RESULT


CHECK(bool isMzToXAxis() const)
	MappingInfo m;
	TEST_EQUAL(m.isMzToXAxis(), true)
RESULT

CHECK(bool isXAxisAsc() const)
	MappingInfo m;
	TEST_EQUAL(m.isXAxisAsc(), true)	
RESULT

CHECK(bool isYAxisAsc() const)
		MappingInfo m;
	TEST_EQUAL(m.isXAxisAsc(), true)
RESULT


CHECK(void setMzToYAxis())
	MappingInfo m;
	m.setMzToYAxis();
	TEST_EQUAL(m.isMzToXAxis(), false)
RESULT

CHECK(void setMzToXAxis())
	MappingInfo m;
	m.setMzToYAxis();
	TEST_EQUAL(m.isMzToXAxis(), false)
	m.setMzToXAxis();
	TEST_EQUAL(m.isMzToXAxis(), true)
RESULT


CHECK(void setXAxisDesc())
	MappingInfo m;
	m.setXAxisDesc();
	TEST_EQUAL(m.isXAxisAsc(), false)
RESULT

CHECK(void setXAxisAsc())
	MappingInfo m;
	m.setXAxisDesc();
	TEST_EQUAL(m.isXAxisAsc(), false)
	m.setXAxisAsc();
	TEST_EQUAL(m.isXAxisAsc(), true)
RESULT


CHECK(void setYAxisDesc())
	MappingInfo m;
	m.setYAxisDesc();
	TEST_EQUAL(m.isYAxisAsc(), false)
RESULT

CHECK(void setYAxisAsc())
	MappingInfo m;
	m.setYAxisDesc();
	TEST_EQUAL(m.isYAxisAsc(), false)
	m.setYAxisAsc();
	TEST_EQUAL(m.isYAxisAsc(), true)
RESULT


CHECK(Param getParam())
	MappingInfo m;
	Param p;
	
	p=m.getParam();
	TEST_EQUAL((string)(p.getValue("MappingOfMzTo")),"X-Axis");
	TEST_EQUAL((string)(p.getValue("X-Axis-Orientation")),"Ascending");
	TEST_EQUAL((string)(p.getValue("Y-Axis-Orientation")),"Ascending");

	m.setMzToYAxis();
	m.setXAxisDesc();
	m.setYAxisDesc();
	p=m.getParam();
	TEST_EQUAL((string)(p.getValue("MappingOfMzTo")),"Y-Axis");
	TEST_EQUAL((string)(p.getValue("X-Axis-Orientation")),"Descending");
	TEST_EQUAL((string)(p.getValue("Y-Axis-Orientation")),"Descending");	
RESULT

CHECK(void setParam(const Param& p))
	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



