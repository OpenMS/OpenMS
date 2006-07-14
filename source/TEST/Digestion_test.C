// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Digestion.h>
#include <sstream>

///////////////////////////

START_TEST(Digestion, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PRECISION(0.001)

// default ctor
Digestion* dv_ptr = 0;
CHECK(Digestion())
	dv_ptr = new Digestion;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK(~Digestion())
	delete dv_ptr;
RESULT

//basic accessors
CHECK(basic accessors)
	Digestion s;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	//get
	TEST_EQUAL(s.getEnzyme(),"TTEST")
	TEST_REAL_EQUAL(s.getDigestionTime(),4711.2)
	TEST_REAL_EQUAL(s.getTemperature(),4711.3)
	TEST_REAL_EQUAL(s.getPh(),4711.4)
RESULT

//getType
CHECK(getType)
	Digestion s;
	TEST_EQUAL(s.getType(),"Digestion")
RESULT

//copy ctr
CHECK(Digestion(const Digestion&))
	Digestion s;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",string("red"));
	
	//copy
	Digestion s2(s);

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_EQUAL(s2.getDigestionTime(),4711.2)
	TEST_REAL_EQUAL(s2.getTemperature(),4711.3)
	TEST_REAL_EQUAL(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
RESULT

//assignment operator
CHECK(operator=(const Digestion&))
	Digestion s,s2;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",string("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_EQUAL(s2.getDigestionTime(),4711.2)
	TEST_REAL_EQUAL(s2.getTemperature(),4711.3)
	TEST_REAL_EQUAL(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
RESULT

//clone
CHECK(operator=(const Digestion&))
	Digestion s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Digestion* dp;
	
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",string("red"));
	
	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Digestion*>(st);
	
	//get
	TEST_EQUAL(dp->getEnzyme(),"TTEST")
	TEST_REAL_EQUAL(dp->getDigestionTime(),4711.2)
	TEST_REAL_EQUAL(dp->getTemperature(),4711.3)
	TEST_REAL_EQUAL(dp->getPh(),4711.4)
	TEST_EQUAL(string(dp->getMetaValue("color")),"red")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
