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

#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Tagging, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PRECISION(0.001)

// default ctor
Tagging* dv_ptr = 0;
CHECK(Tagging())
	dv_ptr = new Tagging;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK(~Tagging())
	delete dv_ptr;
RESULT

CHECK(const IsotopeVariant& getVariant() const)
	Tagging s;
	TEST_EQUAL(s.getVariant(),Tagging::LIGHT)
RESULT

CHECK(float getMassShift() const)
	Tagging s;
	TEST_REAL_EQUAL(s.getMassShift(),0.0)
RESULT

CHECK(void setMassShift(float mass_shift))
	Tagging s;
	s.setMassShift(4711.2);
	TEST_REAL_EQUAL(s.getMassShift(),4711.2)
RESULT

CHECK(void setVariant(const IsotopeVariant& variant))
	Tagging s;
	s.setVariant(Tagging::HEAVY);
	TEST_EQUAL(s.getVariant(),Tagging::HEAVY)
RESULT

//getType
CHECK([EXTRA] getType)
	Tagging s;
	TEST_EQUAL(s.getType(),"Tagging")
RESULT

//copy ctr
CHECK(Tagging(const Tagging&))
	Tagging s;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//copy
	Tagging s2(s);

	//get
	TEST_REAL_EQUAL(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_EQUAL(s2.getMass(),23.4)
RESULT

//assignment operator
CHECK(Tagging& operator=(const Tagging&))
	Tagging s,s2;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//assign
	s2 = s;

	//get
	TEST_REAL_EQUAL(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_EQUAL(s2.getMass(),23.4)
RESULT

//clone
CHECK(SampleTreatment* clone() const)
	Tagging s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Tagging* dp;
	
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Tagging*>(st);
	
	//get
	TEST_REAL_EQUAL(dp->getMassShift(),4711.2)
	TEST_EQUAL(dp->getVariant(),Tagging::LIGHT)
	TEST_REAL_EQUAL(dp->getMass(),23.4)
RESULT

CHECK(bool operator== (const SampleTreatment& rhs) const)
	Tagging empty,edit;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.setMassShift(4711.2);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setVariant(Tagging::HEAVY);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setMass(23.4);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);			

	edit.setMetaValue("color",string("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);	
	
	Modification m;
	TEST_EQUAL(m==empty, false);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
