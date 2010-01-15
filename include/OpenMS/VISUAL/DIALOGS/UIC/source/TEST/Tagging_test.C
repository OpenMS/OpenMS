// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Tagging, "$Id: Tagging_test.C 6135 2009-10-19 16:05:59Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Tagging* dv_ptr = 0;
START_SECTION((Tagging()))
	dv_ptr = new Tagging;
	TEST_NOT_EQUAL(dv_ptr, 0)
END_SECTION

// destructor
START_SECTION((virtual ~Tagging()))
	delete dv_ptr;
END_SECTION

START_SECTION((const IsotopeVariant& getVariant() const))
	Tagging s;
	TEST_EQUAL(s.getVariant(),Tagging::LIGHT)
END_SECTION

START_SECTION((DoubleReal getMassShift() const ))
	Tagging s;
	TEST_REAL_SIMILAR(s.getMassShift(),0.0)
END_SECTION

START_SECTION((void setMassShift(DoubleReal mass_shift)))
	Tagging s;
	s.setMassShift(4711.2);
	TEST_REAL_SIMILAR(s.getMassShift(),4711.2)
END_SECTION

START_SECTION((void setVariant(const IsotopeVariant& variant)))
	Tagging s;
	s.setVariant(Tagging::HEAVY);
	TEST_EQUAL(s.getVariant(),Tagging::HEAVY)
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Tagging s;
	TEST_EQUAL(s.getType(),"Tagging")
END_SECTION

//copy ctr
START_SECTION((Tagging(const Tagging&)))
	Tagging s;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//copy
	Tagging s2(s);

	//get
	TEST_REAL_SIMILAR(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(s2.getMass(),23.4)
END_SECTION

//assignment operator
START_SECTION((Tagging& operator=(const Tagging&)))
	Tagging s,s2;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//assign
	s2 = s;

	//get
	TEST_REAL_SIMILAR(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(s2.getMass(),23.4)
END_SECTION

//clone
START_SECTION((virtual SampleTreatment* clone() const ))
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
	TEST_REAL_SIMILAR(dp->getMassShift(),4711.2)
	TEST_EQUAL(dp->getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(dp->getMass(),23.4)
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
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

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);	
	
	Modification m;
	TEST_EQUAL(m==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
