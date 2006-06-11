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
// $Id: Modification_test.C,v 1.3 2006/05/30 15:46:43 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Modification, "$Id: Modification_test.C,v 1.3 2006/05/30 15:46:43 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PRECISION(0.001)

// default ctor
Modification* dv_ptr = 0;
CHECK(Modification())
	dv_ptr = new Modification;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK(~Modification())
	delete dv_ptr;
RESULT

//basic accessors
CHECK(basic accessors)
	Modification s;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_EQUAL(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
RESULT

//getType
CHECK(getType)
	Modification s;
	TEST_EQUAL(s.getType(),"Modification")
RESULT

//copy ctr
CHECK(Modification(const Modification&))
	Modification s;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",string("red"));

	//copy
	Modification s2(s);

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_EQUAL(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
	
RESULT

//assignment operator
CHECK(operator=(const Modification&))
	Modification s,s2;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",string("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_EQUAL(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
RESULT

//clone
CHECK(operator=(const Modification&))
	Modification s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Modification* dp;
	
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",string("red"));

	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Modification*>(st);

	//get
	TEST_EQUAL(dp->getReagentName(),"TTEST")
	TEST_REAL_EQUAL(dp->getMass(),11.9)
	TEST_EQUAL(dp->getSpecificityType(),Modification::AA)
	TEST_EQUAL(dp->getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(string(dp->getMetaValue("color")),"red")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
