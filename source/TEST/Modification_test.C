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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <sstream>

///////////////////////////

START_TEST(Modification, "$Id$")

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

CHECK(const String& getReagentName() const)
	Modification s;
	TEST_EQUAL(s.getReagentName(),"")
RESULT

CHECK(float getMass() const)
	Modification s;
	TEST_REAL_EQUAL(s.getMass(),0.0)
RESULT

CHECK(const SpecificityType& getSpecificityType() const)
	Modification s;
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
RESULT

CHECK(const String& getAffectedAminoAcids() const)
	Modification s;
	TEST_EQUAL(s.getAffectedAminoAcids(),"")
RESULT

CHECK(void setReagentName(const String& reagent_name))
	Modification s;
	s.setReagentName("TTEST");
	TEST_EQUAL(s.getReagentName(),"TTEST")
RESULT

CHECK(void setMass(float mass))
	Modification s;
	s.setMass(11.9);
	TEST_REAL_EQUAL(s.getMass(),11.9)
RESULT

CHECK(void setSpecificityType(const SpecificityType& specificity_type))
	Modification s;
	s.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(s.getSpecificityType(),Modification::CTERM)
RESULT

CHECK(void setAffectedAminoAcids(const String& affected_amino_acids))
	Modification s;
	s.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
RESULT

//getType
CHECK([EXTRA] getType)
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
	s.setMetaValue("color",String("red"));

	//copy
	Modification s2(s);

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_EQUAL(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
	
RESULT

//assignment operator
CHECK(Modification& operator=(const Modification&))
	Modification s,s2;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",String("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_EQUAL(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
RESULT

//clone
CHECK(SampleTreatment* clone() const)
	Modification s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Modification* dp;
	
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",String("red"));

	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Modification*>(st);

	//get
	TEST_EQUAL(dp->getReagentName(),"TTEST")
	TEST_REAL_EQUAL(dp->getMass(),11.9)
	TEST_EQUAL(dp->getSpecificityType(),Modification::AA)
	TEST_EQUAL(dp->getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(dp->getMetaValue("color")),"red")
RESULT

CHECK(bool operator== (const SampleTreatment& rhs) const)
	Modification empty,edit;
	
	TEST_EQUAL(edit==empty, true);

	edit.setMass(11.9);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);			

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);	
	
	Tagging m;
	TEST_EQUAL(m==empty, false);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
