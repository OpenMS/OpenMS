// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Modification* dv_ptr = nullptr;
Modification* dv_nullPointer = nullptr;
START_SECTION((Modification()))
	dv_ptr = new Modification;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Modification()))
	delete dv_ptr;
END_SECTION

START_SECTION((const String& getReagentName() const))
	Modification s;
	TEST_EQUAL(s.getReagentName(),"")
END_SECTION

START_SECTION((double getMass() const ))
	Modification s;
	TEST_REAL_SIMILAR(s.getMass(),0.0)
END_SECTION

START_SECTION((const SpecificityType& getSpecificityType() const))
	Modification s;
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
END_SECTION

START_SECTION((const String& getAffectedAminoAcids() const))
	Modification s;
	TEST_EQUAL(s.getAffectedAminoAcids(),"")
END_SECTION

START_SECTION((void setReagentName(const String& reagent_name)))
	Modification s;
	s.setReagentName("TTEST");
	TEST_EQUAL(s.getReagentName(),"TTEST")
END_SECTION

START_SECTION((void setMass(double mass)))
	Modification s;
	s.setMass(11.9);
	TEST_REAL_SIMILAR(s.getMass(),11.9)
END_SECTION

START_SECTION((void setSpecificityType(const SpecificityType& specificity_type)))
	Modification s;
	s.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(s.getSpecificityType(),Modification::CTERM)
END_SECTION

START_SECTION((void setAffectedAminoAcids(const String& affected_amino_acids)))
	Modification s;
	s.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Modification s;
	TEST_EQUAL(s.getType(),"Modification")
END_SECTION

//copy ctr
START_SECTION((Modification(const Modification&)))
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
	TEST_REAL_SIMILAR(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
	
END_SECTION

//assignment operator
START_SECTION((Modification& operator=(const Modification&)))
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
	TEST_REAL_SIMILAR(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
END_SECTION

//clone
START_SECTION((virtual SampleTreatment* clone() const ))
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
	TEST_REAL_SIMILAR(dp->getMass(),11.9)
	TEST_EQUAL(dp->getSpecificityType(),Modification::AA)
	TEST_EQUAL(dp->getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(dp->getMetaValue("color")),"red")

	delete st;
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Modification empty,edit;
	
	TEST_TRUE(edit == empty);

	edit.setMass(11.9);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);

	edit.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);		

	edit.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);			

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);	
	
	Tagging m;
	TEST_EQUAL(m==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
