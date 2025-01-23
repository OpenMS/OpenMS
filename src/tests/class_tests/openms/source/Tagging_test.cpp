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

#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Tagging, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Tagging* dv_ptr = nullptr;
Tagging* dv_nullPointer = nullptr;
START_SECTION((Tagging()))
	dv_ptr = new Tagging;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Tagging()))
	delete dv_ptr;
END_SECTION

START_SECTION((const IsotopeVariant& getVariant() const))
	Tagging s;
	TEST_EQUAL(s.getVariant(),Tagging::LIGHT)
END_SECTION

START_SECTION((double getMassShift() const ))
	Tagging s;
	TEST_REAL_SIMILAR(s.getMassShift(),0.0)
END_SECTION

START_SECTION((void setMassShift(double mass_shift)))
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
	delete st;
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Tagging empty,edit;
	
	TEST_TRUE(edit == empty);
	
	edit.setMassShift(4711.2);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);

	edit.setVariant(Tagging::HEAVY);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);		

	edit.setMass(23.4);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);			

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);	
	
	Modification m;
	TEST_EQUAL(m==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
