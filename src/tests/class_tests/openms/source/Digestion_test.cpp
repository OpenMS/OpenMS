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

#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Digestion, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Digestion* dv_ptr = nullptr;
Digestion* dv_nullPointer = nullptr;
START_SECTION((Digestion()))
	dv_ptr = new Digestion;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Digestion()))
	delete dv_ptr;
END_SECTION

//basic accessors
START_SECTION((const String& getEnzyme() const))
	Digestion s;
	TEST_EQUAL(s.getEnzyme(),"")
END_SECTION

//basic accessors
START_SECTION((double getDigestionTime() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getDigestionTime(),0.0)
END_SECTION

//basic accessors
START_SECTION((double getTemperature() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getTemperature(),0.0)
END_SECTION

//basic accessors
START_SECTION((double getPh() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getPh(),0.0)
END_SECTION

//basic accessors
START_SECTION((void setEnzyme(const String& enzyme)))
	Digestion s;
	s.setEnzyme("TTEST");
	TEST_EQUAL(s.getEnzyme(),"TTEST")
END_SECTION

//basic accessors
START_SECTION((void setDigestionTime(double digestion_time)))
	Digestion s;
	//set
	s.setDigestionTime(4711.2);
	TEST_REAL_SIMILAR(s.getDigestionTime(),4711.2)
END_SECTION

//basic accessors
START_SECTION((void setTemperature(double temperature)))
	Digestion s;
	s.setTemperature(4711.3);
	TEST_REAL_SIMILAR(s.getTemperature(),4711.3)
END_SECTION

//basic accessors
START_SECTION((void setPh(double ph)))
	Digestion s;
	s.setPh(4711.4);
	TEST_REAL_SIMILAR(s.getPh(),4711.4)
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Digestion s;
	TEST_EQUAL(s.getType(),"Digestion")
END_SECTION

//copy ctr
START_SECTION((Digestion(const Digestion&)))
	Digestion s;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));
	
	//copy
	Digestion s2(s);

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(s2.getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(s2.getTemperature(),4711.3)
	TEST_REAL_SIMILAR(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
END_SECTION

START_SECTION((Digestion& operator=(const Digestion&)))
	Digestion s,s2;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(s2.getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(s2.getTemperature(),4711.3)
	TEST_REAL_SIMILAR(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
END_SECTION

START_SECTION((virtual SampleTreatment* clone() const ))
	Digestion s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Digestion* dp;
	
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));
	
	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Digestion*>(st);
	
	//get
	TEST_EQUAL(dp->getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(dp->getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(dp->getTemperature(),4711.3)
	TEST_REAL_SIMILAR(dp->getPh(),4711.4)
	TEST_EQUAL(string(dp->getMetaValue("color")),"red")
	delete dp;
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Digestion empty,edit;
	
	TEST_TRUE(edit == empty);
	
	edit.setEnzyme("TTEST");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);

	edit.setDigestionTime(4711.2);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);		

	edit.setTemperature(4711.3);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_TRUE(edit == empty);			

	edit.setPh(4711.4);
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
