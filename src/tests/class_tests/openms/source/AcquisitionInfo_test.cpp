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
#include <OpenMS/METADATA/AcquisitionInfo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AcquisitionInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AcquisitionInfo* ptr = nullptr;
AcquisitionInfo* nullPointer = nullptr;
START_SECTION(AcquisitionInfo())
	ptr = new AcquisitionInfo();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~AcquisitionInfo())
	delete ptr;
END_SECTION

START_SECTION(const String& getMethodOfCombination() const)
  AcquisitionInfo tmp;
  TEST_EQUAL(tmp.getMethodOfCombination(),"");
END_SECTION

START_SECTION(void setMethodOfCombination(const String& method_of_combination))
  AcquisitionInfo tmp;
  tmp.setMethodOfCombination("TEST");
  TEST_EQUAL(tmp.getMethodOfCombination(),"TEST");
END_SECTION

START_SECTION(AcquisitionInfo(const AcquisitionInfo& source))
	AcquisitionInfo tmp;
	Acquisition a;
	a.setIdentifier("4711");
	tmp.push_back(a);
	tmp.setMethodOfCombination("Combo");
	tmp.setMetaValue("bla",4.0);

	AcquisitionInfo tmp2(tmp);
	TEST_EQUAL(tmp2.size(), 1);
	TEST_EQUAL(tmp2[0].getIdentifier(), "4711");  
	TEST_EQUAL(tmp2.getMethodOfCombination(), "Combo");  
	TEST_REAL_SIMILAR((double)(tmp2.getMetaValue("bla")), 4.0)
END_SECTION

START_SECTION(AcquisitionInfo& operator= (const AcquisitionInfo& source))
	AcquisitionInfo tmp;
	Acquisition a;
	a.setIdentifier("4711");
	tmp.push_back(a);
	tmp.setMethodOfCombination("Combo");
	tmp.setMetaValue("bla",4.0);

	//normal assignment
	AcquisitionInfo tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.size(), 1);
	TEST_EQUAL(tmp2[0].getIdentifier(), "4711");  
	TEST_EQUAL(tmp2.getMethodOfCombination(), "Combo");
	TEST_REAL_SIMILAR((double)(tmp2.getMetaValue("bla")), 4.0)
		
	//assignment of a empty value
	tmp2 = AcquisitionInfo();
	TEST_EQUAL(tmp2.size(), 0);
	TEST_EQUAL(tmp2.getMethodOfCombination(), "");
	TEST_EQUAL(tmp2.metaValueExists("bla"), false)
END_SECTION

START_SECTION(bool operator== (const AcquisitionInfo& rhs) const)
  AcquisitionInfo empty,edit;
	TEST_EQUAL(empty==edit,true);
	
	Acquisition a;
	edit.push_back(a);
	TEST_EQUAL(empty==edit,false);

	edit.setMetaValue("bla",4.0);
	TEST_EQUAL(empty==edit,false);

	edit = empty;
	edit.setMethodOfCombination("Combo");
	TEST_EQUAL(empty==edit,false);
END_SECTION

START_SECTION(bool operator!= (const AcquisitionInfo& rhs) const)
  AcquisitionInfo empty,edit;
	TEST_EQUAL(empty!=edit,false);
	
	Acquisition a;
	edit.push_back(a);
	TEST_EQUAL(empty!=edit,true);

	edit.setMetaValue("bla",4.0);
	TEST_EQUAL(empty!=edit,true);

	edit = empty;
	edit.setMethodOfCombination("Combo");
	TEST_EQUAL(empty!=edit,true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



