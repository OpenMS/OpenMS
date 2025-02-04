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
#include <OpenMS/METADATA/ScanWindow.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

// static_assert(OpenMS::Test::fulfills_rule_of_5<ScanWindow>(), "Must fulfill rule of 5");
// static_assert(OpenMS::Test::fulfills_rule_of_6<ScanWindow>(), "Must fulfill rule of 6");
// static_assert(OpenMS::Test::fulfills_fast_vector<ScanWindow>(), "Must have fast vector semantics");
// static_assert(std::is_nothrow_move_constructible<ScanWindow>::value, "Must have nothrow move constructible");

START_TEST(ScanWindow, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ScanWindow* ptr = nullptr;
ScanWindow* nullPointer = nullptr;
START_SECTION((ScanWindow()))
	ptr = new ScanWindow();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ScanWindow()))
	delete ptr;
END_SECTION

START_SECTION((ScanWindow(const ScanWindow& source)))
  ScanWindow tmp;
  tmp.begin = 1.0;
  tmp.end = 2.0;
  tmp.setMetaValue("label",String("label"));
  
  ScanWindow tmp2(tmp);
  TEST_REAL_SIMILAR(tmp2.begin, 1.0)
  TEST_REAL_SIMILAR(tmp2.end, 2.0)
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
END_SECTION

START_SECTION((ScanWindow& operator= (const ScanWindow& source)))
  ScanWindow tmp;
  tmp.begin = 1.0;
  tmp.end = 2.0;
  tmp.setMetaValue("label",String("label"));
  
  ScanWindow tmp2;
  tmp2 = tmp;
  TEST_REAL_SIMILAR(tmp2.begin, 1.0)
  TEST_REAL_SIMILAR(tmp2.end, 2.0)
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
END_SECTION

START_SECTION((bool operator==(const ScanWindow &source) const ))
  ScanWindow edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.begin = 1.0;
  TEST_EQUAL(edit==empty,false);
  
  edit = empty; 
  edit.end = 1.0;
  TEST_EQUAL(edit==empty,false);
  
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!=(const ScanWindow &source) const ))
  ScanWindow edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.begin = 1.0;
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty; 
  edit.end = 1.0;
  TEST_EQUAL(edit!=empty,true);
  
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



