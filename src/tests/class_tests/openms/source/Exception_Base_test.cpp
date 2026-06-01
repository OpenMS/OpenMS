// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

using namespace OpenMS;

START_TEST(Exception::Base, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Exception::BaseException* e_ptr = nullptr;
Exception::BaseException* e_nullPointer = nullptr;
START_SECTION(Base() )
	e_ptr = new Exception::BaseException;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~Base() )
	delete e_ptr;
END_SECTION


START_SECTION(Base(const Base& exception) )
  // ???
END_SECTION

START_SECTION((Base(const char* file, int line, const char* function) ))
  // ???
END_SECTION

START_SECTION((Base(const char* file, int line, const char* function, const std::string& name, const std::string& message) ))
  // ???
END_SECTION

START_SECTION(const char* getFile() const )
  // ???
END_SECTION

START_SECTION(const char* getName() const )
  // ???
END_SECTION

START_SECTION(const char* what() const )
  // ???
END_SECTION

START_SECTION(int getLine() const )
  // ???
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
