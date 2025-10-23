// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/Types.h>
///////////////////////////

#include <clocale>
#include <string>

START_TEST(Types, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION([EXTRA] OpenMS_locale initialization preserves system locale)
{
  // This test verifies that initializing OpenMS doesn't permanently change
  // the system locale, which is important for Python bindings
  
  // Get the current locale
  const char* current = setlocale(LC_ALL, nullptr);
  TEST_NOT_EQUAL(current, nullptr);
  
  std::string current_locale(current);
  
  // OpenMS_locale should contain "C" 
  TEST_STRING_EQUAL(OpenMS::Internal::OpenMS_locale, "C");
  
  // But the system locale should still be what it was
  const char* after = setlocale(LC_ALL, nullptr);
  TEST_NOT_EQUAL(after, nullptr);
  
  std::string after_locale(after);
  
  // The system locale should not have been changed to "C" by OpenMS initialization
  // (unless it was "C" to begin with)
  TEST_EQUAL(current_locale, after_locale);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
