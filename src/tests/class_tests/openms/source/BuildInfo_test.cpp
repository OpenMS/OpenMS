// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/SYSTEM/BuildInfo.h>

#include <set>
#include <string>

///////////////////////////

START_TEST(OpenMSOSInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

OpenMSOSInfo* ptr = nullptr;
OpenMSOSInfo* null_ptr = nullptr;

START_SECTION((OpenMSOSInfo()))
{
  ptr = new OpenMSOSInfo();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // a freshly constructed instance reports everything as "unknown"
  TEST_STRING_EQUAL(ptr->getOSAsString(), "unknown")
  TEST_STRING_EQUAL(ptr->getArchAsString(), "unknown")
  TEST_STRING_EQUAL(ptr->getOSVersionAsString(), "unknown")
}
END_SECTION

START_SECTION((~OpenMSOSInfo()))
{
  delete ptr;
}
END_SECTION

START_SECTION((std::string getOSAsString() const))
{
  TEST_STRING_EQUAL(OpenMSOSInfo().getOSAsString(), "unknown")
}
END_SECTION

START_SECTION((std::string getArchAsString() const))
{
  TEST_STRING_EQUAL(OpenMSOSInfo().getArchAsString(), "unknown")
}
END_SECTION

START_SECTION((std::string getOSVersionAsString() const))
{
  TEST_STRING_EQUAL(OpenMSOSInfo().getOSVersionAsString(), "unknown")
}
END_SECTION

START_SECTION((static std::string getBinaryArchitecture()))
{
  // derived from sizeof(size_t); platform-independent expectation
  std::string expected = (sizeof(size_t) == 8) ? "64 bit"
                       : (sizeof(size_t) == 4) ? "32 bit"
                       : "unknown";
  TEST_STRING_EQUAL(OpenMSOSInfo::getBinaryArchitecture(), expected)
}
END_SECTION

START_SECTION((static OpenMSOSInfo getOSInfo()))
{
  OpenMSOSInfo info = OpenMSOSInfo::getOSInfo();

  // the architecture is derived from the pointer width, so it must agree with
  // getBinaryArchitecture()
  TEST_STRING_EQUAL(info.getArchAsString(), OpenMSOSInfo::getBinaryArchitecture())

  // the running OS must be recognised on any supported build platform
  std::set<std::string> supported_os = {"Windows", "MacOS", "Linux"};
  TEST_EQUAL(supported_os.count(info.getOSAsString()) == 1, true)
  TEST_NOT_EQUAL(info.getOSAsString(), "unknown")

  // a version string is always reported (at worst the "unknown" fallback)
  TEST_EQUAL(info.getOSVersionAsString().empty(), false)
}
END_SECTION

START_SECTION((static std::string getActiveSIMDExtensions()))
{
  // build-dependent content, but a stable (deterministic) string per binary
  std::string simd = OpenMSOSInfo::getActiveSIMDExtensions();
  TEST_EQUAL(simd == OpenMSOSInfo::getActiveSIMDExtensions(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
