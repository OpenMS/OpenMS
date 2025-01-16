// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/openms_package_version.h>
#include <OpenMS/DATASTRUCTURES/String.h>


/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(VersionInfo, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(static String getTime())
{
  String t = VersionInfo::getTime();
  NOT_TESTABLE
}
END_SECTION

START_SECTION(static String getVersion() )
{
  TEST_STRING_EQUAL(VersionInfo::getVersion(), String(OPENMS_PACKAGE_VERSION).trim());
}
END_SECTION

START_SECTION((static VersionDetails getVersionStruct()))
{
  VersionInfo::VersionDetails detail;
  detail.version_major = 3;
  detail.version_minor = 4;
  detail.version_patch = 0;
  TEST_EQUAL(VersionInfo::getVersionStruct().version_major, detail.version_major);
  TEST_EQUAL(VersionInfo::getVersionStruct().version_minor, detail.version_minor);
  TEST_EQUAL(VersionInfo::getVersionStruct().version_patch, detail.version_patch);
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] VersionDetails()))
{
  VersionInfo::VersionDetails detail;
  TEST_EQUAL(detail == VersionInfo::VersionDetails::EMPTY, true)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] VersionDetails(const VersionDetails &other)))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c(detail);
  TEST_EQUAL(c.version_major, detail.version_major)
  TEST_EQUAL(c.version_minor, detail.version_minor)
  TEST_EQUAL(c.version_patch, detail.version_patch)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator<(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail < c, false)
  c.version_patch = 3;
  TEST_EQUAL(detail < c, true)
  c.version_patch = 1;
  TEST_EQUAL(detail < c, false)
  c.version_major = 2;
  TEST_EQUAL(detail < c, true)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator==(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_TRUE(detail == c)
  c.version_patch = 3;
  TEST_EQUAL(detail == c, false)
  c.version_patch = 1;
  TEST_EQUAL(detail == c, false)
  c.version_major = 2;
  TEST_EQUAL(detail == c, false)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator>(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail > c, false)
  c.version_patch = 3;
  TEST_EQUAL(detail > c, false)
  c.version_patch = 1;
  TEST_EQUAL(detail > c, true)
  c.version_patch = 11;
  TEST_EQUAL(detail < c, true)
  c.version_patch = 2;
  TEST_EQUAL(detail > c, false)
  
  // note that any version with a pre-release identifier should be "less than" the release version
  c.pre_release_identifier = "alpha";
  TEST_EQUAL(detail > c, true)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] static VersionDetails create(const String &version)))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_TRUE(detail == c)

  detail = VersionInfo::VersionDetails::create("1.9");
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 0;
  TEST_TRUE(detail == c)

  detail = VersionInfo::VersionDetails::create("1.0");
  c.version_major = 1;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_TRUE(detail == c)

  detail = VersionInfo::VersionDetails::create("somestring");
  c.version_major = 0;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_TRUE(detail == c)

  detail = VersionInfo::VersionDetails::create("1.2a.bla");
  c.version_major = 0;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_TRUE(detail == c)

  detail = VersionInfo::VersionDetails::create("1.2.1-bla");
  c.version_major = 1;
  c.version_minor = 2;
  c.version_patch = 1;
  c.pre_release_identifier = "bla";
  TEST_EQUAL(detail.version_major, c.version_major)
  TEST_EQUAL(detail.version_minor, c.version_minor)
  TEST_EQUAL(detail.version_patch, c.version_patch)
  TEST_EQUAL(detail.pre_release_identifier, c.pre_release_identifier)
  TEST_TRUE(detail == c)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
