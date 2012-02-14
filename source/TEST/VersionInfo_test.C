// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors:  Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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
	TEST_STRING_EQUAL(VersionInfo::getVersion(),String(OPENMS_PACKAGE_VERSION).trim());
}
END_SECTION


START_SECTION((static VersionDetails getVersionStruct()))
{
	VersionInfo::VersionDetails detail;
	detail.version_major = 1;
  detail.version_minor = 10;
	TEST_EQUAL(VersionInfo::getVersionStruct() == detail, true);
}
END_SECTION

START_SECTION(static String getRevision() )
{
	// just to call the method
	STATUS("If you have compiled from an SVN sandbox, then this should print a revision number, or a range of revisions followed by \"M\", or something similar.");
  STATUS("Compare with \"svnversion\" or \"svn info\".")
	STATUS(VersionInfo::getRevision());
	NOT_TESTABLE;
}
END_SECTION


START_SECTION(([VersionInfo::VersionDetails] VersionDetails()))
{
  VersionInfo::VersionDetails detail;
	TEST_EQUAL(detail == VersionInfo::VersionDetails::EMPTY, true)
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
	TEST_EQUAL(detail == c, true)
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
	c.version_major = 2;
	TEST_EQUAL(detail > c, false)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] static VersionDetails create(const String &version)))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
	VersionInfo::VersionDetails c;
	c.version_major = 1;
	c.version_minor = 9;
	c.version_patch = 2;
	TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.9");
	c.version_major = 1;
	c.version_minor = 9;
	c.version_patch = 0;
	TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.0");
	c.version_major = 1;
	c.version_minor = 0;
	c.version_patch = 0;
	TEST_EQUAL(detail == c, true)
	
  detail = VersionInfo::VersionDetails::create("somestring");
	c.version_major = 0;
	c.version_minor = 0;
	c.version_patch = 0;
	TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.2a.bla");
	c.version_major = 0;
	c.version_minor = 0;
	c.version_patch = 0;
	TEST_EQUAL(detail == c, true)

	
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
