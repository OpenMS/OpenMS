// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Chris Bielow $
// $Authors: Marc Sturm, Clemens Groepl $
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
	NOT_TESTABLE
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

START_SECTION(static String getVersion() )
{
	TEST_STRING_EQUAL(VersionInfo::getVersion(),String(OPENMS_PACKAGE_VERSION).trim());
}
END_SECTION

START_SECTION(static Int getMajorVersion())
{
	STATUS("We might need to update this for a new release, oops!");
	TEST_EQUAL(VersionInfo::getMajorVersion(), 1);
}
END_SECTION

START_SECTION(static Int getMinorVersion())
{
	STATUS("We might need to update this for a new release, oops!");
	TEST_EQUAL(VersionInfo::getMinorVersion(), 7);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
