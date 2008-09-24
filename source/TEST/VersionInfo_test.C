// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>


///////////////////////////

START_TEST(VersionInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CHECK(static String getVersionAndTime())
	cout << "'" << VersionInfo::getVersionAndTime() << "'" << endl;
	cout << "'" << PACKAGE_VERSION << "'" << endl;
	TEST_EQUAL(VersionInfo::getVersionAndTime().hasPrefix(String(PACKAGE_VERSION).trim()),true)
RESULT

CHECK(static String getRevision() )
{
	// just to call the method
	STATUS("This should print a number if you have compiled from an SVN sandbox.  Compare with \"svnversion\" or \"svn info\".")
	STATUS(VersionInfo::getRevision());
	NOT_TESTABLE;
}
RESULT

CHECK(static String getVersion() )
{
	STATUS("We will need to update this for a new release, oops!");
	TEST_STRING_EQUAL(VersionInfo::getVersion(),"1.2");
}
RESULT

CHECK(static Int getMajorVersion())
{
	STATUS("We might need to update this for a new release, oops!");
	TEST_EQUAL(VersionInfo::getMajorVersion(), 1);
}
RESULT

CHECK(static Int getMinorVersion())
{
	STATUS("We might need to update this for a new release, oops!");
	TEST_EQUAL(VersionInfo::getMinorVersion(), 2);
}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
