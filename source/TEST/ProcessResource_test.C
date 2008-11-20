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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SYSTEM/ProcessResource.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProcessResource, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProcessResource* ptr = 0;
START_SECTION(ProcessResource())
{
	ptr = new ProcessResource();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ProcessResource())
{
	delete ptr;
}
END_SECTION

START_SECTION((static void LimitCPUTime(const Int &seconds)))
{
  // this is quite impossible to test, as on success, the programm will just terminate without throwing an exception.
	// There is a workaround for linux (involving Sig-Handlers), but none for Windows that I know of

	//but we can test if the function is callable and leave it at that
	ProcessResource::LimitCPUTime(19);
	NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



