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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/SYSTEM/FileWatcher.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ControlledVocabulary, "$Id: FileWatcher_test.C 6114 2009-10-13 21:43:07Z marc_sturm $")

/////////////////////////////////////////////////////////////

FileWatcher* ptr = 0;
START_SECTION(FileWatcher(QObject *parent=0))
	ptr = new FileWatcher();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~FileWatcher())
	delete ptr;
END_SECTION

START_SECTION(void setDelayInSeconds(DoubleReal delay))
	NOT_TESTABLE
END_SECTION

START_SECTION(void addFile(const String& path))
	NOT_TESTABLE
END_SECTION

START_SECTION(void removeFile(const String& path))
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
