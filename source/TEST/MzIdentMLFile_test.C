// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

START_TEST(MzIdentMLFile, "$Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MzIdentMLFile* ptr = 0;
MzIdentMLFile* nullPointer = 0;
START_SECTION((MzIdentMLFile()))
	ptr = new MzIdentMLFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzIdentMLFile()))
	delete ptr;
END_SECTION

START_SECTION((template <typename MapType> void load(const String& filename, Identification& id)))
	MzIdentMLFile file;
	Identification id;
	file.load(OPENMS_GET_TEST_DATA_PATH("Mascot_MSMS_example.mzid"), id);
	// TODO
END_SECTION

START_SECTION((template <typename MapType> void store(const String& filename, const Identification& id) const))
	MzIdentMLFile file;

	// TODO
END_SECTION

START_SECTION(bool isValid(const String& filename, std::ostream& os = std::cerr))
  MzIdentMLFile file;
	// TODO
END_SECTION

START_SECTION(bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
	// TODO
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

