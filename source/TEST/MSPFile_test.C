// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPFile* ptr = 0;
CHECK((MSPFile()))
	ptr = new MSPFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MSPFile()))
	delete ptr;
RESULT

CHECK(template <typename MapType> void load(const String &filename, std::vector< PeptideIdentification > &ids, MapType &exp))
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
	PeakMap exp;
RESULT

CHECK((template <typename MapType> void store(const String &filename, const MapType &map) const throw (Exception::UnableToCreateFile)))

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
