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
// $Maintainer: Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

#include <vector>

///////////////////////////

START_TEST(PTMXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PTMXMLFile* ptr = 0;
PTMXMLFile xml_file;

START_SECTION((PTMXMLFile()))
	ptr = new PTMXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((void load(const String& filename, std::map< String, std::pair< String, String > >& ptm_informations)))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION


START_SECTION((void store(String filename, std::map< String, std::pair< String, String > > &ptm_informations) const))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	String temp_filename;
	NEW_TMP_FILE(temp_filename)
	xml_file.store(temp_filename, ptm_informations);
	ptm_informations.clear();
	xml_file.load(temp_filename, ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
