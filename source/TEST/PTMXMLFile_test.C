// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

#include <vector>

///////////////////////////

START_TEST(PTMXMLFile, "$Id: PTMXMLFile_test.C 1300 2007-01-18 07:27:04Z martinlangwisch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PTMXMLFile* ptr;
PTMXMLFile xml_file;

CHECK((PTMXMLFile()))
	ptr = new PTMXMLFile();
RESULT

CHECK(void load(const String& filename, map< String, pair< String, String > >& ptm_informations) const throw (Exception::FileNotFound, Exception::ParseError))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load("data/PTMs.xml", ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N(2) O(2) C(-1) H(-3)")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
RESULT


CHECK(void store(String filename, const map< String, vector< String> >& ptm_informations) const throw (Exception::UnableToCreateFile))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load("data/PTMs.xml", ptm_informations);
	String temp_filename;
	NEW_TMP_FILE(temp_filename)
	xml_file.store(temp_filename, ptm_informations);
	ptm_informations.clear();
	xml_file.load(temp_filename, ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N(2) O(2) C(-1) H(-3)")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
