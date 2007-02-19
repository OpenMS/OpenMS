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
#include <OpenMS/FORMAT/IsotopeXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/IsotopeXMLHandler.h>

#include <vector>

///////////////////////////

START_TEST(IsotopeXMLFile, "$Id: IsotopeXMLFile_test.C 1300 2007-01-18 07:27:04Z martinlangwisch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IsotopeXMLFile* ptr;
IsotopeXMLFile xml_file;

CHECK((IsotopeXMLFile()))
	ptr = new IsotopeXMLFile();
RESULT

CHECK(void load(const String& filename, map< String, vector< pair< DoubleReal, DoubleReal > > >& isotope_informations) const throw (Exception::FileNotFound, Exception::ParseError))

	map< String, vector< pair< DoubleReal, DoubleReal > > > isotope_informations;
	xml_file.load("data/isotopes.xml", isotope_informations);
	
	TEST_EQUAL(isotope_informations["H"][0].first, 1.0078250321)
	TEST_EQUAL(isotope_informations["H"][0].second, 0.999885)
	TEST_EQUAL(isotope_informations["H"][1].first, 2.014101778)
	TEST_EQUAL(isotope_informations["H"][1].second, 0.000115)
	TEST_EQUAL(isotope_informations["Se"][0].first, 73.9224766)
	TEST_EQUAL(isotope_informations["Se"][0].second, 0.0089)
	TEST_EQUAL(isotope_informations["Se"][1].first, 75.9192141)
	TEST_EQUAL(isotope_informations["Se"][1].second, 0.0937)
	TEST_EQUAL(isotope_informations["Se"][2].first, 76.9199146)
	TEST_EQUAL(isotope_informations["Se"][2].second, 0.0763)
	TEST_EQUAL(isotope_informations["Se"][3].first, 77.917095)
	TEST_EQUAL(isotope_informations["Se"][3].second, 0.2377)
	TEST_EQUAL(isotope_informations["Se"][4].first, 79.9165218)
	TEST_EQUAL(isotope_informations["Se"][4].second, 0.4961)
	TEST_EQUAL(isotope_informations["Se"][5].first, 81.9167)
	TEST_EQUAL(isotope_informations["Se"][5].second, 0.0873)
RESULT


CHECK(void store(String filename, const map< String, vector< String> >& isotope_informations) const throw (Exception::UnableToCreateFile))

	map< String, vector< pair< DoubleReal, DoubleReal > > > isotope_informations;
	xml_file.load("data/isotopes.xml", isotope_informations);
	String temp_filename;
	NEW_TMP_FILE(temp_filename)
	xml_file.store(temp_filename, isotope_informations);
	isotope_informations.clear();
	xml_file.load(temp_filename, isotope_informations);
	/*
	TEST_EQUAL(isotope_informations["H"][0], 1.0078250321, 0.999885))
	TEST_EQUAL(isotope_informations["H"][1], 2.014101778, 0.000115))
	TEST_EQUAL(isotope_informations["Se"][0], 73.9224766, 0.0089))
	TEST_EQUAL(isotope_informations["Se"][1], 75.9192141, 0.0937))
	TEST_EQUAL(isotope_informations["Se"][2], 76.9199146, 0.0763))
	TEST_EQUAL(isotope_informations["Se"][3], 77.917095, 0.2377))
	TEST_EQUAL(isotope_informations["Se"][4], 79.9165218, 0.4961))
	TEST_EQUAL(isotope_informations["Se"][5], 81.9167, 0.0873))*/
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
