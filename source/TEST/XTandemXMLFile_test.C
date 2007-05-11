// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(XTandemXMLFile, "$Id: $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XTandemXMLFile xml_file;
XTandemXMLFile* ptr;
Identification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

CHECK((XTandemXMLFile()))
	ptr = new XTandemXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~XTandemXMLFile())
	delete ptr;
RESULT

ptr = new XTandemXMLFile();

CHECK(void load(const String& filename, Identification& protein_identification, std::vector<PeptideIdentification>& id_data) const throw(Exception::FileNotFound, Exception::ParseError))

	xml_file.load("data/XTandemXMLFile_test_1.xml",	protein_identification, peptide_identifications);

	// TODO
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
