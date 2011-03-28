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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(OMSSAXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

OMSSAXMLFile xml_file;
OMSSAXMLFile* ptr;
OMSSAXMLFile* nullPointer = 0;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

START_SECTION((OMSSAXMLFile()))
	ptr = new OMSSAXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~OMSSAXMLFile())
	delete ptr;
END_SECTION

ptr = new OMSSAXMLFile();

START_SECTION(void setModificationDefinitionsSet(const ModificationDefinitionsSet &rhs))
	ModificationDefinitionsSet mod_set("", "Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)");
	ptr->setModificationDefinitionsSet(mod_set);
	NOT_TESTABLE
END_SECTION

START_SECTION(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, bool load_proteins=true))

	xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications);
	OMSSAXMLFile xml_file;

	TEST_EQUAL(protein_identification.getHits().size(), 4)
	TEST_EQUAL(peptide_identifications.size(), 1)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
