// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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
OMSSAXMLFile* nullPointer = nullptr;
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
	ModificationDefinitionsSet mod_set(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)"));
	ptr->setModificationDefinitionsSet(mod_set);
	NOT_TESTABLE
END_SECTION

START_SECTION(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, bool load_proteins=true, bool load_empty_hits = true))
  // two spectra, first with some hits (mapping to 4 proteins), second is empty
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications);

	TEST_EQUAL(protein_identification.getHits().size(), 4)
	TEST_EQUAL(peptide_identifications.size(), 2)

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications, false);
	TEST_EQUAL(protein_identification.getHits().size(), 0)
	TEST_EQUAL(peptide_identifications.size(), 2)

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications, false, false);
	TEST_EQUAL(protein_identification.getHits().size(), 0)
	TEST_EQUAL(peptide_identifications.size(), 1)

END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
