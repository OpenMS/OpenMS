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

#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(XTandemXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XTandemXMLFile* ptr;
XTandemXMLFile* nullPointer = nullptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications;

START_SECTION((XTandemXMLFile()))
	ptr = new XTandemXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~XTandemXMLFile())
	delete ptr;
END_SECTION

XTandemXMLFile xml_file;

START_SECTION(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, ModificationDefinitionsSet& mod_def_set))
{
	ModificationDefinitionsSet mod_set(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)"));

	xml_file.load(OPENMS_GET_TEST_DATA_PATH("XTandemXMLFile_test.xml"), protein_identification, peptide_identifications, mod_set);
	TEST_EQUAL(peptide_identifications.size(), 303);
	TEST_EQUAL(protein_identification.getHits().size(), 497);
  // should have picked up the default N-terminal modifications:
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 6);
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 0);

  mod_set.setModifications("", "Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)");
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("XTandemXMLFile_test_2.xml"), protein_identification, peptide_identifications, mod_set);
	TEST_EQUAL(peptide_identifications.size(), 2);
	TEST_EQUAL(protein_identification.getHits().size(), 21);
  // no additional modifications in this case:
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 3);
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
