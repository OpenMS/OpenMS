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
#include <OpenMS/FORMAT/OMSSACSVFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(OMSSACSVFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

OMSSACSVFile xml_file;
OMSSACSVFile* ptr;
OMSSACSVFile* nullPointer = nullptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

START_SECTION((OMSSACSVFile()))
	ptr = new OMSSACSVFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~OMSSACSVFile())
	delete ptr;
END_SECTION

ptr = new OMSSACSVFile();

START_SECTION(void load(const String &filename, ProteinIdentification &protein_identification, std::vector< PeptideIdentification > &id_data) const)
	ptr->load(OPENMS_GET_TEST_DATA_PATH("OMSSACSVFile_test_1.csv"), protein_identification, peptide_identifications);
	TEST_EQUAL(protein_identification.getHits().size(), 0)
	TEST_EQUAL(peptide_identifications.size(), 1)
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
