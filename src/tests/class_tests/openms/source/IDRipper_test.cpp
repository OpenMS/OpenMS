// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IDRipper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


///load input data
std::vector< ProteinIdentification > protein_identifications;
std::vector< PeptideIdentification > identifications;
String document_id;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test.idXML"), protein_identifications, identifications, document_id);
PeptideIdentification identification = identifications[0];
ProteinIdentification protein_identification = protein_identifications[0];

IDRipper* ptr = nullptr;
IDRipper* null_ptr = nullptr;
START_SECTION(IDRipper())
{
  ptr = new IDRipper();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IDRipper())
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual ~IDRipper()))
{
  // TODO
}
END_SECTION

START_SECTION((void rip(std::map< String, std::pair< std::vector< ProteinIdentification >, std::vector< PeptideIdentification > > > &ripped, std::vector< ProteinIdentification > &proteins, std::vector< PeptideIdentification > &peptides)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
