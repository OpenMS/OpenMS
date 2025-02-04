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
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
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
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDScoreSwitcherAlgorithm_test_input.idXML"), protein_identifications, identifications, document_id);
PeptideIdentification identification = identifications[0];
ProteinIdentification protein_identification = protein_identifications[0];

IDScoreSwitcherAlgorithm* ptr = nullptr;
IDScoreSwitcherAlgorithm* null_ptr = nullptr;
START_SECTION(IDScoreSwitcherAlgorithm())
{
  ptr = new IDScoreSwitcherAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IDScoreSwitcherAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(switchToGeneralScoreType)
{
  IDScoreSwitcherAlgorithm switcher{};
  Size c(0);
  switcher.switchToGeneralScoreType(identifications, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
  TEST_EQUAL(identifications[0].getScoreType(), "Posterior Error Probability");
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
