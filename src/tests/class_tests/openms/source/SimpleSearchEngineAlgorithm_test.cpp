// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SimpleSearchEngineAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimpleSearchEngineAlgorithm* ptr = 0;
SimpleSearchEngineAlgorithm* null_ptr = 0;
START_SECTION(SimpleSearchEngineAlgorithm())
{
	ptr = new SimpleSearchEngineAlgorithm();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~SimpleSearchEngineAlgorithm())
{
	delete ptr;
}
END_SECTION

START_SECTION((ExitCodes search(const String &in_mzML, const String &in_db, std::vector< ProteinIdentification > &prot_ids, std::vector< PeptideIdentification > &pep_ids) const ))
{
  // tested via tool
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



