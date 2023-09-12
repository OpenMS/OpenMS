// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RNPxlModificationsGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RNPxlModificationsGenerator* ptr = nullptr;
RNPxlModificationsGenerator* null_ptr = nullptr;
START_SECTION(RNPxlModificationsGenerator())
{
	ptr = new RNPxlModificationsGenerator();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RNPxlModificationsGenerator())
{
	delete ptr;
}
END_SECTION

START_SECTION((static RNPxlModificationMassesResult initModificationMassesRNA(StringList target_nucleotides, StringList mappings, StringList restrictions, StringList modifications, String sequence_restriction, bool cysteine_adduct, Int max_length=4)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



