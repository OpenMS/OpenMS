// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Your Name $
// $Authors: Your Name $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HydrophobicityProfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HydrophobicityProfile* ptr = nullptr;
HydrophobicityProfile* null_ptr = nullptr;
START_SECTION(HydrophobicityProfile())
{
  ptr = new HydrophobicityProfile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~HydrophobicityProfile())
{
  delete ptr;
}
END_SECTION

START_SECTION(int testfunktion())
{
    TEST_EQUAL(ptr->testfunktion(),1);
}
END_SECTION

START_SECTION(double computeGRAVY(const AASequence& seq))
{
  AASequence seq = AASequence::fromString("ACDE");
  double gravy = ptr->computeGRAVY(seq);
  TEST_REAL_SIMILAR(gravy, -0.675); // Expected GRAVY value for "ACDE" with Kyte-Doolittle scale
  AASequence seq_2 = AASequence::fromString("");
  double gravy_2 = ptr->computeGRAVY(seq_2);
  TEST_REAL_SIMILAR(gravy_2, 0); // Expected GRAVY value for "ACDE" with Kyte-Doolittle scale
}
END_SECTION




// Add more test sections for each method...

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST