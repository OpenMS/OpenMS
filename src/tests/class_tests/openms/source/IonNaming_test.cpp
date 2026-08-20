// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/IonNaming.h>

#include <climits>
#include <string>

///////////////////////////

START_TEST(IonNaming, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((std::string chargeSuffix(int charge)))
{
  // unknown charge is never spelled out
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(0), "")

  // small charges as a run of signs - what TheoreticalSpectrumGenerator has always written
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(1), "+")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(2), "++")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(5), "+++++")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(8), "++++++++")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(-1), "-")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(-3), "---")

  // beyond MAX_REPEATED_SIGNS the number is written out, which bounds the allocation for a
  // charge that came from a file
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(9), "+9")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(-12), "-12")
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(2000000000), "+2000000000")
  // INT_MIN must not be negated in int (undefined behaviour) and must not size an allocation
  TEST_STRING_EQUAL(IonNaming::chargeSuffix(INT_MIN), "-2147483648")
  TEST_EQUAL(IonNaming::chargeSuffix(INT_MIN).size() < 20, true)
}
END_SECTION

START_SECTION((int chargeFromName(const std::string& ion_name)))
{
  // no charge spelled out -> 0, which is distinguishable from charge 1
  TEST_EQUAL(IonNaming::chargeFromName(""), 0)
  TEST_EQUAL(IonNaming::chargeFromName("y5"), 0)
  TEST_EQUAL(IonNaming::chargeFromName("[alpha|ci$y3]"), 0)
  TEST_EQUAL(IonNaming::chargeFromName("a2-B"), 0)
  TEST_EQUAL(IonNaming::chargeFromName("my note"), 0)
  // a neutral loss ending in a digit is not a charge
  TEST_EQUAL(IonNaming::chargeFromName("y1-H2O1"), 0)
  TEST_EQUAL(IonNaming::chargeFromName("b3-H3N1"), 0)

  // run of signs
  TEST_EQUAL(IonNaming::chargeFromName("y5+"), 1)
  TEST_EQUAL(IonNaming::chargeFromName("b3++"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y1-H2O1+"), 1)
  TEST_EQUAL(IonNaming::chargeFromName("[alpha|ci$y3]++"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("c1-"), -1)
  TEST_EQUAL(IonNaming::chargeFromName("a3-B-"), -1)
  TEST_EQUAL(IonNaming::chargeFromName("w12--"), -2)

  // sign plus number
  TEST_EQUAL(IonNaming::chargeFromName("y5+3"), 3)
  TEST_EQUAL(IonNaming::chargeFromName("y5-3"), -3)

  // mzPAF caret, alone and followed by further mzPAF fields
  TEST_EQUAL(IonNaming::chargeFromName("b2^2"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y3-H2O^2"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y4-H2O1^2/3.2ppm"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y3^2*0.75"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("c1^-1"), -1)
  // mzPAF's mass delta is itself a signed number at the end of the name; the caret must win,
  // otherwise a charge 2 ion is read as charge -1
  TEST_EQUAL(IonNaming::chargeFromName("y4-H2O1^2/-1"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y4^2/-1"), 2)
  // a bare caret without a number is not a charge
  TEST_EQUAL(IonNaming::chargeFromName("note ^ see docs"), 0)

  // only the ion name in the first line is inspected, further lines are free text
  TEST_EQUAL(IonNaming::chargeFromName("y5++\nsome user comment"), 2)
  TEST_EQUAL(IonNaming::chargeFromName("y5\nsome user comment"), 0)

  // absurdly long digit runs are not charges and must not overflow
  TEST_EQUAL(IonNaming::chargeFromName("y5+99999999999999999999"), 0)
  TEST_EQUAL(IonNaming::chargeFromName("y5^99999999999999999999"), 0)
}
END_SECTION

START_SECTION(([EXTRA] chargeSuffix and chargeFromName round trip))
{
  // what a producer writes is what a consumer reads back, for every charge a spectrum can carry
  for (int z = -20; z <= 20; ++z)
  {
    if (z == 0) { continue; }
    const std::string name = "y5" + IonNaming::chargeSuffix(z);
    TEST_EQUAL(IonNaming::chargeFromName(name), z)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
