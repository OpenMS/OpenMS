// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/DecoyGenerator.h>

///////////////////////////

START_TEST(DecoyGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DecoyGenerator* dg = nullptr;
DecoyGenerator* nullPointer = nullptr;

START_SECTION((DecoyGenerator()))
{
  dg = new DecoyGenerator();
  TEST_NOT_EQUAL(dg, nullPointer)
}
END_SECTION

START_SECTION((~DecoyGenerator()))
{
  delete dg;
}
END_SECTION

dg = new DecoyGenerator();
dg->setSeed(4711);

START_SECTION((AASequence reverseProtein(const AASequence& protein)))
  TEST_EQUAL(dg->reverseProtein(AASequence::fromString("PRTEINE")).toString(), "ENIETRP")
END_SECTION

START_SECTION((AASequence reversePeptide(const AASequence& protein, const String& protease)))
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"EDITPEPTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"TSETRTPEPREDI")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(),"TPEPRTSETREDI")
END_SECTION

START_SECTION((AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_atempts, int seed)))
  // 1. no cutting site (full peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"DIESETEPTTP")
  // 2. cutting site after "TESTR", "RPEPTR" and "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"ETTSRTPEPRIED")
  // 3. cutting site after "TESTRPEPTR", "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(), "SPERTETTPRIED")
END_SECTION

delete dg;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
