// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

START_SECTION((AASequence reversePeptides(const AASequence& protein, const String& protease)))
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

START_SECTION((std::vector<AASequence> shuffle(const AASequence& protein, const String& protease, int decoy_factor)))
  // Uses shufflePeptides with fixed seed (4711) for each peptide
  // Returns a vector with decoy_factor entries, each containing a complete shuffled protein
  
  // 1. no cutting site (full peptide is shuffled to minimize sequence identity)
  auto result1 = dg->shuffle(AASequence::fromString("TESTPEPTIDE"), "Trypsin");
  TEST_EQUAL(result1.size(), 1)
  TEST_EQUAL(result1[0].toString(),"DIESETEPTTP")
  
  // 2. cutting site after "TESTR", "RPEPTR" and "IDE" (each peptide is shuffled to minimize sequence identity)
  auto result2 = dg->shuffle(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P");
  TEST_EQUAL(result2.size(), 1)
  TEST_EQUAL(result2[0].toString(),"ETTSRTPRPEDIE")
  
  // 3. cutting site after "TESTRPEPTR", "IDE" (each peptide is shuffled to minimize sequence identity)
  auto result3 = dg->shuffle(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin");
  TEST_EQUAL(result3.size(), 1)
  TEST_EQUAL(result3[0].toString(), "ETEPTSRRTPDIE")
  
  // 4. decoy_factor=2 should generate 2 complete decoy proteins in the vector
  // Each variant gets a fresh DecoyGenerator with seed 4711, producing different shuffles
  auto result4 = dg->shuffle(AASequence::fromString("TESTPEPTIDE"), "Trypsin", 2);
  TEST_EQUAL(result4.size(), 2)
  TEST_EQUAL(result4[0].toString(), "DIESETEPTTP")
  TEST_EQUAL(result4[1].toString(), "PTEPIDEETTS")
  
  // 5. Top-down use-case: no protease cleavage (entire protein shuffled as one peptide)
  auto result5 = dg->shuffle(AASequence::fromString("LONGERPROTEINSEQUENCEFORTOPDOWN"), "no cleavage");
  TEST_EQUAL(result5.size(), 1)
  TEST_EQUAL(result5[0].toString(), "UPTECWLFREONNNORNOOEOESQEIDRGPT")
  
  // 6. Top-down with decoy_factor=2 (two different shuffled proteins in vector)
  auto result6 = dg->shuffle(AASequence::fromString("LONGERPROTEINSEQUENCEFORTOPDOWN"), "no cleavage", 2);
  TEST_EQUAL(result6.size(), 2)
  TEST_EQUAL(result6[0].toString(), "UPTECWLFREONNNORNOOEOESQEIDRGPT")
  TEST_EQUAL(result6[1].toString(), "FEIRPNEEDOROERNNCNETQGLOWOOPTUS")
END_SECTION

delete dg;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
