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
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <limits>
#include <map>

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

START_SECTION((AASequence reversePeptides(const AASequence& protein, const std::string& protease)))
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"EDITPEPTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"TSETRTPEPREDI")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(),"TPEPRTSETREDI")
END_SECTION

START_SECTION((AASequence shufflePeptides(const AASequence& aas, const std::string& protease, const int max_atempts, int seed)))
  // 1. no cutting site (full peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"DIESETEPTTP")
  // 2. cutting site after "TESTR", "RPEPTR" and "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"ETTSRTPEPRIED")
  // 3. cutting site after "TESTRPEPTR", "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(), "SPERTETTPRIED")
END_SECTION

START_SECTION((std::vector<AASequence> shuffle(const AASequence& protein, const std::string& protease, int decoy_factor)))
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

START_SECTION((NASequence reversePeptides(const NASequence& rna, const std::string& rnase, bool keep_five_prime, bool keep_three_prime) const))
{
  const auto original = NASequence::fromString("pACGU[Um]ACGp");
  RNaseDigestion digestion;
  digestion.setEnzyme("RNase_T1");
  const auto reversed = dg->reversePeptides(original, "RNase_T1");
  TEST_EQUAL(reversed.toString(), "pCAGCA[Um]UGp")
  TEST_TRUE(digestion.getFragmentPositions(reversed) == digestion.getFragmentPositions(original))
  TEST_TRUE(reversed.getFormula() == original.getFormula())
  TEST_EQUAL(reversed.getFivePrimeMod(), original.getFivePrimeMod())
  TEST_EQUAL(reversed.getThreePrimeMod(), original.getThreePrimeMod())
  TEST_EXCEPTION(Exception::IllegalArgument, dg->reversePeptides(original, "unspecific cleavage"))
  TEST_EQUAL(dg->reversePeptides(NASequence::fromString("pAUGCp"), "no cleavage").toString(), "pCGUAp")
  // Also preserve before-cut and multi-nucleotide contexts, not just T1's 3' G.
  const auto contextual = NASequence::fromString("pACAUCGACAACGUACAp");
  for (const std::string rnase : {"RNase_MC1", "mazF", "RNase_4", "colicin_E5"})
  {
    digestion.setEnzyme(rnase);
    const auto transformed = dg->reversePeptides(contextual, rnase);
    TEST_TRUE(digestion.getFragmentPositions(transformed) == digestion.getFragmentPositions(contextual))
    TEST_TRUE(transformed.getFormula() == contextual.getFormula())
  }
}
END_SECTION

START_SECTION((NASequence shufflePeptides(const NASequence& rna, const std::string& rnase, int max_attempts, double identity_threshold, bool keep_five_prime, bool keep_three_prime)))
{
  const auto original = NASequence::fromString("pACUU[Um]ACGACUACGp");
  RNaseDigestion digestion;
  digestion.setEnzyme("RNase_T1");
  dg->setSeed(42);
  const auto shuffled = dg->shufflePeptides(original, "RNase_T1", 30, 0.5, true, true);
  const auto repeated = dg->shufflePeptides(original, "RNase_T1", 30, 0.5, true, true);
  TEST_EQUAL(shuffled, repeated)
  TEST_TRUE(shuffled.getFormula() == original.getFormula())
  TEST_REAL_SIMILAR(shuffled.getMonoWeight(), original.getMonoWeight())
  TEST_TRUE(digestion.getFragmentPositions(shuffled) == digestion.getFragmentPositions(original))
  TEST_EQUAL(shuffled[0], original[0])
  TEST_EQUAL(shuffled[shuffled.size() - 1], original[original.size() - 1])
  std::map<std::string, Size> before, after;
  for (Size i = 0; i < original.size(); ++i)
  {
    ++before[original[i]->getCode()];
    ++after[shuffled[i]->getCode()];
  }
  TEST_TRUE(before == after)
  const auto immutable = NASequence::fromString("p[Um][Um][Um]p");
  TEST_EQUAL(dg->shufflePeptides(immutable, "no cleavage", 1, 0.0), immutable)
  TEST_EQUAL(dg->shufflePeptides(original, "RNase_T1", 0, 0.0), original)
  TEST_EXCEPTION(Exception::InvalidParameter, dg->shufflePeptides(original, "RNase_T1", -1))
  TEST_EXCEPTION(Exception::InvalidParameter, dg->shufflePeptides(original, "RNase_T1", 1, std::numeric_limits<double>::quiet_NaN()))
}
END_SECTION

delete dg;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
