// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, kg290 $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/test_config.h>

#include <algorithm>

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

START_SECTION((void startDeBruijn(Size k, UInt64 seed, const std::string& keep_residues); void addProteinToDeBruijn(const AASequence&); void finalizeDeBruijn(); AASequence deBruijn(const AASequence&) const))
  const auto protein_a = AASequence::fromString("APEPTIDEG");
  const auto protein_b = AASequence::fromString("KPEPTIDEA");
  DecoyGenerator debruijn;
  debruijn.startDeBruijn(2, 4711);
  debruijn.addProteinToDeBruijn(protein_a);
  debruijn.addProteinToDeBruijn(protein_b);
  debruijn.finalizeDeBruijn();

  const auto decoy_a = debruijn.deBruijn(protein_a).toUnmodifiedString();
  const auto decoy_b = debruijn.deBruijn(protein_b).toUnmodifiedString();
  TEST_STRING_EQUAL(decoy_a, "ETKEPEPIP")
  TEST_EQUAL(decoy_a.size(), protein_a.size())
  TEST_EQUAL(decoy_b.size(), protein_b.size())
  TEST_EQUAL(decoy_a.substr(3, 5), decoy_b.substr(3, 5)) // shared seven-residue repeat remains shared for L-k residues

  DecoyGenerator seeded_again;
  seeded_again.startDeBruijn(2, 4711);
  seeded_again.addProteinToDeBruijn(protein_a);
  seeded_again.addProteinToDeBruijn(protein_b);
  seeded_again.finalizeDeBruijn();
  TEST_EQUAL(seeded_again.deBruijn(protein_a).toUnmodifiedString(), decoy_a)

  DecoyGenerator keep_cleavage;
  keep_cleavage.startDeBruijn(2, 4711, "K");
  keep_cleavage.addProteinToDeBruijn(AASequence::fromString("MPEPTIDEK"));
  keep_cleavage.finalizeDeBruijn();
  TEST_EQUAL(keep_cleavage.deBruijn(AASequence::fromString("MPEPTIDEK")).toUnmodifiedString().back(), 'K')

  TEST_EXCEPTION(Exception::Precondition, debruijn.deBruijn(AASequence::fromString("WWWWWWW")))
  TEST_EXCEPTION(Exception::InvalidParameter, debruijn.deBruijn(AASequence::fromString("M(Oxidation)PEPTIDEK")))

  DecoyGenerator not_started;
  TEST_EXCEPTION(Exception::Precondition, not_started.addProteinToDeBruijn(protein_a))
  TEST_EXCEPTION(Exception::Precondition, not_started.finalizeDeBruijn())
  TEST_EXCEPTION(Exception::Precondition, not_started.deBruijn(protein_a))
  TEST_EXCEPTION(Exception::InvalidParameter, not_started.startDeBruijn(0, 4711))

  DecoyGenerator finalized;
  finalized.startDeBruijn(2, 4711);
  finalized.addProteinToDeBruijn(protein_a);
  finalized.finalizeDeBruijn();
  TEST_EXCEPTION(Exception::Precondition, finalized.addProteinToDeBruijn(protein_b))
  TEST_EXCEPTION(Exception::Precondition, finalized.finalizeDeBruijn())
END_SECTION

START_SECTION((AASequence reverseProtein(const AASequence& protein)))
  TEST_EQUAL(dg->reverseProtein(AASequence::fromString("PRTEINE")).toString(), "ENIETRP")
END_SECTION

START_SECTION((AASequence reversePeptides(const AASequence& protein, const std::string& protease)))
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"EDITPEPTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"TSETRTPEPREDI")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(),"TPEPRTSETREDI")
  // no enzymatic termini: the whole protein is reversed (its digest would be every subsequence of it)
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "unspecific cleavage").toString(),"EDIRTPEPRTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence::fromString("TESTRPEPTRIDE"), "no cleavage").toString(),"EDIRTPEPRTSET")
  TEST_EQUAL(dg->reversePeptides(AASequence(), "Trypsin").toString(), "")
END_SECTION

START_SECTION((AASequence shufflePeptides(const AASequence& aas, const std::string& protease, const int max_atempts, int seed)))
  // 1. no cutting site (full peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTPEPTIDE"), "Trypsin").toString(),"DIESETEPTTP")
  // 2. cutting site after "TESTR", "RPEPTR" and "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin/P").toString(),"ETTSRTPEPRIED")
  // 3. cutting site after "TESTRPEPTR", "IDE" (each peptide is shuffled to minimize sequence identity)
  TEST_EQUAL(dg->shufflePeptides(AASequence::fromString("TESTRPEPTRIDE"), "Trypsin").toString(), "SPERTETTPRIED")
  // 4. no enzymatic termini: the whole protein is shuffled as one sequence, as with 'no cleavage'
  //    (its digest would be every subsequence of it, and the decoy their concatenation)
  {
    DecoyGenerator dg_unspecific, dg_no_cleavage;
    dg_unspecific.setSeed(4711);
    dg_no_cleavage.setSeed(4711);
    const AASequence protein = AASequence::fromString("TESTRPEPTRIDE");
    const std::string decoy = dg_unspecific.shufflePeptides(protein, "unspecific cleavage").toString();
    TEST_EQUAL(decoy.size(), protein.size())
    TEST_EQUAL(decoy, dg_no_cleavage.shufflePeptides(protein, "no cleavage").toString())
    TEST_NOT_EQUAL(decoy, protein.toString())
    std::string sorted_decoy = decoy, sorted_target = protein.toString();
    std::sort(sorted_decoy.begin(), sorted_decoy.end());
    std::sort(sorted_target.begin(), sorted_target.end());
    TEST_EQUAL(sorted_decoy, sorted_target) // same residues
  }
  // 5. empty protein
  TEST_EQUAL(dg->shufflePeptides(AASequence(), "Trypsin").toString(), "")
  TEST_EQUAL(dg->shufflePeptides(AASequence(), "unspecific cleavage").toString(), "")
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

  // 7. unspecific cleavage: no enzymatic termini, the entire protein is shuffled as one peptide (as with 'no cleavage')
  auto result7 = dg->shuffle(AASequence::fromString("LONGERPROTEINSEQUENCEFORTOPDOWN"), "unspecific cleavage", 2);
  TEST_EQUAL(result7.size(), 2)
  TEST_EQUAL(result7[0].toString(), result6[0].toString())
  TEST_EQUAL(result7[1].toString(), result6[1].toString())

  // 8. empty protein
  auto result8 = dg->shuffle(AASequence(), "Trypsin");
  TEST_EQUAL(result8.size(), 1)
  TEST_EQUAL(result8[0].toString(), "")
END_SECTION

delete dg;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
