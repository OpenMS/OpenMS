// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#ifdef _OPENMP
#include <omp.h>
#endif

///////////////////////////
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedPeptideGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModifiedPeptideGenerator* ptr = nullptr;
ModifiedPeptideGenerator* null_ptr = nullptr;
START_SECTION(ModifiedPeptideGenerator())
{
  ptr = new ModifiedPeptideGenerator();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ModifiedPeptideGenerator())
{
  delete ptr;
}
END_SECTION

START_SECTION((static void applyFixedModifications(const ModifiedPeptideGenerator::MapToResidueType& fixed_mods,
                AASequence& peptide)))
{
  // query modification of interest from ModificationsDB
  StringList modNames;
  modNames << "Carbamidomethyl (C)";
  ModifiedPeptideGenerator::MapToResidueType fixed_mods = ModifiedPeptideGenerator::getModifications(modNames);

  AASequence seq0 = AASequence::fromString("AAAACAAAA"); // exactly one target site
  AASequence seq1 = AASequence::fromString("AAAAAAAAA"); // no target site
  AASequence seq2 = AASequence::fromString("CCCCCCCCC"); // all target sites
  AASequence seq3 = AASequence::fromString("AAAACAAC(Carbamidomethyl)AAA"); // one of two target sites already modified
  AASequence seq4 = AASequence::fromString("AAAACAAC(Oxidation)AAA"); // one of two target sites already modified

  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq0);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq1);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq2);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq3);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq4);

  TEST_EQUAL(seq0.toString(), "AAAAC(Carbamidomethyl)AAAA");
  TEST_EQUAL(seq1.toString(), "AAAAAAAAA");
  TEST_EQUAL(seq2.toString(), "C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)");
  TEST_EQUAL(seq3.toString(), "AAAAC(Carbamidomethyl)AAC(Carbamidomethyl)AAA");
  TEST_EQUAL(seq4.toString(), "AAAAC(Carbamidomethyl)AAC(Oxidation)AAA");

  // test terminal modifications
  modNames.clear();
  modNames << "Carbamyl (N-term)";

  fixed_mods.val.clear();
  fixed_mods = ModifiedPeptideGenerator::getModifications(modNames);

  seq0 = AASequence::fromString("KAAAAAAAA"); // exactly one target site
  seq1 = AASequence::fromString("K(Carbamyl)AAAAAAAA"); // ambiguous case: is mod Carbamyl (K) or (N-Term)?
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq0);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, seq1);
  TEST_EQUAL(seq0.toString(), ".(Carbamyl)KAAAAAAAA");
  TEST_EQUAL(seq1.toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAA");
 
}
END_SECTION

START_SECTION((static void applyVariableModifications(const ModifiedPeptideGenerator::MapToResidueType& var_mods, const AASequence& peptide, Size max_variable_mods_per_peptide, std::vector< AASequence > &all_modified_peptides, bool keep_unmodified=true)))
{
  // query modification of interest from ModificationsDB
  StringList modNames;
  modNames << "Oxidation (M)";

  ModifiedPeptideGenerator::MapToResidueType variable_mods = ModifiedPeptideGenerator::getModifications(modNames);

  vector<AASequence> modified_peptides;

  // test behavior if sequence empty
  AASequence seq;
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test behavior if peptide empty
  seq = AASequence::fromString("");
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test behavior if no target site in sequence
  seq = AASequence::fromString("AAAAAAAAA");
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test flag to preserve passed peptide
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides[0], AASequence::fromString("AAAAAAAAA")); // only the original peptide
  modified_peptides.clear();

  // test behavior if one target site in sequence and different number of maximum variable modifications are choosen
  seq = AASequence::fromString("AAAAMAAAA"); // exactly one target site
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // one generated peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // also only one generated peptide as there is only one target site
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  // test again keeping of the original peptides
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 2); // original and modified peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAMAAAA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();

  // test behavior if two target sites are in peptide and we need some combinatorics
  // only one additional variable modification
  seq = AASequence::fromString("AAMAAAMAA"); // two target site
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 2); // two modified peptides each modified at a different site
  TEST_EQUAL(modified_peptides[0].toString(), "AAMAAAM(Oxidation)AA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAM(Oxidation)AAAMAA");
  modified_peptides.clear();

  // up to two variable modifications per peptide
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 3); // three modified peptides, 1 modified at first M, 1 modified at second M and both M modified
  TEST_EQUAL(modified_peptides[0].toString(), "AAMAAAM(Oxidation)AA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAM(Oxidation)AAAMAA");
  TEST_EQUAL(modified_peptides[2].toString(), "AAM(Oxidation)AAAM(Oxidation)AA");
  modified_peptides.clear();

  // two different modifications with same target site
  modNames.clear();
  modNames << "Glutathione (C)" << "Carbamidomethyl (C)";
  variable_mods.val.clear();
  variable_mods = ModifiedPeptideGenerator::getModifications(modNames);

  seq = AASequence::fromString("ACAACAACA"); // three target sites
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, false);

  TEST_EQUAL(modified_peptides.size(), 6); // six modified peptides, as both C modifications can occur at all 3 sites
  std::sort(modified_peptides.begin(), modified_peptides.end(),
    [](const AASequence & a, const AASequence & b) -> bool
    { 
      return a.toString() < b.toString(); 
    });

  TEST_EQUAL(modified_peptides[0].toString(), "AC(Carbamidomethyl)AACAACA");
  TEST_EQUAL(modified_peptides[1].toString(), "AC(Glutathione)AACAACA");
  TEST_EQUAL(modified_peptides[2].toString(), "ACAAC(Carbamidomethyl)AACA");
  TEST_EQUAL(modified_peptides[3].toString(), "ACAAC(Glutathione)AACA");
  TEST_EQUAL(modified_peptides[4].toString(), "ACAACAAC(Carbamidomethyl)A");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 7); // same as above but +1 unmodified peptide

  modified_peptides.clear();

  seq = AASequence::fromString("ACAACAACA"); // three target sites and maximum of three occurrences of the two modifications Glutathione and Carbamidomethyl
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 3, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 3*3*3-1); // three sites with 3 possibilities (none, Glut., Carb.) each - but we need to subtract (none, none, none). The unmodified peptide

  modified_peptides.clear();
  seq = AASequence::fromString("AAAAC(Carbamidomethyl)AAAA"); // target site already occupied
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 3, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide as target site was already occupied

  // three different modifications
  modNames.clear();
  modNames << "Glutathione (C)" << "Carbamidomethyl (C)" << "Oxidation (M)";
  variable_mods.val.clear();
  variable_mods = ModifiedPeptideGenerator::getModifications(modNames);

  modified_peptides.clear();

  seq = AASequence::fromString("ACMACMACA"); // three target sites (C) and two target sites (M)
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 3, modified_peptides, false);

  // for exactly one mod: 2*3 + 2 = 8 (two possibilities each for each C sites and 1 for each of the two M sites)
  // for exactly two mods: 1 (iff two Ox) + 6 * 2 (iff one Ox at two pos.) + (3 choose 2) * 4 (iff no Ox) = 25
  // for exactly three mod: 6 (iff two Ox) + ((3 choose2) * 4) * 2 (iff one Ox at two pos.) + 8 (iff no Ox) = 38
  TEST_EQUAL(modified_peptides.size(), 8 + 25 + 38);

  // test terminal modifications
  modNames.clear();
  modNames << "Carbamyl (N-term)" << "Oxidation (M)";
  
  variable_mods.val.clear();
  variable_mods = ModifiedPeptideGenerator::getModifications(modNames);

  modified_peptides.clear();
  seq = AASequence::fromString("KAAAAAAAMA"); // exactly one target site
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 2, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 4);
  std::sort(modified_peptides.begin(), modified_peptides.end(),
    [](const AASequence & a, const AASequence & b) -> bool
    { 
      return a.toString() < b.toString(); 
    });

  TEST_EQUAL(modified_peptides[0].toString(), ".(Carbamyl)KAAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[1].toString(), ".(Carbamyl)KAAAAAAAMA");
  TEST_EQUAL(modified_peptides[2].toString(), "KAAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[3].toString(), "KAAAAAAAMA");
  
  modified_peptides.clear();
  seq = AASequence::fromString("K(Carbamyl)AAAAAAAMA"); 
  ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 2, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 4);
  std::sort(modified_peptides.begin(), modified_peptides.end(),
    [](const AASequence & a, const AASequence & b) -> bool
    { 
      return a.toString() < b.toString(); 
    });

  TEST_EQUAL(modified_peptides[0].toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[1].toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAMA");
  TEST_EQUAL(modified_peptides[2].toString(), "K(Carbamyl)AAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[3].toString(), "K(Carbamyl)AAAAAAAMA");
}
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{
// only do this in release, since MP errors are unlikely to occur in Debug mode anyway and it takes 5min to run the test in Debug
#ifdef NDEBUG
  int nr_iterations (1e5);
  int test = 0;
  // many modifications
  vector<String> all_mods = {"Carbamidomethyl (C)", "Oxidation (M)", "Carbamyl (M)", "Phospho (S)", "Phospho (T)", "Carbamyl (T)", "Phospho (Y)", "Carbamyl (K)", "Carbamyl (N-term)"};

  ModifiedPeptideGenerator::MapToResidueType variable_mods = ModifiedPeptideGenerator::getModifications(all_mods);

  const std::string s("ACDEFGHIKLMNPQRSTVWY");
  const AASequence seq = AASequence::fromString(s);

#pragma omp parallel for reduction (+: test)
  for (int i = 0; i < nr_iterations; ++i)
  {
    vector<AASequence> modified_peptides;
    modified_peptides.reserve(199);
    ModifiedPeptideGenerator::applyVariableModifications(variable_mods, seq, 4, modified_peptides, true);
    test += modified_peptides.size();
  }
  TEST_EQUAL(test, 199 * nr_iterations)
#endif
}
END_SECTION

START_SECTION([EXTRA] ModifiedPeptideGenerator::getModifications)
{
  StringList mods;
  mods << "Pyro-carbamidomethyl (N-term C)" << "Carbamidomethyl (C)" << "Carbamyl (N-term)" << "Deamidated (Protein N-term F)";
  auto mod_map = ModifiedPeptideGenerator::getModifications(mods);
  TEST_EQUAL(mods.size(), mod_map.val.size());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


