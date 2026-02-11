// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLPresets.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <unordered_set>

using namespace std;

///////////////////////////

#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace OpenMS;
using namespace std;

START_TEST(NuXLParameterParsing, "$Id$")

/////////////////////////////////////////////////////////////

NuXLParameterParsing* ptr = nullptr;
NuXLParameterParsing* nullPointer = nullptr;

START_SECTION(NuXLParameterParsing())
  ptr = new NuXLParameterParsing();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~NuXLParameterParsing())
  delete ptr;
END_SECTION

START_SECTION((static std::vector<ResidueModification> getModifications(StringList modNames)))
{
  StringList mods;
  mods.push_back("Phospho (S)");
  mods.push_back("Oxidation (M)");
  mods.push_back("Acetyl (N-term)");
  
  vector<ResidueModification> result = NuXLParameterParsing::getModifications(mods);
  TEST_EQUAL(result.size(), 3)
  TEST_EQUAL(result[0].getId(), "Phospho")
  TEST_EQUAL(result[1].getId(), "Oxidation")
  TEST_EQUAL(result[2].getId(), "Acetyl")
  TEST_EQUAL(result[2].getTermSpecificity(), ResidueModification::N_TERM)
}
END_SECTION

/*
START_SECTION((static NucleotideToFragmentAdductMap getTargetNucleotideToFragmentAdducts(StringList fragment_adducts)))
{
  StringList fragment_adducts;
  fragment_adducts.push_back("U:H-2O-1;U-H2O");  // Format is nucleotide:formula;name
  fragment_adducts.push_back("U:H-3P-1O-4;U-H3PO4");
  
  NuXLParameterParsing::NucleotideToFragmentAdductMap result =
    NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);
  
  TEST_EQUAL(result.size(), 1)  // Only U nucleotide
  TEST_EQUAL(result['U'].size(), 2)  // Two adducts for U
  
  // Test invalid format
  StringList invalid_adducts;
  invalid_adducts.push_back("Invalid format");
  NuXLParameterParsing::NucleotideToFragmentAdductMap empty_result =
    NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(invalid_adducts);
  TEST_EQUAL(empty_result.empty(), true)
}
END_SECTION

START_SECTION((static MS2AdductsOfSinglePrecursorAdduct getFeasibleFragmentAdducts(const String& exp_pc_adduct, const String& exp_pc_formula, const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts, const std::set<char>& can_xl, const bool always_add_default_marker_ions, const bool default_marker_ions_RNA)))
{
  // Create test data
  StringList fragment_adducts;
  fragment_adducts.push_back("U:H-2O-1;U-H2O");
  NuXLParameterParsing::NucleotideToFragmentAdductMap nuc_to_frag =
    NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);
  
  set<char> can_xl = {'U'};
  
  // Test monomer case
  MS2AdductsOfSinglePrecursorAdduct result =
    NuXLParameterParsing::getFeasibleFragmentAdducts(
      "U-H2O",  // precursor adduct
      "C9H11N2O7P",  // precursor formula
      nuc_to_frag,
      can_xl,
      true,
      true);
  
  TEST_EQUAL(result.feasible_adducts.size(), 1)  // One nucleotide
  TEST_NOT_EQUAL(result.marker_ions.size(), 0)  // Should have marker ions
  
  // Test empty precursor formula
  MS2AdductsOfSinglePrecursorAdduct empty_result =
    NuXLParameterParsing::getFeasibleFragmentAdducts(
      "U-H2O",
      "",
      nuc_to_frag,
      can_xl,
      true,
      true);
  
  TEST_EQUAL(empty_result.feasible_adducts.empty(), true)
  TEST_EQUAL(empty_result.marker_ions.empty(), true)
}
END_SECTION

START_SECTION((static std::vector<NuXLFragmentAdductDefinition> getMarkerIonsMassSet(const PrecursorsToMS2Adducts& pc2adducts)))
{
  // First create some test data
  StringList fragment_adducts;
  fragment_adducts.push_back("U:H-2O-1;U-H2O");
  NuXLParameterParsing::NucleotideToFragmentAdductMap nuc_to_frag =
    NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);
  
  set<char> can_xl = {'U'};
  
  MS2AdductsOfSinglePrecursorAdduct ms2_adducts =
    NuXLParameterParsing::getFeasibleFragmentAdducts(
      "U-H2O",
      "C9H11N2O7P",
      nuc_to_frag,
      can_xl,
      true,
      true);
  
  NuXLParameterParsing::PrecursorsToMS2Adducts pc2adducts;
  pc2adducts["U-H2O"] = ms2_adducts;
  
  vector<NuXLFragmentAdductDefinition> result =
    NuXLParameterParsing::getMarkerIonsMassSet(pc2adducts);
  
  TEST_NOT_EQUAL(result.size(), 0)  // Should have marker ions
  
  // Test empty input
  NuXLParameterParsing::PrecursorsToMS2Adducts empty_pc2adducts;
  vector<NuXLFragmentAdductDefinition> empty_result =
    NuXLParameterParsing::getMarkerIonsMassSet(empty_pc2adducts);
  
  TEST_EQUAL(empty_result.empty(), true)
}
END_SECTION
*/
START_SECTION((static PrecursorsToMS2Adducts getAllFeasibleFragmentAdducts(const NuXLModificationMassesResult& precursor_adducts, const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts, const std::set<char>& can_xl, const bool always_add_default_marker_ions, const bool default_marker_ions_RNA)))
{
  // Retrieve preset for RNA-UV (U). 
  // This preset contains all necessary information to generate the precursor adducts and fragment adducts.
  
  StringList modifications; // The precursor adducts are then used to generate all feasible fragment adducts.  
  StringList fragment_adducts; // these are responsible for shifted fragment ions. Their fragment adducts thus determine which shifts will be observed on b-,a-,y-ions
  String can_cross_link; // nucleotides that can directly cross-link
  // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
  StringList target_nucleotides;
  // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
  StringList mappings;

  NuXLPresets::getPresets("RNA-UV (U)", target_nucleotides, mappings, modifications, fragment_adducts, can_cross_link);

  // test preset strings
  TEST_EQUAL(target_nucleotides.size(), 4)
  StringList expected_target_nucleotides = {
    "A=C10H14N5O7P",
    "C=C9H14N3O8P",
    "G=C10H14N5O8P",
    "U=C9H13N2O9P"
  };
  TEST_TRUE(target_nucleotides == expected_target_nucleotides);
  if (target_nucleotides != expected_target_nucleotides)
  {
    for (const auto& t : target_nucleotides) { std::cerr << t << endl; }
  }

  // list all precursor adducts (including nucleotides that can't cross-link)
  StringList expected_modifications  = {
      "U:",
      "U:-H2O",
      "C:",
      "C:-H2O",
      "C:-NH3",
      "G:",
      "G:-H2O",
      "G:-NH3",
      "A:",
      "A:-NH3"
  };
  TEST_TRUE(modifications == expected_modifications);
  if (modifications != expected_modifications)
  {
    for (const auto& m : modifications) { std::cerr << m << endl; }
  }
  
  StringList expected_fragment_adducts = {
      "U:C3O;C3O",
      "U:C4H4N2O2;U'",
      "U:C4H2N2O1;U'-H2O",
      "U:C9H13N2O9P1;U",
      "U:C9H11N2O8P1;U-H2O",
      "U:C9H12N2O6;U-HPO3",
      "U:C9H10N2O5;U-H3PO4",
      "C:C4H5N3O;C'",
      "C:C4H3N3;C'-H2O",
      "C:C4H2N2O;C'-NH3",
      "C:C9H14N3O8P;C",
      "C:C9H11N2O8P;C-NH3",
      "C:C9H12N3O7P;C-H2O",
      "C:C9H9N2O7P;C-NH3-H2O",
      "C:C9H13N3O5;C-HPO3",
      "C:C9H11N3O4;C-H3PO4",
      "C:C9H10N2O5;C-NH3-HPO3",
      "C:C9H8N2O4;C-NH3-H3PO4",
      "G:C5H5N5O;G'",
      "G:C5H3N5;G'-H2O",
      "G:C5H2N4O;G'-NH3",
      "G:C10H14N5O8P;G",
      "G:C10H12N5O7P;G-H2O",
      "G:C10H11N4O8P;G-NH3",
      "G:C10H9N4O7P;G-NH3-H2O",
      "G:C10H13N5O5;G-HPO3",
      "G:C10H11N5O4;G-H3PO4",
      "G:C10H10N4O5;G-NH3-HPO3",
      "G:C10H8N4O4;G-NH3-H3PO4",
      "A:C5H5N5;A'",
      "A:C5H2N4;A'-NH3",
      "A:C10H14N5O7P;A",
      "A:C10H12N5O6P;A-H2O",
      "A:C10H11N4O7P;A-NH3",
      "A:C10H9N4O6P;A-NH3-H2O",
      "A:C10H13N5O4;A-HPO3",
      "A:C10H11N5O3;A-H3PO4",
      "A:C10H10N5O4;A-NH3-HPO3",
      "A:C10H8N5O3;A-NH3-H3PO4"
  };

  TEST_TRUE(fragment_adducts == expected_fragment_adducts);

  TEST_EQUAL(can_cross_link, "U")

  TEST_EQUAL(mappings.size(), 4)
  StringList expected_mapping = {
    "A->A",
    "C->C",
    "G->G",
    "U->U"
  };
  TEST_EQUAL(mappings.size(), 4)
  TEST_TRUE(mappings == expected_mapping);
  if (mappings != expected_mapping)
  {
    for (const auto& m : mappings) { std::cerr << m << endl; }
  }

  // convert string to set
  set<char> can_xl;
  for (const auto& c : can_cross_link) { can_xl.insert(c); } // sort and make unique

  NuXLModificationMassesResult  mm = NuXLModificationsGenerator::initModificationMassesNA(
            target_nucleotides,
            StringList(), 
            can_xl,
            mappings,
            modifications, 
            "", 
            false, 
            2);

  mm.formula2mass[""] = 0; // insert "null" modification otherwise peptides without NA will not be searched
  mm.mod_combinations[""].insert("none");

  // check content of formula2mass and ensure order is correct
  StringList formula = {
      "",
      "C18H22N4O16P2",
      "C18H23N5O15P2",
      "C18H24N4O17P2",
      "C18H25N5O16P2",
      "C19H22N6O15P2",
      "C19H22N6O16P2",
      "C19H23N7O14P2",
      "C19H23N7O15P2",
      "C19H25N7O15P2",
      "C19H25N7O16P2",
      "C9H11N2O8P1",
      "C9H13N2O9P1"
  };

  DoubleList mass = {
      0,
      612.051,
      611.067,
      630.061,
      629.077,
      636.062,
      652.057,
      635.078,
      651.073,
      653.088,
      669.083,
      306.025,
      324.036
  };

  // compare using approximate double comparison due to floating point precision
  for (Size i = 0; i != formula.size(); ++i)
  {
    TEST_REAL_SIMILAR(mm.formula2mass[formula[i]], mass[i])
  }

  std::vector<std::string> expected_precursor_formula = {
      "",
      "C18H22N4O16P2",
      "C18H23N5O15P2",
      "C18H24N4O17P2",
      "C18H25N5O16P2",
      "C19H22N6O15P2",
      "C19H22N6O16P2",
      "C19H23N7O14P2",
      "C19H23N7O15P2",
      "C19H25N7O15P2",
      "C19H25N7O16P2",
      "C9H11N2O8P1",
      "C9H13N2O9P1"
  };

  std::vector<std::vector<std::string>> expected_precursors = {
      {"none"},
      {"UC-H3N1", "UU-H2O1"},
      {"CU-H2O1"},
      {"UU"},
      {"CU"},
      {"UA-H3N1"},
      {"UG-H3N1"},
      {"AU-H2O1"},
      {"GU-H2O1"},
      {"AU"},
      {"GU"},
      {"U-H2O1"},
      {"U"}
  };

  // check content of mod_combinations (order matters)
  for (Size i = 0; i != expected_precursor_formula.size(); ++i)
  {
    TEST_EQUAL(mm.mod_combinations[expected_precursor_formula[i]].size(), expected_precursors[i].size())
    for (Size j = 0; j != expected_precursors[i].size(); ++j)
    {
      auto it = mm.mod_combinations[expected_precursor_formula[i]].begin();
      std::advance(it, j);
      TEST_EQUAL(*it, expected_precursors[i][j])
    }
  }

  // first, we determine which fragments adducts can be generated from a single nucleotide (that has no losses)
  NuXLParameterParsing::NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts = NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);

  // check if all nucleotides are covered by the map (required to generate fragment adducts)

  std::vector<std::pair<char, std::vector<std::string>>> expected_nucleotides_to_fragments = {
      {'A', {
          "A'-NH3",
          "A'",
          "A-NH3-H3PO4",
          "A-H3PO4",
          "A-NH3-HPO3",
          "A-HPO3",
          "A-NH3-H2O",
          "A-H2O",
          "A-NH3",
          "A"
      }},
      {'C', {
          "C'-H2O",
          "C'-NH3",
          "C'",
          "C-NH3-H3PO4",
          "C-H3PO4",
          "C-NH3-HPO3",
          "C-HPO3",
          "C-NH3-H2O",
          "C-H2O",
          "C-NH3",
          "C"
      }},
      {'G', {
          "G'-H2O",
          "G'-NH3",
          "G'",
          "G-NH3-H3PO4",
          "G-H3PO4",
          "G-NH3-HPO3",
          "G-HPO3",
          "G-NH3-H2O",
          "G-H2O",
          "G-NH3",
          "G"
      }},
      {'U', {
          "C3O",
          "U'-H2O",
          "U'",
          "U-H3PO4",
          "U-HPO3",
          "U-H2O",
          "U"
      }}
  };

  TEST_EQUAL(nucleotide_to_fragment_adducts.size(), expected_nucleotides_to_fragments.size())
  
  for (size_t nt = 0; nt != expected_nucleotides_to_fragments.size(); ++nt)
  {
    // test target nucleotide formula
    auto it = nucleotide_to_fragment_adducts.begin();
    std::advance(it, nt);
    TEST_EQUAL(expected_nucleotides_to_fragments[nt].first, it->first); // order should match the expected target nucleotides
    TEST_EQUAL(expected_nucleotides_to_fragments[nt].second.size(), it->second.size())
    if (expected_nucleotides_to_fragments[nt].second.size() != it->second.size())
    {
      for (const auto& f : it->second) { std::cerr << f.name << endl; }
    }
    for (size_t f = 0; f != expected_nucleotides_to_fragments[nt].second.size(); ++f)
    {
      auto fit = it->second.begin();
      std::advance(fit, f);
      // order should match exactly the expected fragments
      TEST_TRUE(expected_nucleotides_to_fragments[nt].second[f] == fit->name)
    }
  }


  // calculate all feasible fragment adducts from all possible precursor adducts
  NuXLParameterParsing::PrecursorsToMS2Adducts all_feasible_fragment_adducts = NuXLParameterParsing::getAllFeasibleFragmentAdducts(mm, nucleotide_to_fragment_adducts, can_xl, true, true);

  // print all chemically feasible fragment adducts
  std::vector<std::pair<std::string, std::pair<char, std::vector<std::string>>>> expected_pc2nuc2fragment_adducts = {
      {"AU", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"AU-H2O1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"CU", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"CU-H2O1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"GU", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"GU-H2O1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"U", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"U-H2O1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-H2O"}}},
      {"UA-H3N1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"UC-H3N1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"UG-H3N1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"UU", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"UU-H2O1", {'U', {"C3O", "U'-H2O", "U'", "U-H3PO4", "U-HPO3", "U-H2O", "U"}}},
      {"none", {'U', {}}}
  };

  TEST_EQUAL(all_feasible_fragment_adducts.size(), expected_pc2nuc2fragment_adducts.size())

  // Compare each precursor and its associated fragments in order
  auto it_expected = expected_pc2nuc2fragment_adducts.begin();
  auto it_actual = all_feasible_fragment_adducts.begin();

  while (it_expected != expected_pc2nuc2fragment_adducts.end() && it_actual != all_feasible_fragment_adducts.end())
  {
    const String& expected_precursor = it_expected->first;
    const String& actual_precursor = it_actual->first;

    // Ensure the order of precursors is preserved
    TEST_EQUAL(expected_precursor, actual_precursor)

    const auto& expected_adducts = it_expected->second;
    const MS2AdductsOfSinglePrecursorAdduct& ms2adducts = it_actual->second;
    const std::vector<NucleotideToFeasibleFragmentAdducts>& nt2feasible = ms2adducts.feasible_adducts;

    if (!nt2feasible.empty())
    {
      const auto& n2fsa = nt2feasible[0];

      // Check nucleotide
      TEST_EQUAL(n2fsa.first, expected_adducts.first)

      // Check fragment adducts
      TEST_EQUAL(n2fsa.second.size(), expected_adducts.second.size())
      for (Size i = 0; i < expected_adducts.second.size(); ++i)
      {
        if (i < n2fsa.second.size())
        {
          TEST_EQUAL(n2fsa.second[i].name, expected_adducts.second[i])
        }
      }
    }

    ++it_expected;
    ++it_actual;
  }
}
END_SECTION

START_SECTION(([EXTRA] std::hash<NuXLFragmentAdductDefinition>))
{
  // Test that equal definitions have equal hashes
  NuXLFragmentAdductDefinition fad1(EmpiricalFormula("C4H4N2O2"), "U'", 112.027);
  NuXLFragmentAdductDefinition fad2(EmpiricalFormula("C4H4N2O2"), "U'", 112.027);

  std::hash<NuXLFragmentAdductDefinition> hasher;
  TEST_EQUAL(hasher(fad1), hasher(fad2))

  // Test that different formulas produce different hashes
  NuXLFragmentAdductDefinition fad3(EmpiricalFormula("C4H2N2O1"), "U'-H2O", 94.016);
  TEST_NOT_EQUAL(hasher(fad1), hasher(fad3))

  // Test that different names produce different hashes (same formula)
  NuXLFragmentAdductDefinition fad4(EmpiricalFormula("C4H4N2O2"), "different_name", 112.027);
  TEST_NOT_EQUAL(hasher(fad1), hasher(fad4))

  // Note: mass is not used in operator== so it should not affect the hash
  NuXLFragmentAdductDefinition fad5(EmpiricalFormula("C4H4N2O2"), "U'", 999.999);
  TEST_EQUAL(hasher(fad1), hasher(fad5))

  // Test empty definitions
  NuXLFragmentAdductDefinition fad_empty1;
  NuXLFragmentAdductDefinition fad_empty2;
  TEST_EQUAL(hasher(fad_empty1), hasher(fad_empty2))

  // Test that definitions work in unordered_set
  std::unordered_set<NuXLFragmentAdductDefinition> fad_set;
  fad_set.insert(fad1);
  fad_set.insert(fad2);  // Duplicate, should not increase size
  fad_set.insert(fad3);  // Different definition
  TEST_EQUAL(fad_set.size(), 2)
  TEST_EQUAL(fad_set.count(fad1), 1)
  TEST_EQUAL(fad_set.count(fad2), 1)  // Same as fad1
  TEST_EQUAL(fad_set.count(fad3), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST