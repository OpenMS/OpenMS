// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/CometModification.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace OpenMS;
using namespace std;

START_TEST(CometModification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CometModification* ptr = nullptr;
CometModification* null_ptr = nullptr;

START_SECTION(CometModification())
{
  ptr = new CometModification();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~CometModification())
{
  delete ptr;
}
END_SECTION

START_SECTION(CometModification(const ResidueModification* mod, int binary_grp, int max_mods))
{
  // Test construction from a regular amino acid modification (Oxidation on M)
  const ResidueModification* ox = ModificationsDB::getInstance()->getModification("Oxidation", "M", ResidueModification::ANYWHERE);
  CometModification cm(ox, 0, 3);
  TEST_REAL_SIMILAR(cm.mass, 15.994915)
  TEST_STRING_EQUAL(cm.residues, "M")
  TEST_EQUAL(cm.binary_group, 0)
  TEST_EQUAL(cm.max_mods_per_peptide, 3)
  TEST_EQUAL(cm.term_distance, -1)
  TEST_EQUAL(cm.nc_term, 0)
  TEST_EQUAL(cm.required, false)

  // Test construction from a peptide N-term modification (Acetyl)
  const ResidueModification* ac_nterm = ModificationsDB::getInstance()->getModification("Acetyl", "", ResidueModification::N_TERM);
  CometModification cm_nterm(ac_nterm, 0, 5);
  TEST_REAL_SIMILAR(cm_nterm.mass, 42.010565)
  TEST_STRING_EQUAL(cm_nterm.residues, "n")
  TEST_EQUAL(cm_nterm.term_distance, 0)
  TEST_EQUAL(cm_nterm.nc_term, 2)

  // Test construction from a peptide C-term modification
  const ResidueModification* am_cterm = ModificationsDB::getInstance()->getModification("Amidated", "", ResidueModification::C_TERM);
  CometModification cm_cterm(am_cterm, 0, 5);
  TEST_STRING_EQUAL(cm_cterm.residues, "c")
  TEST_EQUAL(cm_cterm.term_distance, 0)
  TEST_EQUAL(cm_cterm.nc_term, 3)

  // Test construction with binary group
  CometModification cm_bin(ox, 2, 5);
  TEST_EQUAL(cm_bin.binary_group, 2)
}
END_SECTION

START_SECTION(bool isMergeableWith(const CometModification& other) const)
{
  // Two amino acid mods with same mass => mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";
    a.term_distance = -1;
    a.nc_term = 0;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "R";
    b.term_distance = -1;
    b.nc_term = 0;

    TEST_EQUAL(a.isMergeableWith(b), true)
  }

  // Different mass => not mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";

    CometModification b;
    b.mass = 15.9949;
    b.residues = "M";

    TEST_EQUAL(a.isMergeableWith(b), false)
  }

  // Different binary group => not mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";
    a.binary_group = 1;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "R";
    b.binary_group = 2;

    TEST_EQUAL(a.isMergeableWith(b), false)
  }

  // Peptide N-term + amino acid => mergeable
  {
    CometModification nterm;
    nterm.mass = 43.0058;
    nterm.residues = "n";
    nterm.term_distance = 0;
    nterm.nc_term = 2;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    TEST_EQUAL(nterm.isMergeableWith(aa), true)
  }

  // Peptide C-term + amino acid => mergeable
  {
    CometModification cterm;
    cterm.mass = 43.0058;
    cterm.residues = "c";
    cterm.term_distance = 0;
    cterm.nc_term = 3;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    TEST_EQUAL(cterm.isMergeableWith(aa), true)
  }

  // Peptide N-term + peptide C-term => NOT mergeable
  {
    CometModification nterm;
    nterm.mass = 43.0058;
    nterm.residues = "n";
    nterm.term_distance = 0;
    nterm.nc_term = 2;

    CometModification cterm;
    cterm.mass = 43.0058;
    cterm.residues = "c";
    cterm.term_distance = 0;
    cterm.nc_term = 3;

    TEST_EQUAL(nterm.isMergeableWith(cterm), false)
  }

  // Protein N-term + amino acid => NOT mergeable
  {
    CometModification prot_nterm;
    prot_nterm.mass = 43.0058;
    prot_nterm.residues = "n";
    prot_nterm.term_distance = 0;
    prot_nterm.nc_term = 0;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    TEST_EQUAL(prot_nterm.isMergeableWith(aa), false)
  }

  // Protein C-term + amino acid => NOT mergeable
  {
    CometModification prot_cterm;
    prot_cterm.mass = 43.0058;
    prot_cterm.residues = "c";
    prot_cterm.term_distance = 0;
    prot_cterm.nc_term = 1;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    TEST_EQUAL(prot_cterm.isMergeableWith(aa), false)
  }

  // Protein N-term + protein N-term (same type) => mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "n";
    a.term_distance = 0;
    a.nc_term = 0;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "K";
    b.term_distance = 0;
    b.nc_term = 0;

    TEST_EQUAL(a.isMergeableWith(b), true)
  }

  // Protein N-term + protein C-term => NOT mergeable
  {
    CometModification prot_nterm;
    prot_nterm.mass = 43.0058;
    prot_nterm.residues = "n";
    prot_nterm.term_distance = 0;
    prot_nterm.nc_term = 0;

    CometModification prot_cterm;
    prot_cterm.mass = 43.0058;
    prot_cterm.residues = "c";
    prot_cterm.term_distance = 0;
    prot_cterm.nc_term = 1;

    TEST_EQUAL(prot_nterm.isMergeableWith(prot_cterm), false)
  }

  // Mass within tolerance => mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";

    CometModification b;
    b.mass = 43.0058 + 5e-7; // within 1e-6
    b.residues = "R";

    TEST_EQUAL(a.isMergeableWith(b), true)
  }

  // Mass outside tolerance => not mergeable
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";

    CometModification b;
    b.mass = 43.0058 + 2e-6; // outside 1e-6
    b.residues = "R";

    TEST_EQUAL(a.isMergeableWith(b), false)
  }
}
END_SECTION

START_SECTION(void merge(const CometModification& other))
{
  // Merge two amino acid mods
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";
    a.term_distance = -1;
    a.nc_term = 0;
    a.max_mods_per_peptide = 3;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "R";
    b.term_distance = -1;
    b.nc_term = 0;
    b.max_mods_per_peptide = 5;

    a.merge(b);
    TEST_STRING_EQUAL(a.residues, "KR")
    TEST_EQUAL(a.term_distance, -1)
    TEST_EQUAL(a.nc_term, 0)
    TEST_EQUAL(a.max_mods_per_peptide, 5) // takes max
  }

  // Merge three amino acid mods (K, R, S)
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";
    a.term_distance = -1;
    a.nc_term = 0;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "R";
    b.term_distance = -1;
    b.nc_term = 0;

    CometModification c;
    c.mass = 43.0058;
    c.residues = "S";
    c.term_distance = -1;
    c.nc_term = 0;

    a.merge(b);
    a.merge(c);
    TEST_STRING_EQUAL(a.residues, "KRS")
  }

  // Merge peptide N-term with amino acid
  // Expected: "nK" with term_distance=-1 and nc_term=0 (Comet convention)
  {
    CometModification nterm;
    nterm.mass = 43.0058;
    nterm.residues = "n";
    nterm.term_distance = 0;
    nterm.nc_term = 2;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    nterm.merge(aa);
    TEST_STRING_EQUAL(nterm.residues, "nK")
    TEST_EQUAL(nterm.term_distance, -1) // merged: no distance constraint
    TEST_EQUAL(nterm.nc_term, 0)        // nc_term=0 per Comet convention when term_distance=-1
  }

  // Merge amino acid with peptide N-term (reverse order)
  {
    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    CometModification nterm;
    nterm.mass = 43.0058;
    nterm.residues = "n";
    nterm.term_distance = 0;
    nterm.nc_term = 2;

    aa.merge(nterm);
    TEST_STRING_EQUAL(aa.residues, "Kn")
    TEST_EQUAL(aa.term_distance, -1)
    TEST_EQUAL(aa.nc_term, 0) // nc_term=0 per Comet convention
  }

  // Merge peptide N-term with multiple amino acids
  // Matches the example from the issue discussion: "nRST"
  {
    CometModification nterm;
    nterm.mass = 43.0058;
    nterm.residues = "n";
    nterm.term_distance = 0;
    nterm.nc_term = 2;

    CometModification r;
    r.mass = 43.0058;
    r.residues = "R";
    r.term_distance = -1;

    CometModification s;
    s.mass = 43.0058;
    s.residues = "S";
    s.term_distance = -1;

    CometModification t;
    t.mass = 43.0058;
    t.residues = "T";
    t.term_distance = -1;

    nterm.merge(r);
    nterm.merge(s);
    nterm.merge(t);
    TEST_STRING_EQUAL(nterm.residues, "nRST")
    TEST_EQUAL(nterm.term_distance, -1)
    TEST_EQUAL(nterm.nc_term, 0)
  }

  // Merge peptide C-term with amino acid
  {
    CometModification cterm;
    cterm.mass = 43.0058;
    cterm.residues = "c";
    cterm.term_distance = 0;
    cterm.nc_term = 3;

    CometModification aa;
    aa.mass = 43.0058;
    aa.residues = "K";
    aa.term_distance = -1;
    aa.nc_term = 0;

    cterm.merge(aa);
    TEST_STRING_EQUAL(cterm.residues, "cK")
    TEST_EQUAL(cterm.term_distance, -1)
    TEST_EQUAL(cterm.nc_term, 0)
  }

  // Duplicate residues should not be added twice
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "KR";

    CometModification b;
    b.mass = 43.0058;
    b.residues = "K";

    a.merge(b);
    TEST_STRING_EQUAL(a.residues, "KR") // K not duplicated
  }

  // Required flag: OR logic
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";
    a.required = false;

    CometModification b;
    b.mass = 43.0058;
    b.residues = "R";
    b.required = true;

    a.merge(b);
    TEST_EQUAL(a.required, true)
  }
}
END_SECTION

START_SECTION(String toCometString(Size index) const)
{
  // Basic format check
  {
    CometModification mod;
    mod.mass = 15.9949;
    mod.residues = "M";
    mod.binary_group = 0;
    mod.max_mods_per_peptide = 3;
    mod.term_distance = -1;
    mod.nc_term = 0;
    mod.required = false;

    String result = mod.toCometString(1);
    TEST_STRING_EQUAL(result, "variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0")
  }

  // Index zero-padding: single digit
  {
    CometModification mod;
    mod.mass = 79.966331;
    mod.residues = "STY";
    mod.binary_group = 0;
    mod.max_mods_per_peptide = 3;
    mod.term_distance = -1;
    mod.nc_term = 0;

    String result = mod.toCometString(9);
    TEST_EQUAL(result.hasPrefix("variable_mod09"), true)
  }

  // Index zero-padding: double digit
  {
    CometModification mod;
    mod.mass = 79.966331;
    mod.residues = "STY";
    mod.binary_group = 0;
    mod.max_mods_per_peptide = 3;
    mod.term_distance = -1;
    mod.nc_term = 0;

    String result = mod.toCometString(10);
    TEST_EQUAL(result.hasPrefix("variable_mod10"), true)
  }

  // Terminal modification with amino acids
  // Matches Comet documentation example: acetylation on K and N-terminus
  {
    CometModification mod;
    mod.mass = 42.010565;
    mod.residues = "nK";
    mod.binary_group = 0;
    mod.max_mods_per_peptide = 3;
    mod.term_distance = -1;
    mod.nc_term = 0;

    String result = mod.toCometString(1);
    TEST_STRING_EQUAL(result, "variable_mod01 = 42.0106 nK 0 3 -1 0 0 0.0")
  }

  // Required modification
  {
    CometModification mod;
    mod.mass = 15.9949;
    mod.residues = "M";
    mod.binary_group = 0;
    mod.max_mods_per_peptide = 3;
    mod.term_distance = -1;
    mod.nc_term = 0;
    mod.required = true;

    String result = mod.toCometString(1);
    TEST_STRING_EQUAL(result, "variable_mod01 = 15.9949 M 0 3 -1 0 1 0.0")
  }
}
END_SECTION

START_SECTION(static std::vector<CometModification> mergeModifications(const std::vector<CometModification>& mods))
{
  // Three mods with same mass => merge to one
  {
    CometModification k;
    k.mass = 43.0058;
    k.residues = "K";
    k.term_distance = -1;
    k.nc_term = 0;

    CometModification r;
    r.mass = 43.0058;
    r.residues = "R";
    r.term_distance = -1;
    r.nc_term = 0;

    CometModification s;
    s.mass = 43.0058;
    s.residues = "S";
    s.term_distance = -1;
    s.nc_term = 0;

    vector<CometModification> input = {k, r, s};
    vector<CometModification> result = CometModification::mergeModifications(input);

    TEST_EQUAL(result.size(), 1)
    TEST_STRING_EQUAL(result[0].residues, "KRS")
  }

  // Different masses => no merge
  {
    CometModification a;
    a.mass = 43.0058;
    a.residues = "K";

    CometModification b;
    b.mass = 15.9949;
    b.residues = "M";

    vector<CometModification> input = {a, b};
    vector<CometModification> result = CometModification::mergeModifications(input);

    TEST_EQUAL(result.size(), 2)
  }

  // N-term + amino acids merge; protein N-term stays separate
  // Input: Protein N-term, Peptide N-term, K, R, S (all same mass)
  // Expected: 2 entries (protein N-term separate, rest merged into "nKRS")
  {
    CometModification prot_nterm;
    prot_nterm.mass = 43.0058;
    prot_nterm.residues = "n";
    prot_nterm.term_distance = 0;
    prot_nterm.nc_term = 0;

    CometModification pep_nterm;
    pep_nterm.mass = 43.0058;
    pep_nterm.residues = "n";
    pep_nterm.term_distance = 0;
    pep_nterm.nc_term = 2;

    CometModification k;
    k.mass = 43.0058;
    k.residues = "K";
    k.term_distance = -1;
    k.nc_term = 0;

    CometModification r;
    r.mass = 43.0058;
    r.residues = "R";
    r.term_distance = -1;
    r.nc_term = 0;

    CometModification s;
    s.mass = 43.0058;
    s.residues = "S";
    s.term_distance = -1;
    s.nc_term = 0;

    vector<CometModification> input = {prot_nterm, pep_nterm, k, r, s};
    vector<CometModification> result = CometModification::mergeModifications(input);

    TEST_EQUAL(result.size(), 2)

    // First entry: protein N-term (stays separate)
    TEST_STRING_EQUAL(result[0].residues, "n")
    TEST_EQUAL(result[0].term_distance, 0)
    TEST_EQUAL(result[0].nc_term, 0)

    // Second entry: merged peptide N-term + amino acids
    TEST_EQUAL(result[1].residues.hasSubstring("n"), true)
    TEST_EQUAL(result[1].residues.hasSubstring("K"), true)
    TEST_EQUAL(result[1].residues.hasSubstring("R"), true)
    TEST_EQUAL(result[1].residues.hasSubstring("S"), true)
    TEST_EQUAL(result[1].term_distance, -1)
    TEST_EQUAL(result[1].nc_term, 0)
  }

  // Empty input => empty output
  {
    vector<CometModification> input;
    vector<CometModification> result = CometModification::mergeModifications(input);
    TEST_EQUAL(result.size(), 0)
  }

  // Single input => same output
  {
    CometModification a;
    a.mass = 15.9949;
    a.residues = "M";

    vector<CometModification> input = {a};
    vector<CometModification> result = CometModification::mergeModifications(input);

    TEST_EQUAL(result.size(), 1)
    TEST_STRING_EQUAL(result[0].residues, "M")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
