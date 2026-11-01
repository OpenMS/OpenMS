// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLSiteLocalization.h"
///////////////////////////

#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  // per-residue intensities of shifted / unshifted prefix (b) and suffix (y) ions
  struct Evidence
  {
    explicit Evidence(size_t length) : n_shifts(length, 0.0), n_noshifts(length, 0.0), c_shifts(length, 0.0), c_noshifts(length, 0.0) {}
    vector<double> n_shifts, n_noshifts, c_shifts, c_noshifts;
  };

  // b ions place the cross-link at residues first..last: b(first) unshifted, b(last + 1) shifted
  Evidence prefixEvidence(size_t length, int first, int last)
  {
    Evidence e(length);
    if (first > 0) { e.n_noshifts[first - 1] = 1.0; }
    e.n_shifts[last] = 1.0;
    return e;
  }

  NuXLSiteLocalization::Site localize(const string& sequence, bool reactive_protein_n_term, const Evidence& e)
  {
    return NuXLSiteLocalization::localizeNHS(sequence, reactive_protein_n_term, e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
  }
}

START_TEST(NuXLSiteLocalization, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::pair<int, int> supportedRange(const std::vector<double>& n_shifts, const std::vector<double>& n_noshifts, const std::vector<double>& c_shifts, const std::vector<double>& c_noshifts)))
{
  // no shifted ion: no evidence for any site
  {
    Evidence e(8);
    e.n_noshifts[2] = 1.0;
    e.c_noshifts[5] = 1.0;
    TEST_EQUAL(NuXLSiteLocalization::supportedRange(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == make_pair(-1, -1), true)
  }

  // b3 unshifted, b5 shifted: residues 3 and 4 (0-based 2-4 without b4)
  {
    const Evidence e = prefixEvidence(8, 2, 4);
    TEST_EQUAL(NuXLSiteLocalization::supportedRange(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == make_pair(2, 4), true)
  }

  // only b2 shifted: without b1 and y7, the first two residues can't be told apart (KgAIVTGK in MSFragger)
  {
    const Evidence e = prefixEvidence(8, 0, 1);
    TEST_EQUAL(NuXLSiteLocalization::supportedRange(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == make_pair(0, 1), true)
  }

  // y7 shifted, y6 unshifted: the second residue, not the first
  {
    Evidence e(8);
    e.c_shifts[1] = 1.0;
    e.c_noshifts[2] = 1.0;
    TEST_EQUAL(NuXLSiteLocalization::supportedRange(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == make_pair(1, 1), true)
  }

  // contradicting ions: the stronger ones decide
  {
    Evidence e = prefixEvidence(8, 2, 4);
    e.c_shifts[6] = 0.2; // a weak y2 that would put the site at the C-terminus
    TEST_EQUAL(NuXLSiteLocalization::supportedRange(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == make_pair(2, 4), true)
  }
}
END_SECTION

START_SECTION((static Site localizeNHS(const std::string& sequence, bool reactive_protein_n_term, const std::vector<double>& n_shifts, const std::vector<double>& n_noshifts, const std::vector<double>& c_shifts, const std::vector<double>& c_noshifts)))
{
  // lysine at the peptide N-terminus (KgAIVTGK)
  {
    const NuXLSiteLocalization::Site s = localize("KGAIVTGK", false, prefixEvidence(8, 0, 1));
    TEST_EQUAL(s.position, 0)
    TEST_STRING_EQUAL(s.method, "K")
    TEST_EQUAL(s.range_begin, 0)
    TEST_EQUAL(s.range_end, 1)
  }

  // the free protein N-terminus comes first
  {
    const NuXLSiteLocalization::Site s = localize("SAKAAR", true, prefixEvidence(6, 0, 2));
    TEST_EQUAL(s.position, 0)
    TEST_STRING_EQUAL(s.method, "N-term")
  }

  // ... but only if the fragment ions allow the first residue
  {
    const NuXLSiteLocalization::Site s = localize("SAKAAR", true, prefixEvidence(6, 2, 3));
    TEST_EQUAL(s.position, 2)
    TEST_STRING_EQUAL(s.method, "K")
  }

  // of two lysines, the one closer to the middle of the range
  {
    const NuXLSiteLocalization::Site s = localize("AAAKAAKAAR", false, prefixEvidence(10, 2, 7));
    TEST_EQUAL(s.position, 3)
    TEST_STRING_EQUAL(s.method, "K")
  }

  // equally close: the one closer to the N-terminus
  {
    const NuXLSiteLocalization::Site s = localize("AAKAKAAR", false, prefixEvidence(8, 2, 4));
    TEST_EQUAL(s.position, 2)
  }

  // no lysine: serine next to a histidine
  {
    const NuXLSiteLocalization::Site s = localize("LLSHAAGR", false, prefixEvidence(8, 1, 2));
    TEST_EQUAL(s.position, 2)
    TEST_STRING_EQUAL(s.method, "STY+H")
  }

  // serine without a histidine within four residues: unlocalized
  {
    const NuXLSiteLocalization::Site s = localize("LLSAAAAHR", false, prefixEvidence(9, 1, 2));
    TEST_EQUAL(s.position, -1)
    TEST_STRING_EQUAL(s.method, "unlocalized")
    TEST_EQUAL(s.range_begin, 1)
    TEST_EQUAL(s.range_end, 2)
  }

  // the fragment ions point to a residue NHS esters don't react with: unlocalized (TKphVNVGTIGHVDHGK)
  {
    const NuXLSiteLocalization::Site s = localize("TKPHVNVGTIGHVDHGK", false, prefixEvidence(17, 2, 3));
    TEST_EQUAL(s.position, -1)
    TEST_STRING_EQUAL(s.method, "unlocalized")
  }

  // no fragment evidence: unlocalized without a range
  {
    const NuXLSiteLocalization::Site s = localize("KGAIVTGK", true, Evidence(8));
    TEST_EQUAL(s.position, -1)
    TEST_EQUAL(s.range_begin, -1)
    TEST_EQUAL(s.range_end, -1)
  }
}
END_SECTION

START_SECTION((static std::vector<std::string> proteinNTermini(const std::vector<std::string_view>& protein_sequences)))
{
  const vector<string> n_termini = NuXLSiteLocalization::proteinNTermini({"PEPTIDEK", "MSAEKLLR"});
  TEST_EQUAL(n_termini.size(), 3)
  ABORT_IF(n_termini.size() != 3)
  TEST_STRING_EQUAL(n_termini[0], "MSAEKLLR")
  TEST_STRING_EQUAL(n_termini[1], "PEPTIDEK")
  TEST_STRING_EQUAL(n_termini[2], "SAEKLLR") // initiator methionine removed
}
END_SECTION

START_SECTION((static bool isProteinNTerminal(const std::vector<std::string>& protein_n_termini, std::string_view peptide)))
{
  const vector<string> n_termini = NuXLSiteLocalization::proteinNTermini({"PEPTIDEK", "MSAEKLLR"});
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "MSAEK"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "SAEK"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "PEPTIDEK"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "AEKLLR"), false)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "EPTIDEK"), false)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(n_termini, "PEPTIDEKR"), false) // longer than the protein
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal({}, "PEPTIDEK"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
