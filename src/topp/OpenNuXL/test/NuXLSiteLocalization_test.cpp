// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

#include <set>
#include <string>
#include <string_view>
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

  // residues first..last; empty for first < 0
  vector<int> residues(int first, int last)
  {
    vector<int> r;
    for (int i = first; first >= 0 && i <= last; ++i) { r.push_back(i); }
    return r;
  }

  NuXLSiteLocalization::Site localize(const string& sequence, bool reactive_protein_n_term, const Evidence& e)
  {
    return NuXLSiteLocalization::localizeNHS(sequence, reactive_protein_n_term, e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
  }
}

START_TEST(NuXLSiteLocalization, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::vector<int> supportedResidues(const std::vector<double>& n_shifts, const std::vector<double>& n_noshifts, const std::vector<double>& c_shifts, const std::vector<double>& c_noshifts, double tolerance = DEFAULT_TOLERANCE)))
{
  // no fragment ion: no evidence for any site
  {
    Evidence e(8);
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(-1, -1), true)
  }

  // unshifted ions alone: b3 and y3 unshifted, the cross-link is on residues 4-5 (0-based 3-4)
  {
    Evidence e(8);
    e.n_noshifts[2] = 1.0;
    e.c_noshifts[5] = 1.0;
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(3, 4), true)
  }

  // an unshifted b1 is ignored (b1 of lysine coincides with common lysine ions): kLLDEGR with b2 shifted
  {
    Evidence e = prefixEvidence(7, 0, 1);
    e.n_noshifts[0] = 5.0;
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(0, 1), true)
  }

  // b3 unshifted, b5 shifted: residues 3 and 4 (0-based 2-4 without b4)
  {
    const Evidence e = prefixEvidence(8, 2, 4);
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(2, 4), true)
  }

  // only b2 shifted: without y7, the first two residues can't be told apart (KgAIVTGK in MSFragger)
  {
    const Evidence e = prefixEvidence(8, 0, 1);
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(0, 1), true)
  }

  // y7 shifted, y6 unshifted: the second residue, not the first
  {
    Evidence e(8);
    e.c_shifts[1] = 1.0;
    e.c_noshifts[2] = 1.0;
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(1, 1), true)
  }

  // a weak contradicting ion (below the tolerance) doesn't decide: NDKNGIIHTTIGK with a weak a2+SDA
  {
    Evidence e(13);
    e.n_shifts[1] = 0.02; // a2+SDA: site at residue 1-2
    e.n_shifts[2] = 1.0;  // b3+SDA: site at residue 1-3
    e.c_noshifts[3] = 1.0; // y10 unshifted: site at residue 1-3
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(0, 2), true)
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts, 0.0) == residues(0, 1), true)
    const NuXLSiteLocalization::Site s = localize("NDKNGIIHTTIGK", false, e);
    TEST_EQUAL(s.position, 2)
    TEST_STRING_EQUAL(s.method, "K")
  }

  // contradicting ions of equal intensity support both ends of the peptide, but not the residues between them
  {
    Evidence e(8);
    e.n_shifts[1] = 1.0; // b2 shifted: residue 1-2
    e.c_shifts[6] = 1.0; // y2 shifted: residue 7-8
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == vector<int>({0, 1, 6, 7}), true)
    const NuXLSiteLocalization::Site s = localize("AAAKAAAR", false, e); // the lysine is not supported
    TEST_EQUAL(s.position, -1)
    TEST_EQUAL(s.range_begin, 0)
    TEST_EQUAL(s.range_end, 7)
  }

  // contradicting ions: the stronger ones decide
  {
    Evidence e = prefixEvidence(8, 2, 4);
    e.c_shifts[6] = 0.2; // a weak y2 that would put the site at the C-terminus
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(2, 4), true)
  }
}
END_SECTION

START_SECTION((static Site localizeNHS(const std::string& sequence, bool reactive_protein_n_term, const std::vector<double>& n_shifts, const std::vector<double>& n_noshifts, const std::vector<double>& c_shifts, const std::vector<double>& c_noshifts, double tolerance = DEFAULT_TOLERANCE, const std::vector<bool>& blocked = {})))
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

  // a lysine whose amine is blocked (e.g. acetylated) is skipped
  {
    const Evidence e = prefixEvidence(10, 2, 7);
    vector<bool> blocked(10, false);
    blocked[3] = true;
    const NuXLSiteLocalization::Site s = NuXLSiteLocalization::localizeNHS("AAAKAAKAAR", false, e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts,
      NuXLSiteLocalization::DEFAULT_TOLERANCE, blocked);
    TEST_EQUAL(s.position, 6)
    TEST_STRING_EQUAL(s.method, "K")
  }

  // a blocked C-terminal lysine (not cleaved by trypsin if modified): the threonine next to the histidine
  {
    const Evidence e = prefixEvidence(7, 3, 6);
    vector<bool> blocked(7, false);
    blocked[6] = true;
    const NuXLSiteLocalization::Site s = NuXLSiteLocalization::localizeNHS("HFAATGK", false, e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts,
      NuXLSiteLocalization::DEFAULT_TOLERANCE, blocked);
    TEST_EQUAL(s.position, 4)
    TEST_STRING_EQUAL(s.method, "STY+H")
    TEST_EQUAL(s.range_begin, 3)
    TEST_EQUAL(s.range_end, 6)
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
    TEST_EQUAL(s.range_begin, 0) // the unshifted b1 doesn't count
    TEST_EQUAL(s.range_end, 2)
  }

  // the fragment ions point to a residue NHS esters don't react with: unlocalized (TKphVNVGTIGHVDHGK)
  {
    const NuXLSiteLocalization::Site s = localize("TKPHVNVGTIGHVDHGK", false, prefixEvidence(17, 2, 3));
    TEST_EQUAL(s.position, -1)
    TEST_STRING_EQUAL(s.method, "unlocalized")
  }

  // lysine found by the unshifted y ions alone (vtkIFVDEGPSMK: y10 unshifted)
  {
    Evidence e(13);
    e.c_noshifts[3] = 1.0;
    const NuXLSiteLocalization::Site s = localize("VTKIFVDEGPSMK", false, e);
    TEST_EQUAL(s.position, 2)
    TEST_STRING_EQUAL(s.method, "K")
    TEST_EQUAL(s.range_begin, 0)
    TEST_EQUAL(s.range_end, 2)
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

START_SECTION((static std::vector<double> consistencyScores(const std::vector<double>& n_shifts, const std::vector<double>& n_noshifts, const std::vector<double>& c_shifts, const std::vector<double>& c_noshifts, double unshifted_weight = DEFAULT_UNSHIFTED_WEIGHT)))
{
  // no fragment ion: no evidence for any site
  {
    const Evidence e(8);
    TEST_EQUAL(NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts).empty(), true)
  }

  // b2 unshifted, b5 shifted: b5 counts for residues 0-4 and against 5-7, b2 (weight 0.5) against 0-1
  {
    const Evidence e = prefixEvidence(8, 2, 4);
    const vector<double> s = NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
    TEST_EQUAL(s.size(), 8)
    TEST_REAL_SIMILAR(s[0], 0.5 / 1.5)
    TEST_REAL_SIMILAR(s[1], 0.5 / 1.5)
    TEST_REAL_SIMILAR(s[2], 1.0 / 1.5)
    TEST_REAL_SIMILAR(s[4], 1.0 / 1.5)
    TEST_REAL_SIMILAR(s[5], -1.0 / 1.5)
    TEST_REAL_SIMILAR(s[7], -1.0 / 1.5)
  }

  // an unshifted ion of a rung that also has a shifted one lost the nucleic acid: it doesn't count
  {
    Evidence e(8);
    e.c_shifts[3] = 1.0;   // y5 shifted: residues 3-7
    e.c_noshifts[4] = 1.0; // y4 unshifted (no shifted y4): against residues 4-7
    const vector<double> without_loss = NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
    e.c_noshifts[3] = 2.0; // y5 also unshifted, even more intense
    const vector<double> with_loss = NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
    TEST_EQUAL(with_loss == without_loss, true)
    TEST_REAL_SIMILAR(with_loss[2], -1.0 / 1.5)
    TEST_REAL_SIMILAR(with_loss[3], 1.0 / 1.5)
    TEST_REAL_SIMILAR(with_loss[4], 0.5 / 1.5)
    // supportedResidues takes y5 for evidence against residue 3
    TEST_EQUAL(NuXLSiteLocalization::supportedResidues(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts) == residues(0, 2), true)
  }

  // unshifted ions alone: b3 and y3 unshifted, residues 3-4 contradicted by none
  {
    Evidence e(8);
    e.n_noshifts[2] = 1.0;
    e.c_noshifts[5] = 1.0;
    const vector<double> s = NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts);
    TEST_EQUAL(s == vector<double>({-0.5, -0.5, -0.5, 0.0, 0.0, -0.5, -0.5, -0.5}), true)
    // ... count only with a weight above 0
    TEST_EQUAL(NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts, 0.0).empty(), true)
  }

  // an unshifted b1 is ignored (b1 of lysine coincides with common lysine ions)
  {
    Evidence e(7);
    e.n_noshifts[0] = 5.0;
    TEST_EQUAL(NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts).empty(), true)
  }
}
END_SECTION

START_SECTION((static Site localizeByConsistency(const std::vector<double>& scores, double tolerance = DEFAULT_TOLERANCE, const std::vector<bool>& blocked = {})))
{
  const auto locate = [](const Evidence& e, const vector<bool>& blocked = {}, double unshifted_weight = NuXLSiteLocalization::DEFAULT_UNSHIFTED_WEIGHT)
  {
    return NuXLSiteLocalization::localizeByConsistency(
      NuXLSiteLocalization::consistencyScores(e.n_shifts, e.n_noshifts, e.c_shifts, e.c_noshifts, unshifted_weight),
      NuXLSiteLocalization::DEFAULT_TOLERANCE, blocked);
  };

  // no fragment evidence: unlocalized without a range
  {
    const NuXLSiteLocalization::Site s = NuXLSiteLocalization::localizeByConsistency({});
    TEST_EQUAL(s.position, -1)
    TEST_EQUAL(s.range_begin, -1)
    TEST_EQUAL(s.range_end, -1)
    TEST_STRING_EQUAL(s.method, "unlocalized")
    TEST_REAL_SIMILAR(s.delta, 0.0)
  }

  // residues no ion separates tie: the one closer to the N-terminus, both in the range, no delta
  {
    const NuXLSiteLocalization::Site s = locate(prefixEvidence(8, 2, 4));
    TEST_EQUAL(s.position, 2)
    TEST_EQUAL(s.range_begin, 2)
    TEST_EQUAL(s.range_end, 4)
    TEST_STRING_EQUAL(s.method, "consistency")
    TEST_REAL_SIMILAR(s.delta, 0.0)
  }

  // a weak misassigned shifted ion doesn't pull the site towards the terminus (as when summing all longer shifted ions)
  {
    Evidence e(10);
    e.n_noshifts[2] = 1.0; // b3 unshifted: against residues 0-2 (weight 0.5)
    e.n_shifts[3] = 1.0;   // b4 shifted: residues 0-3
    e.c_shifts[3] = 1.0;   // y7 shifted: residues 3-9
    e.c_noshifts[4] = 1.0; // y6 unshifted: against residues 4-9 (weight 0.5)
    e.c_shifts[8] = 0.05;  // weak y2 shifted: residues 8-9
    const NuXLSiteLocalization::Site s = locate(e);
    TEST_EQUAL(s.position, 3)
    TEST_EQUAL(s.range_begin, 3)
    TEST_EQUAL(s.range_end, 3)
    TEST_REAL_SIMILAR(s.delta, (1.95 - -0.45) / 3.05) // residue 3: 1 + 1 - 0.05, residues 8 and 9: -1 + 1 - 0.5 + 0.05
  }

  // a blocked residue (e.g. a C-terminal lysine the enzyme would have cleaved after) is no candidate: only y2 shifted
  {
    Evidence e(6);
    e.c_shifts[4] = 1.0;
    NuXLSiteLocalization::Site s = locate(e);
    TEST_EQUAL(s.position, 4)
    TEST_EQUAL(s.range_end, 5)
    TEST_REAL_SIMILAR(s.delta, 0.0)
    vector<bool> blocked(6, false);
    blocked[5] = true;
    s = locate(e, blocked);
    TEST_EQUAL(s.position, 4)
    TEST_EQUAL(s.range_begin, 4)
    TEST_EQUAL(s.range_end, 4)
    TEST_REAL_SIMILAR(s.delta, 2.0)
    s = locate(e, vector<bool>(6, true));
    TEST_EQUAL(s.position, -1)
    TEST_STRING_EQUAL(s.method, "unlocalized")
  }

  // FDQKIIDK+U-H2O (TOPP test OpenNuXL_1): shifted b ions start at K3 and shifted y ions end there. The neighbour rules
  // of NuXL:localization 'fragments' reject every ladder site there, so that the shifted lysine immonium ions alone
  // decide and their tie goes to the C-terminal K7.
  {
    Evidence e(8);
    e.n_shifts = {0.0, 0.0, 0.0, 0.146, 0.241, 0.081, 0.087, 0.0};
    e.n_noshifts = {0.271, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    e.c_shifts = {0.0, 0.213, 0.056, 0.064, 0.0, 0.0, 0.0, 0.0};
    e.c_noshifts = {0.0, 0.0, 0.0, 0.0, 0.0, 0.124, 0.279, 0.256};
    const NuXLSiteLocalization::Site s = locate(e);
    TEST_EQUAL(s.position, 3)
    // only the shifted y5 separates K3 from Q2; all ions: 0.888 shifted, 0.125 + 0.659 unshifted (weight 0.5)
    TEST_REAL_SIMILAR(s.delta, 2 * 0.064 / 1.28)
  }

  // DEPMITLPR+U (TOPP test OpenNuXL_1): the shifted b3 and b4 (0.129) put the cross-link before I4, the about as
  // intense shifted y6 and y5 (0.128) after P2: P2-I4 can't be told apart
  {
    Evidence e(9);
    e.n_shifts = {0.0, 0.0, 0.100, 0.029, 0.021, 0.0, 0.0, 0.022, 0.0};
    e.n_noshifts = {0.0, 0.496, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    e.c_shifts = {0.0, 0.0, 0.534, 0.090, 0.038, 0.0, 0.0, 0.0, 0.0};
    e.c_noshifts = {0.0, 0.0, 0.038, 0.269, 0.483, 0.700, 0.262, 0.503, 0.147};
    const NuXLSiteLocalization::Site s = locate(e);
    TEST_EQUAL(s.position, 2)
    TEST_EQUAL(s.range_begin, 2)
    TEST_EQUAL(s.range_end, 4)
    TEST_EQUAL(s.delta < 0.01, true)
  }

  // MVFCLSAGK+U-H2O (TOPP data OpenNuXL_1_input): shifted y7 and y6, unshifted y5-y1. Without counting unshifted ions,
  // C3-G7 can't be told apart (the C-terminal lysine blocked); with them, the unshifted y5 and y4 contradict L4 and S5.
  {
    Evidence e(9);
    e.n_noshifts = {0.0, 0.847, 0.059, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    e.c_shifts = {0.0, 0.0, 0.341, 0.195, 0.0, 0.0, 0.0, 0.0, 0.0};
    e.c_noshifts = {0.0, 0.0, 0.0, 0.0, 0.123, 0.197, 0.130, 0.293, 0.216};
    vector<bool> blocked(9, false);
    blocked[8] = true;
    NuXLSiteLocalization::Site s = locate(e, blocked, 0.0);
    TEST_EQUAL(s.position, 3)
    TEST_EQUAL(s.range_begin, 3)
    TEST_EQUAL(s.range_end, 7)
    TEST_REAL_SIMILAR(s.delta, 0.0)
    s = locate(e, blocked);
    TEST_EQUAL(s.position, 3)
    TEST_EQUAL(s.range_begin, 3)
    TEST_EQUAL(s.range_end, 4)
    TEST_EQUAL(s.delta > 0.0, true)
  }
}
END_SECTION

START_SECTION((static void digestProtein(const EnzymaticDigestion& digestion, std::string_view protein, size_t min_length, size_t max_length, bool clip_nterm_methionine, std::vector<std::string_view>& peptides)))
{
  ProteaseDigestion trypsin;
  trypsin.setEnzyme("Trypsin");
  trypsin.setMissedCleavages(1);
  const string protein = "MSAEKLLRPEK";
  vector<string_view> peptides;

  NuXLSiteLocalization::digestProtein(trypsin, protein, 3, 0, false, peptides);
  set<string> digest(peptides.begin(), peptides.end());
  TEST_EQUAL(digest == set<string>({"MSAEK", "LLRPEK", "MSAEKLLRPEK"}), true)

  // initial methionine removed: also the N-terminal peptides of SAEKLLRPEK, as views into the protein
  NuXLSiteLocalization::digestProtein(trypsin, protein, 3, 0, true, peptides);
  digest = set<string>(peptides.begin(), peptides.end());
  TEST_EQUAL(digest == set<string>({"MSAEK", "LLRPEK", "MSAEKLLRPEK", "SAEK", "SAEKLLRPEK"}), true)
  for (const string_view& p : peptides) { TEST_EQUAL(p.data() >= protein.data() && p.data() + p.size() <= protein.data() + protein.size(), true) }

  // no initial methionine: nothing to clip
  NuXLSiteLocalization::digestProtein(trypsin, "PEPTIDEKAAR", 3, 0, true, peptides);
  TEST_EQUAL(peptides.size(), 3)
}
END_SECTION

START_SECTION((static ProteinTermini proteinTermini(const std::vector<std::string_view>& proteins, const EnzymaticDigestion& digestion, size_t min_length, size_t max_length, bool clip_nterm_methionine)))
{
  ProteaseDigestion trypsin;
  trypsin.setEnzyme("Trypsin");
  trypsin.setMissedCleavages(0);
  const NuXLSiteLocalization::ProteinTermini termini = NuXLSiteLocalization::proteinTermini({"MACDEKPGGR", "RACDEKAAAAK"}, trypsin, 3, 0, true);
  // MACDEKPGGR: no cleavage before P, so its N-terminal peptides are MACDEKPGGR and, without the methionine, ACDEKPGGR
  TEST_EQUAL(termini.n_terminal == vector<string>({"ACDEKPGGR", "MACDEKPGGR"}), true)
  TEST_EQUAL(termini.c_terminal == vector<string>({"AAAAK", "ACDEKPGGR", "MACDEKPGGR"}), true) // the clipped protein ends it, too
}
END_SECTION

START_SECTION((static bool isProteinNTerminal(const ProteinTermini& termini, std::string_view peptide)))
{
  ProteaseDigestion trypsin;
  trypsin.setEnzyme("Trypsin");
  trypsin.setMissedCleavages(0);
  const NuXLSiteLocalization::ProteinTermini termini = NuXLSiteLocalization::proteinTermini({"MACDEKPGGR", "RACDEKAAAAK"}, trypsin, 3, 0, true);
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(termini, "MACDEKPGGR"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(termini, "ACDEKPGGR"), true)
  // internal peptide of RACDEKAAAAK: a prefix of the clipped MACDEKPGGR, but not a product of its N-terminus
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(termini, "ACDEK"), false)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(termini, "AAAAK"), false)
  TEST_EQUAL(NuXLSiteLocalization::isProteinNTerminal(NuXLSiteLocalization::ProteinTermini(), "MACDEKPGGR"), false)
}
END_SECTION

START_SECTION((static bool isProteinCTerminal(const ProteinTermini& termini, std::string_view peptide)))
{
  ProteaseDigestion trypsin;
  trypsin.setEnzyme("Trypsin");
  trypsin.setMissedCleavages(0);
  const NuXLSiteLocalization::ProteinTermini termini = NuXLSiteLocalization::proteinTermini({"AKRGKSVEEILGK", "SVEEILGKAAR"}, trypsin, 3, 0, true);
  TEST_EQUAL(NuXLSiteLocalization::isProteinCTerminal(termini, "SVEEILGK"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinCTerminal(termini, "AAR"), true)
  TEST_EQUAL(NuXLSiteLocalization::isProteinCTerminal(termini, "GKSVEEILGK"), false) // a suffix, but not a product
  TEST_EQUAL(NuXLSiteLocalization::isProteinCTerminal(NuXLSiteLocalization::ProteinTermini(), "SVEEILGK"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
