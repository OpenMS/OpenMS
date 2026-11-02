// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLFDR.h"
///////////////////////////

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <map>
#include <random>
#include <string>

using namespace OpenMS;
using namespace std;

namespace
{
  const string LONG_XL = "PEPTIDEPEPTIDEK"; // 15 residues
  const string SHORT_XL = "PEPTIDK";        // 7 residues

  // identification of one spectrum with a single hit
  PeptideIdentification identification(const string& ref, const string& sequence, double score, bool decoy, bool xl, int na_length)
  {
    PeptideIdentification pi;
    pi.setSpectrumReference(ref);
    pi.setHigherScoreBetter(true);
    pi.setScoreType("NuXLScore");
    PeptideHit ph;
    ph.setSequence(AASequence::fromString(sequence));
    ph.setScore(score);
    ph.setMetaValue("target_decoy", decoy ? "decoy" : "target");
    ph.setMetaValue("NuXL:isXL", xl ? 1 : 0);
    ph.setMetaValue("NuXL:NA_length", na_length);
    pi.setHits({ph});
    return pi;
  }

  /**
    Cross-links: 300 long targets with the best scores, 100 short targets below them with @p short_decoys short decoys
    among them (one after every 8 short targets; all pass 10% FDR over all cross-links), 5 long decoys at the end.
    Peptides: 50 targets and 2 decoys. Shuffled, so that the merge has to restore the input order.
  */
  PeptideIdentificationList example(int short_decoys)
  {
    PeptideIdentificationList ids;
    int n = 0;
    auto add = [&](const string& seq, double score, bool decoy, bool xl, int na_length)
    { ids.push_back(identification("scan=" + to_string(n++), seq, score, decoy, xl, na_length)); };
    for (int i = 0; i < 300; ++i) { add(LONG_XL, 1000.0 - i, false, true, 2); }
    for (int i = 0; i < 100; ++i) { add(SHORT_XL, 700.0 - 3.0 * i, false, true, 2); }
    for (int k = 0; k < short_decoys; ++k) { add(SHORT_XL, 700.0 - 3.0 * (8 * k) - 1.0, true, true, 2); }
    for (int k = 0; k < 5; ++k) { add(LONG_XL, 100.0 - k, true, true, 2); }
    for (int i = 0; i < 50; ++i) { add("SAMPLER", 500.0 - i, false, false, 0); }
    for (int k = 0; k < 2; ++k) { add("SAMPLER", 460.5 - 10.0 * k, true, false, 0); }
    std::shuffle(ids.begin(), ids.end(), std::mt19937(7));
    return ids;
  }

  // q-values of FalseDiscoveryRate (with NuXLFDR's parameters) applied to @p ids, by spectrum reference
  map<string, double> qValues(PeptideIdentificationList ids)
  {
    FalseDiscoveryRate fdr;
    Param p = fdr.getParameters();
    p.setValue("add_decoy_proteins", "true");
    p.setValue("add_decoy_peptides", "true");
    fdr.setParameters(p);
    fdr.apply(ids, true);
    map<string, double> q;
    for (const auto& pi : ids) { q[pi.getSpectrumReference()] = pi.getHits()[0].getScore(); }
    return q;
  }

  PeptideIdentificationList select(const PeptideIdentificationList& ids, bool xl, int short_xl /* -1: any */)
  {
    PeptideIdentificationList out;
    for (const auto& pi : ids)
    {
      const PeptideHit& h = pi.getHits()[0];
      if ((static_cast<int>(h.getMetaValue("NuXL:isXL")) != 0) != xl) continue;
      if (short_xl != -1 && NuXLFDR::isShortXL(h) != (short_xl == 1)) continue;
      out.push_back(pi);
    }
    return out;
  }
}

START_TEST(NuXLFDR, "$Id$")

START_SECTION((static XLGroups xlGroupsFromString(const std::string& s)))
{
  TEST_EQUAL(NuXLFDR::xlGroupsFromString("auto") == NuXLFDR::XLGroups::AUTO, true)
  TEST_EQUAL(NuXLFDR::xlGroupsFromString("on") == NuXLFDR::XLGroups::ON, true)
  TEST_EQUAL(NuXLFDR::xlGroupsFromString("off") == NuXLFDR::XLGroups::OFF, true)
  TEST_EXCEPTION(Exception::InvalidValue, NuXLFDR::xlGroupsFromString("yes"))
}
END_SECTION

START_SECTION((static bool isShortXL(const PeptideHit& hit)))
{
  auto hit = [](const string& seq, bool xl, int na_length) { return identification("", seq, 1.0, false, xl, na_length).getHits()[0]; };
  TEST_EQUAL(NuXLFDR::isShortXL(hit("PEPTIDK", true, 2)), true)
  TEST_EQUAL(NuXLFDR::isShortXL(hit("PEPTIDEKR", true, 3)), true) // 9 residues: still short
  TEST_EQUAL(NuXLFDR::isShortXL(hit("PEPTIDEKRR", true, 2)), false) // 10 residues
  TEST_EQUAL(NuXLFDR::isShortXL(hit("PEPTIDK", true, 1)), false) // mononucleotide
  TEST_EQUAL(NuXLFDR::isShortXL(hit("PEPTIDK", false, 0)), false) // not a cross-link
  PeptideHit modified = hit("PEPM(Oxidation)IDK", true, 2); // modifications don't count as residues
  TEST_EQUAL(NuXLFDR::isShortXL(modified), true)
}
END_SECTION

START_SECTION((bool useShortXLGroup(const PeptideIdentificationList& xl_pi, Size& short_decoys) const))
{
  NuXLFDR off(1, NuXLFDR::XLGroups::OFF), on(1, NuXLFDR::XLGroups::ON), autom(1, NuXLFDR::XLGroups::AUTO);
  Size short_decoys = 99;
  const PeptideIdentificationList xl12 = select(example(12), true, -1), xl9 = select(example(9), true, -1);
  TEST_EQUAL(off.useShortXLGroup(xl12, short_decoys), false)
  TEST_EQUAL(short_decoys, 0)
  TEST_EQUAL(on.useShortXLGroup(xl9, short_decoys), true)
  TEST_EQUAL(on.useShortXLGroup(select(example(9), true, 0), short_decoys), false) // no short cross-links
  // auto: 12 short decoys pass 10% FDR over all cross-links -> own group; with 9 not
  TEST_EQUAL(autom.useShortXLGroup(xl12, short_decoys), true)
  TEST_EQUAL(short_decoys, 12)
  TEST_EQUAL(autom.useShortXLGroup(xl9, short_decoys), false)
  TEST_EQUAL(short_decoys, 9)
  TEST_EQUAL(autom.useShortXLGroup(PeptideIdentificationList(), short_decoys), false)
}
END_SECTION

START_SECTION((void calculatePeptideAndXLQValueAtPSMLevel(const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep_pi, PeptideIdentificationList& xl_pi) const))
{
  const PeptideIdentificationList ids = example(12);
  vector<string> xl_order;
  for (const auto& pi : ids)
  {
    if (static_cast<int>(pi.getHits()[0].getMetaValue("NuXL:isXL")) != 0) { xl_order.push_back(pi.getSpectrumReference()); }
  }
  const map<string, double> q_peptides = qValues(select(ids, false, -1));
  const map<string, double> q_xl = qValues(select(ids, true, -1));
  map<string, double> q_groups = qValues(select(ids, true, 1));
  for (const auto& [ref, q] : qValues(select(ids, true, 0))) { q_groups[ref] = q; }

  for (auto mode : {NuXLFDR::XLGroups::OFF, NuXLFDR::XLGroups::AUTO})
  {
    PeptideIdentificationList pep_pi, xl_pi;
    NuXLFDR(1, mode).calculatePeptideAndXLQValueAtPSMLevel(ids, pep_pi, xl_pi);
    const bool groups = mode == NuXLFDR::XLGroups::AUTO;
    const map<string, double>& expected = groups ? q_groups : q_xl;
    TEST_EQUAL(pep_pi.size(), 52)
    TEST_EQUAL(xl_pi.size(), xl_order.size())
    Size mismatches = 0, order_mismatches = 0, annotation_mismatches = 0, differ_from_class_wide = 0;
    for (Size i = 0; i < xl_pi.size(); ++i)
    {
      const string ref = xl_pi[i].getSpectrumReference();
      const PeptideHit& h = xl_pi[i].getHits()[0];
      if (ref != xl_order[i]) ++order_mismatches;
      if (std::fabs(h.getScore() - expected.at(ref)) > 1e-12) ++mismatches;
      if (std::fabs(h.getScore() - q_xl.at(ref)) > 1e-12) ++differ_from_class_wide;
      const std::string group = h.getMetaValue("NuXL:fdr_group", "").toString();
      if (group != (groups ? (NuXLFDR::isShortXL(h) ? "short_XL" : "XL") : "")) ++annotation_mismatches;
    }
    TEST_EQUAL(mismatches, 0)
    TEST_EQUAL(order_mismatches, 0) // input order kept
    TEST_EQUAL(annotation_mismatches, 0)
    TEST_EQUAL(differ_from_class_wide > 0, groups) // the groups change the q-values
    // peptides are not affected
    Size peptide_mismatches = 0;
    for (const auto& pi : pep_pi)
    {
      if (std::fabs(pi.getHits()[0].getScore() - q_peptides.at(pi.getSpectrumReference())) > 1e-12) ++peptide_mismatches;
    }
    TEST_EQUAL(peptide_mismatches, 0)
  }
  // the short cross-links are worse than the others: in their own group, their q-values are higher (none lower)
  {
    PeptideIdentificationList pep_pi, xl_pi;
    NuXLFDR(1, NuXLFDR::XLGroups::ON).calculatePeptideAndXLQValueAtPSMLevel(ids, pep_pi, xl_pi);
    Size higher = 0, lower = 0, short_targets = 0;
    for (const auto& pi : xl_pi)
    {
      const PeptideHit& h = pi.getHits()[0];
      if (!NuXLFDR::isShortXL(h) || h.isDecoy()) continue;
      ++short_targets;
      const double class_wide = q_xl.at(pi.getSpectrumReference());
      if (h.getScore() > class_wide + 1e-12) ++higher;
      if (h.getScore() < class_wide - 1e-12) ++lower;
    }
    TEST_EQUAL(short_targets, 100)
    TEST_EQUAL(lower, 0)
    TEST_EQUAL(higher >= 95, true)
  }
}
END_SECTION

END_TEST
