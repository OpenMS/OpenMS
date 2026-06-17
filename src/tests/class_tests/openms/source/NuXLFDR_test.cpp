// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/NUXL/NuXLFDR.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <cmath>
#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NuXLFDR, "$Id$")

/////////////////////////////////////////////////////////////

NuXLFDR* ptr = nullptr;
NuXLFDR* nullPointer = nullptr;

START_SECTION((NuXLFDR(size_t report_top_hits)))
{
  ptr = new NuXLFDR(2);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~NuXLFDR()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void splitIntoPeptidesAndXLs(const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep_pi, PeptideIdentificationList& xl_pi) const))
{
  // issue #9488, ANAL-42: splitIntoPeptidesAndXLs must keep the first hit of EACH
  // class (plain peptide and cross-link) - not just the first hit overall.

  NuXLFDR fdr(2);

  // --- Case 1: plain-peptide hit first, cross-link hit second ---
  {
    PeptideHit plain;
    plain.setSequence(AASequence::fromString("PEPTIDE"));
    plain.setScore(20.0);
    plain.setMetaValue("NuXL:isXL", 0);

    PeptideHit xl;
    xl.setSequence(AASequence::fromString("PEPTIDER"));
    xl.setScore(10.0);
    xl.setMetaValue("NuXL:isXL", 1);

    PeptideIdentification id;
    id.setHits({plain, xl});

    PeptideIdentificationList ids;
    ids.push_back(id);

    PeptideIdentificationList pep_pi, xl_pi;
    fdr.splitIntoPeptidesAndXLs(ids, pep_pi, xl_pi);

    // both classes must receive the identification (before the fix xl_pi was empty)
    TEST_EQUAL(pep_pi.size(), 1)
    TEST_EQUAL(xl_pi.size(), 1)
    TEST_EQUAL(pep_pi[0].getHits().size(), 1)
    TEST_EQUAL(xl_pi[0].getHits().size(), 1)
    // verify class routing
    TEST_EQUAL((int)pep_pi[0].getHits()[0].getMetaValue("NuXL:isXL"), 0)
    TEST_EQUAL((int)xl_pi[0].getHits()[0].getMetaValue("NuXL:isXL") != 0, true)
  }

  // --- Case 2 (symmetric): cross-link hit first, plain-peptide hit second ---
  {
    PeptideHit xl;
    xl.setSequence(AASequence::fromString("PEPTIDER"));
    xl.setScore(20.0);
    xl.setMetaValue("NuXL:isXL", 1);

    PeptideHit plain;
    plain.setSequence(AASequence::fromString("PEPTIDE"));
    plain.setScore(10.0);
    plain.setMetaValue("NuXL:isXL", 0);

    PeptideIdentification id;
    id.setHits({xl, plain});

    PeptideIdentificationList ids;
    ids.push_back(id);

    PeptideIdentificationList pep_pi, xl_pi;
    fdr.splitIntoPeptidesAndXLs(ids, pep_pi, xl_pi);

    // plain hit must still be captured even though the XL hit ranked first
    TEST_EQUAL(pep_pi.size(), 1)
    TEST_EQUAL(xl_pi.size(), 1)
    TEST_EQUAL(pep_pi[0].getHits().size(), 1)
    TEST_EQUAL(xl_pi[0].getHits().size(), 1)
    TEST_EQUAL((int)pep_pi[0].getHits()[0].getMetaValue("NuXL:isXL"), 0)
    TEST_EQUAL((int)xl_pi[0].getHits()[0].getMetaValue("NuXL:isXL") != 0, true)
  }

  // --- Case 3: only one hit of a single class -> only that list is filled ---
  {
    PeptideHit plain;
    plain.setSequence(AASequence::fromString("PEPTIDE"));
    plain.setScore(20.0);
    plain.setMetaValue("NuXL:isXL", 0);

    PeptideIdentification id;
    id.setHits({plain});

    PeptideIdentificationList ids;
    ids.push_back(id);

    PeptideIdentificationList pep_pi, xl_pi;
    fdr.splitIntoPeptidesAndXLs(ids, pep_pi, xl_pi);

    TEST_EQUAL(pep_pi.size(), 1)
    TEST_EQUAL(xl_pi.size(), 0)
  }
}
END_SECTION

START_SECTION((void calculatePeptideAndXLQValueAndFilterAtPSMLevel(const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep, double peptide_PSM_qvalue_threshold, double peptide_peptide_qvalue_threshold, PeptideIdentificationList& xl_pi, std::vector<double> xl_PSM_qvalue_thresholds, std::vector<double> xl_peptidelevel_qvalue_thresholds, const std::string& out_idxml, int decoy_factor) const))
{
  // issue #9488, ANAL-43: when all XL hits share the same score the tie-breaker
  // divides by a zero score range and (without the guard) writes NaN into every
  // XL hit score. After the fix the scores must remain finite.

  NuXLFDR fdr(1);

  // one (hit-less after pruning) protein identification so tmp_prots[0] is valid
  ProteinIdentification prot;
  prot.setIdentifier("run0");
  ProteinHit ph_prot;
  ph_prot.setAccession("PROT0");
  prot.insertHit(ph_prot);
  vector<ProteinIdentification> protein_ids;
  protein_ids.push_back(prot);

  // two cross-link PSMs with IDENTICAL scores and IDENTICAL "NuXL:score" -> zero range.
  // No peptide evidence is attached so removeUnreferencedProteins prunes the protein hit
  // and the later coverage computation iterates over an empty hit list (safe).
  PeptideIdentificationList peptide_ids;
  for (int i = 0; i < 2; ++i)
  {
    PeptideHit ph;
    ph.setSequence(AASequence::fromString("PEPTIDEK"));
    ph.setScore(10.0);
    ph.setMetaValue("NuXL:isXL", 1);
    ph.setMetaValue("NuXL:score", 5.0); // identical for all hits -> degenerate range
    ph.setMetaValue("target_decoy", "target");

    PeptideIdentification id;
    id.setIdentifier("run0");
    id.setScoreType("NuXL");
    id.setHigherScoreBetter(true);
    id.setSpectrumReference("spec" + std::to_string(i));
    id.setHits({ph});
    peptide_ids.push_back(id);
  }

  PeptideIdentificationList pep_pi, xl_pi;
  std::string tmp_prefix;
  NEW_TMP_FILE(tmp_prefix); // unique writable prefix; the method appends suffixes to it
  std::string out_prefix = tmp_prefix;
  // disabled thresholds (outside the (0,1) open interval)
  vector<double> xl_psm_thresholds = {1.0};
  vector<double> xl_peptide_thresholds = {1.0};

  fdr.calculatePeptideAndXLQValueAndFilterAtPSMLevel(
      protein_ids,
      peptide_ids,
      pep_pi,
      -1.0,   // peptide PSM q-value threshold (disabled)
      -1.0,   // peptide peptide-level q-value threshold (disabled)
      xl_pi,
      xl_psm_thresholds,
      xl_peptide_thresholds,
      out_prefix,
      1);     // decoy_factor

  // the two XL PSMs must survive and carry finite scores (NaN without the fix)
  TEST_EQUAL(xl_pi.size(), 2)
  for (const auto& pi : xl_pi)
  {
    for (const auto& h : pi.getHits())
    {
      TEST_EQUAL(std::isfinite(h.getScore()), true)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
