// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/IdentificationSummary.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <vector>

///////////////////////////

START_TEST(IdentificationSummary, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// Helper: build a PeptideIdentification with one hit of the given sequence.
auto mkPep = [](const std::string& seq, const std::string& score_type, float sig) -> PeptideIdentification
{
  PeptideIdentification pid;
  PeptideHit h;
  h.setSequence(AASequence::fromString(seq));
  pid.setHits({h});
  pid.setScoreType(score_type);
  pid.setSignificanceThreshold(sig);
  return pid;
};

// Helper: build a single ProteinIdentification with trypsin digestion and two hits.
auto mkProtIds = [](const std::string& score_type, float sig) -> std::vector<ProteinIdentification>
{
  std::vector<ProteinIdentification> prot_ids(1);
  prot_ids[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  prot_ids[0].getSearchParameters().missed_cleavages = 4;
  prot_ids[0].setScoreType(score_type);
  prot_ids[0].setSignificanceThreshold(sig);
  std::vector<ProteinHit> phs;
  ProteinHit p1; p1.setAccession("PROT1"); p1.setScore(0.9); phs.push_back(p1);
  ProteinHit p2; p2.setAccession("PROT2"); p2.setScore(0.5); phs.push_back(p2);
  prot_ids[0].setHits(phs);
  return prot_ids;
};

IdentificationSummary* ptr = nullptr;
IdentificationSummary* null_ptr = nullptr;

START_SECTION((IdentificationSummary()))
{
  ptr = new IdentificationSummary();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~IdentificationSummary()))
{
  delete ptr;
}
END_SECTION

START_SECTION((Result compute(std::vector<ProteinIdentification>& prot_ids, PeptideIdentificationList& pep_ids)))
{
  // FDR score type -> fdr thresholds reported
  {
    std::vector<ProteinIdentification> prot_ids = mkProtIds("FDR", 0.01f);
    PeptideIdentificationList pep_ids;
    pep_ids.push_back(mkPep("AAAAAAAAAAAAAK", "FDR", 0.05f));  // length 14, 0 missed cleavages
    pep_ids.push_back(mkPep("AAAAKAAAAAR", "FDR", 0.05f));     // length 11, 1 missed cleavage
    pep_ids.push_back(mkPep("AAAKAAARAAAKAAAR", "FDR", 0.05f));// length 16, 3 missed cleavages

    IdentificationSummary is;
    IdentificationSummary::Result r = is.compute(prot_ids, pep_ids);

    TEST_EQUAL(r.peptide_spectrum_matches, 3)
    TEST_EQUAL(r.unique_peptides.count, 3)
    TEST_EQUAL(r.unique_proteins.count, 2)
    TEST_REAL_SIMILAR(r.unique_peptides.fdr_threshold, 0.05)
    TEST_REAL_SIMILAR(r.unique_proteins.fdr_threshold, 0.01)
    // (0 + 1 + 3) missed cleavages over 3 peptides
    TEST_REAL_SIMILAR(r.missed_cleavages_mean, 4.0 / 3.0)
    // (0.9 + 0.5) / 2 protein hits
    TEST_REAL_SIMILAR(r.protein_hit_scores_mean, 0.7)
    // (14 + 11 + 16) / 3 unique peptides
    TEST_REAL_SIMILAR(r.peptide_length_mean, 41.0 / 3.0)
  }

  // Non-FDR score type -> fdr thresholds stay at the -1 sentinel
  {
    std::vector<ProteinIdentification> prot_ids = mkProtIds("Posterior Error Probability", 0.0f);
    PeptideIdentificationList pep_ids;
    pep_ids.push_back(mkPep("AAAAAAAAAAAAAK", "Posterior Error Probability", 0.0f));
    pep_ids.push_back(mkPep("AAAAKAAAAAR", "Posterior Error Probability", 0.0f));

    IdentificationSummary is;
    IdentificationSummary::Result r = is.compute(prot_ids, pep_ids);

    TEST_EQUAL(r.peptide_spectrum_matches, 2)
    TEST_EQUAL(r.unique_peptides.count, 2)
    TEST_REAL_SIMILAR(r.unique_peptides.fdr_threshold, -1.0)
    TEST_REAL_SIMILAR(r.unique_proteins.fdr_threshold, -1.0)
  }
}
END_SECTION

START_SECTION(([IdentificationSummary::UniqueID] defaults))
{
  IdentificationSummary::UniqueID u;
  TEST_EQUAL(u.count, 0)
  TEST_REAL_SIMILAR(u.fdr_threshold, -1.0)
}
END_SECTION

START_SECTION(([IdentificationSummary::Result] bool operator==(const Result& rhs) const))
{
  IdentificationSummary::Result a, b;
  TEST_EQUAL(a == b, true)

  a.peptide_spectrum_matches = 7;
  TEST_EQUAL(a == b, false)
  b.peptide_spectrum_matches = 7;
  TEST_EQUAL(a == b, true)

  a.unique_peptides.count = 3;
  TEST_EQUAL(a == b, false)
  b.unique_peptides.count = 3;
  TEST_EQUAL(a == b, true)

  a.protein_hit_scores_mean = 0.7;
  TEST_EQUAL(a == b, false)
  b.protein_hit_scores_mean = 0.7;
  TEST_EQUAL(a == b, true)
}
END_SECTION

START_SECTION((const std::string& getName() const override))
{
  IdentificationSummary is;
  TEST_STRING_EQUAL(is.getName(), "Summary of detected Proteins and Peptides from idXML file")
}
END_SECTION

START_SECTION((QCBase::Status requirements() const override))
{
  IdentificationSummary is;
  TEST_EQUAL(is.requirements() == QCBase::Status(QCBase::Requires::ID), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
