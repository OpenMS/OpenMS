// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

namespace
{
  /// One ProteinIdentification with a single indistinguishable group over @p paths MS runs.
  ProteinIdentification makeIdOnlyRun(const StringList& paths)
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("PI_0");
    prot_id.setScoreType("q-value");
    prot_id.setHigherScoreBetter(false);
    prot_id.setPrimaryMSRunPath(paths);

    ProteinHit ph1;
    ph1.setAccession("PROT_A");
    ph1.setScore(0.01);
    ph1.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
    ProteinHit ph2;
    ph2.setAccession("PROT_B");
    ph2.setScore(0.02);
    ph2.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
    prot_id.setHits({ph1, ph2});

    ProteinIdentification::ProteinGroup group;
    group.accessions = {"PROT_A", "PROT_B"};
    group.probability = 0.01;
    prot_id.insertIndistinguishableProteins(group);
    return prot_id;
  }

  PeptideIdentificationList makePeptides()
  {
    PeptideIdentification pep_id;
    pep_id.setIdentifier("PI_0");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDEK"));
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pep_id.setHits({hit});
    return PeptideIdentificationList{pep_id};
  }
}

START_TEST(ProteinGroupArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::shared_ptr<arrow::Table> exportToArrow(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  // Single MS run: run_file_name is the origin file without path or extension, matching what the
  // psm and feature tables emit for the same run.
  auto prot_id = makeIdOnlyRun({"/data/SimpleSearchEngine_1.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "SimpleSearchEngine_1")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a merged run has no single origin file per protein group))
{
  // A protein group inferred across several runs has no per-row origin file, and the pg table has
  // no place to record more than one. Emit the file's "unknown" convention ("") plus a warning
  // rather than stamping every row with the run's first file.
  // TODO(#9817): provisional - the QPX spec requires a non-null run_file_name here and quantms
  // drops empty-valued rows, so this may become a non-empty token instead. See the matching
  // comment in ProteinGroupArrowExport.cpp; this assertion pins current behaviour, not a contract.
  auto prot_id = makeIdOnlyRun({"/data/runA.mzML", "/data/runB.mzML"});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "")
}
END_SECTION

START_SECTION(([EXTRA] exportToArrow - a run without spectra_data yields an empty run_file_name))
{
  auto prot_id = makeIdOnlyRun({});
  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, makePeptides());
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
