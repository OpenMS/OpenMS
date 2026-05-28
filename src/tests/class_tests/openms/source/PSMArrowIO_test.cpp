// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/PSMArrowIO.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

namespace
{
  void buildMinimalIds(std::vector<ProteinIdentification>& prot_ids,
                      PeptideIdentificationList& pep_ids)
  {
    ProteinIdentification prot;
    prot.setIdentifier("run_1");
    prot.setScoreType("score");
    prot.setHigherScoreBetter(true);
    prot.setPrimaryMSRunPath({"/tmp/foo.mzML"});
    prot.getSearchParameters().digestion_enzyme.setName("Trypsin");
    prot.getSearchParameters().setMetaValue("extra_features", "COMET:deltaCn,MS:1002049");
    prot_ids = {prot};

    PeptideIdentification pid;
    pid.setIdentifier("run_1");
    pid.setScoreType("score");
    pid.setHigherScoreBetter(true);
    pid.setRT(123.4);
    pid.setMZ(567.89);
    pid.setSpectrumReference("scan=42");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    hit.setScore(0.95);
    hit.setMetaValue("target_decoy", "target");
    hit.setMetaValue("COMET:deltaCn", 0.5);
    PeptideEvidence ev;
    ev.setProteinAccession("sp|P12345|EXAMPLE");
    ev.setAABefore('K');
    ev.setAAAfter('R');
    ev.setStart(42);
    ev.setEnd(48);
    hit.addPeptideEvidence(ev);
    pid.getHits().push_back(hit);
    pep_ids.push_back(pid);
  }
}

START_TEST(PSMArrowIO, "$Id$")

START_SECTION(([EXTRA] export_then_import_round_trip))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  TEST_TRUE(File::exists(dir + "/psms.parquet"));
  TEST_TRUE(File::exists(dir + "/proteins.parquet"));
  TEST_TRUE(File::exists(dir + "/protein_groups.parquet"));
  TEST_TRUE(File::exists(dir + "/search_params.parquet"));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  TEST_EQUAL(prot_ids_in.size(), 1);
  // Identifier is synthesized on load per IdXMLFile.cpp:530 parity — the stored
  // "run_1" is replaced with `<search_engine>_<date>_<UniqueIdGenerator>`. The
  // pep_id collection is re-stamped in lock-step so both sides agree.
  TEST_NOT_EQUAL(prot_ids_in[0].getIdentifier(), "");
  TEST_NOT_EQUAL(prot_ids_in[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(pep_ids_in[0].getIdentifier(), prot_ids_in[0].getIdentifier());
  TEST_STRING_EQUAL(prot_ids_in[0].getSearchParameters().digestion_enzyme.getName(), "Trypsin");
  TEST_STRING_EQUAL(String(prot_ids_in[0].getSearchParameters().getMetaValue("extra_features")),
                    "COMET:deltaCn,MS:1002049");

  // Primary MS run path round-trips (spectra_data UserParam).
  StringList msrun_in;
  prot_ids_in[0].getPrimaryMSRunPath(msrun_in);
  TEST_EQUAL(msrun_in.size(), 1);
  TEST_STRING_EQUAL(msrun_in.front(), "/tmp/foo.mzML");

  TEST_EQUAL(pep_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in[0].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_in[0].getHits()[0].getSequence().toString(), "PEPTIDE");
  TEST_REAL_SIMILAR(double(pep_ids_in[0].getHits()[0].getMetaValue("COMET:deltaCn")), 0.5);

  // PeptideEvidence positional fields round-trip (.idparquet schema extension).
  TEST_EQUAL(pep_ids_in[0].getHits()[0].getPeptideEvidences().size(), 1);
  const auto& ev_in = pep_ids_in[0].getHits()[0].getPeptideEvidences().front();
  TEST_STRING_EQUAL(ev_in.getProteinAccession(), "sp|P12345|EXAMPLE");
  TEST_EQUAL(ev_in.getAABefore(), 'K');
  TEST_EQUAL(ev_in.getAAAfter(),  'R');
  TEST_EQUAL(ev_in.getStart(), 42);
  TEST_EQUAL(ev_in.getEnd(),   48);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] import_missing_subfile_returns_false))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";
  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  // Delete one subfile.
  File::remove(dir + "/psms.parquet");

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_FALSE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] export_target_is_regular_file_returns_false))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String path;
  NEW_TMP_FILE(path)
  path += ".idparquet";
  // Create a regular file at the target path.
  std::ofstream f(path); f << "not a directory"; f.close();

  TEST_FALSE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, path));

  File::remove(path);
}
END_SECTION

START_SECTION(([EXTRA] empty_psms_round_trips))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);
  pep_ids.clear();

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));
  TEST_EQUAL(pep_ids_in.size(), 0);
  TEST_EQUAL(prot_ids_in.size(), 1);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] decoy_round_trips))
{
  // Build a single decoy PSM (target_decoy=decoy) and verify is_decoy round-trips.
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;

  ProteinIdentification prot;
  prot.setIdentifier("run_decoy");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  prot_ids.push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_decoy");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("DECOYSEQ"));
  hit.setCharge(2);
  hit.setScore(0.1);
  hit.setMetaValue("target_decoy", "decoy");
  PeptideEvidence ev;
  ev.setProteinAccession("DECOY_sp|P99999|FAKE");
  hit.addPeptideEvidence(ev);
  pid.getHits().push_back(hit);
  pep_ids.push_back(pid);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  TEST_EQUAL(pep_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in[0].getHits().size(), 1);
  TEST_STRING_EQUAL(String(pep_ids_in[0].getHits()[0].getMetaValue("target_decoy")), "decoy");

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] multi_rank_round_trips))
{
  // Build a PID with three rank-ordered PSMs and verify all three survive round-trip
  // (this is the new default behaviour after the keep_all_passing -> best_per_spectrum_only flip).
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;

  ProteinIdentification prot;
  prot.setIdentifier("run_multirank");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  prot_ids.push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_multirank");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);

  for (int rank = 1; rank <= 3; ++rank)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(std::string(rank, 'P') + "EPTIDE"));
    hit.setCharge(2);
    hit.setScore(1.0 / rank);
    hit.setRank(static_cast<UInt>(rank));
    PeptideEvidence ev;
    ev.setProteinAccession("sp|P" + String(rank) + "|RANK");
    hit.addPeptideEvidence(ev);
    pid.getHits().push_back(hit);
  }
  pep_ids.push_back(pid);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir, /*export_all_psms=*/true));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  TEST_EQUAL(pep_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in[0].getHits().size(), 3);
  for (size_t i = 0; i < 3; ++i)
  {
    const auto& h = pep_ids_in[0].getHits()[i];
    TEST_EQUAL(h.getRank(), static_cast<UInt>(i + 1));
  }
  TEST_REAL_SIMILAR(pep_ids_in[0].getHits()[0].getScore(), 1.0);
  TEST_REAL_SIMILAR(pep_ids_in[0].getHits()[1].getScore(), 0.5);
  TEST_REAL_SIMILAR(pep_ids_in[0].getHits()[2].getScore(), 1.0/3.0);

  File::removeDirRecursively(dir);
}
END_SECTION

END_TEST
