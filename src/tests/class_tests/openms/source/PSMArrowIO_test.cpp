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
  TEST_STRING_EQUAL(prot_ids_in[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(prot_ids_in[0].getSearchParameters().digestion_enzyme.getName(), "Trypsin");
  TEST_STRING_EQUAL(String(prot_ids_in[0].getSearchParameters().getMetaValue("extra_features")),
                    "COMET:deltaCn,MS:1002049");

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

END_TEST
