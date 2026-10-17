// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/test_config.h>
#include <filesystem>
#include <fstream>

using namespace OpenMS;
using ID = IdentificationData;

START_TEST(IdentificationDataFileHandler, "$Id$")

START_SECTION((owning FileHandler writes native bundles and recognizes their manifest without an extension))
{
  ID data;
  auto& run = data.addRun("run");
  ID::SourceFile source;
  source.path = "input.mzML";
  const auto sid = run.addSource(source);
  ID::ScoreDefinition definition;
  definition.name = "search score";
  const auto score = run.addScore(definition);
  const auto query = run.addIdentification(sid, {});
  ID::MatchData peptide;
  peptide.representation = "PEPTIDE";
  peptide.charge = 2;
  const auto match = run.addMatch(query, peptide, {9.0});
  run.setPrimaryScore(score);
  const auto uuid = run.getUuid();
  std::string path;
  NEW_TMP_FILE(path)
  FileHandler().storeIdentifications(path, data, {FileTypes::IDPARQUET});
  TEST_TRUE(IdentificationDataFile::isNativeFile(path))
  ID loaded;
  FileHandler().loadIdentifications(path, loaded);
  TEST_EQUAL(loaded.getRuns().size(), 1)
  TEST_EQUAL(loaded.getRun("run").getUuid(), uuid)
  TEST_REAL_SIMILAR(*loaded.getRun("run").getScore(match, loaded.getRun("run").getScoreId(0)), 9.0)
  TEST_EXCEPTION(Exception::InvalidFileType, FileHandler().loadIdentifications(path, loaded, {FileTypes::IDXML}))
  TEST_EQUAL(loaded.getRun("run").getUuid(), uuid)
  // The manifest identifies the bundle without an extension.
  TEST_EQUAL(FileHandler::getType(path), FileTypes::IDPARQUET)
  // A tool rerun replaces its previous native bundle, like any other output format.
  peptide.representation = "SEQVENCE";
  data.getRun("run").addMatch(query, peptide, {7.0});
  FileHandler().storeIdentifications(path, data, {FileTypes::IDPARQUET});
  FileHandler().loadIdentifications(path, loaded);
  TEST_EQUAL(loaded.getRun("run").getNumberOfMatches(), 2)
  TEST_EQUAL(loaded.getRun("run").getUuid(), uuid)
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  FileHandler().loadIdentifications(path, proteins, peptides, {FileTypes::IDPARQUET});
  TEST_EQUAL(peptides.size(), 1)
  TEST_EQUAL(peptides.front().getHits().front().getSequence().toString(), "PEPTIDE")
  TEST_REAL_SIMILAR(peptides.front().getHits().front().getScore(), 9.0)
  // .idparquet is native only: established identifications are imported and replace the bundle, like a tool rerun.
  FileHandler().storeIdentifications(path, proteins, peptides, {FileTypes::IDPARQUET});
  TEST_TRUE(IdentificationDataFile::isNativeFile(path))
  TEST_FALSE(std::filesystem::exists(std::filesystem::path(path) / "psms.parquet"))
  FileHandler().loadIdentifications(path, loaded);
  ABORT_IF(loaded.getRuns().size() != 1)
  TEST_EQUAL(loaded.getRuns()[0].getNumberOfMatches(), 2)
  // The import is a new run (the legacy load named it like IdXMLFile does).
  TEST_NOT_EQUAL(loaded.getRuns()[0].getUuid(), uuid)
  std::filesystem::remove_all(path);

  // The four-table layout of OpenMS 3.6 is not read, by either overload.
  const auto legacy_bundle = path + ".idparquet";
  std::filesystem::create_directories(legacy_bundle);
  std::ofstream(std::filesystem::path(legacy_bundle) / "psms.parquet") << "PAR1";
  TEST_FALSE(IdentificationDataFile::isNativeFile(legacy_bundle))
  TEST_EXCEPTION(Exception::InvalidFileType, FileHandler().loadIdentifications(legacy_bundle, proteins, peptides, {FileTypes::IDPARQUET}))
  TEST_EXCEPTION(Exception::InvalidFileType, FileHandler().loadIdentifications(legacy_bundle, loaded))
  std::filesystem::remove_all(legacy_bundle);
}
END_SECTION

START_SECTION((established idXML input imports into the owning model))
{
  ID data;
  TEST_EXCEPTION(Exception::InvalidValue,
                 FileHandler().loadIdentifications(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), data, {FileTypes::IDXML}))
  TEST_TRUE(data.getRuns().empty())
  FileHandler().loadIdentifications(OPENMS_GET_TEST_DATA_PATH("IDScoreSwitcherAlgorithm_test_input.idXML"), data, {FileTypes::IDXML});
  TEST_FALSE(data.getRuns().empty())
  Size matches = 0;
  for (const auto& run : data.getRuns())
    matches += run.getNumberOfMatches();
  TEST_NOT_EQUAL(matches, 0)
}
END_SECTION

END_TEST
