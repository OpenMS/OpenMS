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
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/type.h>

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
// Helper function to create a test ConsensusMap
/////////////////////////////////////////////////////////////

ConsensusMap createTestConsensusMap()
{
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");

  // Set up column headers
  ConsensusMap::ColumnHeaders headers;
  ConsensusMap::ColumnHeader ch0;
  ch0.filename = "sample1.mzML";
  ch0.label = "";
  headers[0] = ch0;

  ConsensusMap::ColumnHeader ch1;
  ch1.filename = "sample2.mzML";
  ch1.label = "";
  headers[1] = ch1;

  cmap.setColumnHeaders(headers);

  // Set up protein identification
  ProteinIdentification prot_id;
  prot_id.setIdentifier("PI_0");

  ProteinHit ph1;
  ph1.setAccession("PROT_A");
  ph1.setScore(0.99);
  ProteinHit ph2;
  ph2.setAccession("PROT_B");
  ph2.setScore(0.95);
  prot_id.setHits({ph1, ph2});

  // Add protein group
  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.01;  // q-value
  pg.accessions = {"PROT_A", "PROT_B"};
  prot_id.insertProteinGroup(pg);

  cmap.setProteinIdentifications({prot_id});

  // Create ConsensusFeature 1 - identified
  ConsensusFeature cf1;
  cf1.setMZ(500.25);
  cf1.setRT(100.5);
  cf1.setCharge(2);
  cf1.setQuality(0.8);

  // Add feature handles using BaseFeature
  BaseFeature bf1;
  bf1.setIntensity(1000.0f);
  bf1.setMZ(500.25);
  bf1.setRT(100.5);
  cf1.insert(0, bf1);

  BaseFeature bf2;
  bf2.setIntensity(2000.0f);
  bf2.setMZ(500.25);
  bf2.setRT(101.0);
  cf1.insert(1, bf2);

  // Add peptide identification
  PeptideIdentification pep_id1;
  pep_id1.setRT(100.5);
  pep_id1.setMZ(500.25);
  pep_id1.setScoreType("Posterior Error Probability");
  pep_id1.setHigherScoreBetter(false);

  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDE"));
  hit1.setCharge(2);
  hit1.setScore(0.01);
  hit1.setMetaValue("target_decoy", "target");

  PeptideEvidence ev1;
  ev1.setProteinAccession("PROT_A");
  hit1.setPeptideEvidences({ev1});

  pep_id1.setHits({hit1});
  cf1.setPeptideIdentifications({pep_id1});

  cmap.push_back(cf1);

  // Create ConsensusFeature 2 - identified, multiple proteins
  ConsensusFeature cf2;
  cf2.setMZ(600.30);
  cf2.setRT(200.0);
  cf2.setCharge(3);
  cf2.setQuality(0.6);

  BaseFeature bf3;
  bf3.setIntensity(500.0f);
  bf3.setMZ(600.30);
  bf3.setRT(200.0);
  cf2.insert(0, bf3);

  PeptideIdentification pep_id2;
  pep_id2.setRT(200.0);
  pep_id2.setMZ(600.30);
  pep_id2.setScoreType("Posterior Error Probability");
  pep_id2.setHigherScoreBetter(false);

  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setCharge(3);
  hit2.setScore(0.05);
  hit2.setMetaValue("target_decoy", "decoy");

  PeptideEvidence ev2a;
  ev2a.setProteinAccession("PROT_A");
  PeptideEvidence ev2b;
  ev2b.setProteinAccession("PROT_B");
  hit2.setPeptideEvidences({ev2a, ev2b});

  pep_id2.setHits({hit2});
  cf2.setPeptideIdentifications({pep_id2});

  cmap.push_back(cf2);

  // Create ConsensusFeature 3 - unidentified
  ConsensusFeature cf3;
  cf3.setMZ(400.0);
  cf3.setRT(50.0);
  cf3.setCharge(1);
  cf3.setQuality(0.3);

  BaseFeature bf4;
  bf4.setIntensity(100.0f);
  bf4.setMZ(400.0);
  bf4.setRT(50.0);
  cf3.insert(0, bf4);

  cmap.push_back(cf3);

  return cmap;
}

ConsensusMap createEmptyConsensusMap()
{
  return ConsensusMap();
}

/////////////////////////////////////////////////////////////

START_TEST(ConsensusMapArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(exportToArrow - empty consensus map)
{
  ConsensusMap cmap = createEmptyConsensusMap();

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  // Schema should still be present
  TEST_EQUAL(table->num_columns() > 0, true)
}
END_SECTION

START_SECTION(exportToArrow - basic export)
{
  ConsensusMap cmap = createTestConsensusMap();

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 3)

  // Check schema has expected columns
  auto schema = table->schema();
  TEST_EQUAL(schema->GetFieldIndex("sequence") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("peptidoform") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("observed_mz") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("rt") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("charge") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("intensities") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("pg_accessions") >= 0, true)
}
END_SECTION

START_SECTION(exportToArrow - column types)
{
  ConsensusMap cmap = createTestConsensusMap();

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  auto schema = table->schema();

  // Check string columns
  auto seq_field = schema->GetFieldByName("sequence");
  TEST_EQUAL(seq_field->type()->id(), arrow::Type::STRING)

  auto peptidoform_field = schema->GetFieldByName("peptidoform");
  TEST_EQUAL(peptidoform_field->type()->id(), arrow::Type::STRING)

  // Check numeric columns
  auto rt_field = schema->GetFieldByName("rt");
  TEST_EQUAL(rt_field->type()->id(), arrow::Type::FLOAT)

  auto mz_field = schema->GetFieldByName("observed_mz");
  TEST_EQUAL(mz_field->type()->id(), arrow::Type::FLOAT)

  auto charge_field = schema->GetFieldByName("charge");
  TEST_EQUAL(charge_field->type()->id(), arrow::Type::INT16)

  // Check list columns
  auto intensities_field = schema->GetFieldByName("intensities");
  TEST_EQUAL(intensities_field->type()->id(), arrow::Type::LIST)

  auto pg_acc_field = schema->GetFieldByName("pg_accessions");
  TEST_EQUAL(pg_acc_field->type()->id(), arrow::Type::LIST)
}
END_SECTION

START_SECTION(exportToArrow - modifications column)
{
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");

  ConsensusFeature cf;
  cf.setMZ(500.0);
  cf.setRT(100.0);
  cf.setCharge(2);

  PeptideIdentification pep_id;
  pep_id.setScoreType("score");

  PeptideHit hit;
  // Create a modified sequence
  hit.setSequence(AASequence::fromString("PEPTM(Oxidation)IDE"));
  hit.setCharge(2);
  pep_id.setHits({hit});
  cf.setPeptideIdentifications({pep_id});

  cmap.push_back(cf);

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)

  auto schema = table->schema();
  auto mod_field = schema->GetFieldByName("modifications");
  TEST_EQUAL(mod_field->type()->id(), arrow::Type::LIST)
}
END_SECTION

START_SECTION(exportToParquet - basic export)
{
  ConsensusMap cmap = createTestConsensusMap();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
  TEST_EQUAL(File::empty(filename), false)
}
END_SECTION

START_SECTION(exportToParquet - empty consensus map)
{
  ConsensusMap cmap;

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename);
  TEST_EQUAL(success, true)

  // Verify file was created (even for empty data)
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportToParquet - with compression options)
{
  ConsensusMap cmap = createTestConsensusMap();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  // Test with different compression settings
  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::SNAPPY;

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename, pq_config);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportToParquet - ZSTD compression)
{
  ConsensusMap cmap = createTestConsensusMap();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::ZSTD;
  pq_config.compression_level = 9;

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename, pq_config);
  TEST_EQUAL(success, true)

  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportToParquet - no compression)
{
  ConsensusMap cmap = createTestConsensusMap();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::NONE;

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename, pq_config);
  TEST_EQUAL(success, true)

  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

END_TEST
