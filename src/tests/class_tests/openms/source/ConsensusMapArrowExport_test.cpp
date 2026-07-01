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
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

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

  std::string filename;
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

  std::string filename;
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

  std::string filename;
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

  std::string filename;
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

  std::string filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::NONE;

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename, pq_config);
  TEST_EQUAL(success, true)

  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION((static bool exportToParquetStreaming(const ConsensusMap& cmap, const std::string& filename, size_t batch_size, const ParquetWriteConfig& config)))
{
  arrow::MemoryPool* pool = arrow::default_memory_pool();

  // Read a parquet file into a single-chunk (combined) Arrow table. Streaming writes
  // one row group per batch, so columns come back multi-chunk; combine before row indexing.
  auto read_combined = [pool](const std::string& path) -> std::shared_ptr<arrow::Table>
  {
    auto open_res = arrow::io::ReadableFile::Open(path.c_str());
    if (!open_res.ok()) { return nullptr; }
    auto reader_res = parquet::arrow::OpenFile(open_res.ValueOrDie(), pool);
    if (!reader_res.ok()) { return nullptr; }
    auto reader = std::move(reader_res).ValueOrDie();
    std::shared_ptr<arrow::Table> t;
    if (!reader->ReadTable(&t).ok()) { return nullptr; }
    auto comb = t->CombineChunks(pool);
    if (!comb.ok()) { return nullptr; }
    return comb.ValueOrDie();
  };

  // Build a ConsensusMap with many features and a UNIQUE rt per feature so the rt
  // check below is order-sensitive (catches any batch-boundary reordering). Mixes
  // identified/unidentified, modified/unmodified, target/decoy.
  ConsensusMap cmap;
  ConsensusMap::ColumnHeaders headers;
  ConsensusMap::ColumnHeader ch0; ch0.filename = "sample1.mzML"; headers[0] = ch0;
  cmap.setColumnHeaders(headers);

  ProteinIdentification prot_id;
  prot_id.setIdentifier("PI_0");
  ProteinHit ph; ph.setAccession("PROT_A"); ph.setScore(0.99);
  prot_id.setHits({ph});
  ProteinIdentification::ProteinGroup pg; pg.probability = 0.01; pg.accessions = {"PROT_A"};
  prot_id.insertProteinGroup(pg);
  cmap.setProteinIdentifications({prot_id});

  const size_t N = 2500;
  for (size_t i = 0; i < N; ++i)
  {
    ConsensusFeature cf;
    cf.setMZ(400.0 + (i % 500));
    cf.setRT(100.0 + i);          // unique per row -> order-sensitive
    cf.setCharge(2 + (i % 3));
    BaseFeature bf; bf.setIntensity(1000.0f + i); bf.setMZ(cf.getMZ()); bf.setRT(cf.getRT());
    cf.insert(0, bf);
    if (i % 4 != 0) // most identified, every 4th left unidentified
    {
      PeptideIdentification pid;
      pid.setRT(cf.getRT()); pid.setMZ(cf.getMZ());
      pid.setScoreType("Posterior Error Probability"); pid.setHigherScoreBetter(false);
      PeptideHit hit;
      hit.setSequence(AASequence::fromString((i % 3 == 0) ? "PEM(Oxidation)TIDER" : "PEPTIDEK"));
      hit.setCharge(cf.getCharge());
      hit.setScore(0.01 + (i % 100) * 0.001);
      hit.setMetaValue("target_decoy", (i % 2 == 0) ? "target" : "decoy");
      PeptideEvidence ev; ev.setProteinAccession("PROT_A");
      hit.setPeptideEvidences({ev});
      pid.setHits({hit});
      cf.setPeptideIdentifications({pid});
    }
    cmap.push_back(cf);
  }

  // --- Streaming write with a small batch (multiple row groups: 2500 / 1000) ---
  std::string stream_file; NEW_TMP_FILE(stream_file)
  TEST_EQUAL(ConsensusMapArrowExport::exportToParquetStreaming(cmap, stream_file, 1000), true)
  auto st_table = read_combined(stream_file);
  TEST_NOT_EQUAL(st_table, nullptr)
  TEST_EQUAL(st_table->num_rows(), (int64_t)N)

  // --- One-shot reference ---
  std::string ref_file; NEW_TMP_FILE(ref_file)
  TEST_EQUAL(ConsensusMapArrowExport::exportToParquet(cmap, ref_file), true)
  auto ref_table = read_combined(ref_file);
  TEST_NOT_EQUAL(ref_table, nullptr)
  TEST_EQUAL(ref_table->num_rows(), st_table->num_rows())
  TEST_EQUAL(st_table->num_columns(), ref_table->num_columns())

  // --- Column-wise, order-sensitive equivalence on representative columns ---
  {
    auto seq_s = std::static_pointer_cast<arrow::StringArray>(st_table->GetColumnByName("sequence")->chunk(0));
    auto seq_r = std::static_pointer_cast<arrow::StringArray>(ref_table->GetColumnByName("sequence")->chunk(0));
    auto pf_s  = std::static_pointer_cast<arrow::StringArray>(st_table->GetColumnByName("peptidoform")->chunk(0));
    auto pf_r  = std::static_pointer_cast<arrow::StringArray>(ref_table->GetColumnByName("peptidoform")->chunk(0));
    auto chg_s = std::static_pointer_cast<arrow::Int16Array>(st_table->GetColumnByName("charge")->chunk(0));
    auto chg_r = std::static_pointer_cast<arrow::Int16Array>(ref_table->GetColumnByName("charge")->chunk(0));
    auto rt_s  = std::static_pointer_cast<arrow::FloatArray>(st_table->GetColumnByName("rt")->chunk(0));
    auto rt_r  = std::static_pointer_cast<arrow::FloatArray>(ref_table->GetColumnByName("rt")->chunk(0));
    auto dec_s = std::static_pointer_cast<arrow::BooleanArray>(st_table->GetColumnByName("is_decoy")->chunk(0));
    auto dec_r = std::static_pointer_cast<arrow::BooleanArray>(ref_table->GetColumnByName("is_decoy")->chunk(0));
    bool all_eq = true;
    for (int64_t r = 0; r < st_table->num_rows(); ++r)
    {
      if (seq_s->GetString(r) != seq_r->GetString(r)) { all_eq = false; break; }
      if (pf_s->GetString(r)  != pf_r->GetString(r))  { all_eq = false; break; }
      if (chg_s->Value(r)     != chg_r->Value(r))     { all_eq = false; break; }
      if (rt_s->Value(r)      != rt_r->Value(r))      { all_eq = false; break; } // order-sensitive
      if (dec_s->Value(r)     != dec_r->Value(r))     { all_eq = false; break; }
    }
    TEST_EQUAL(all_eq, true)
  }

  // --- Full logical equivalence: also covers nested/list columns (modifications,
  // intensities, pg_accessions, gg_*) and lookup-derived columns (pg_global_qvalue).
  // Both tables were CombineChunks()ed and round-tripped through parquet the same way,
  // so schema (incl. nested types) and values must match; ignore schema metadata. ---
  TEST_EQUAL(st_table->Equals(*ref_table), true)

  // --- batch_size >= N: single batch, still correct ---
  {
    std::string one_batch; NEW_TMP_FILE(one_batch)
    TEST_EQUAL(ConsensusMapArrowExport::exportToParquetStreaming(cmap, one_batch, N + 1000), true)
    auto t = read_combined(one_batch);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), (int64_t)N)
  }

  // --- batch_size == 0 guard: treated as default, valid file ---
  {
    std::string zero_batch; NEW_TMP_FILE(zero_batch)
    TEST_EQUAL(ConsensusMapArrowExport::exportToParquetStreaming(cmap, zero_batch, 0), true)
    auto t = read_combined(zero_batch);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), (int64_t)N)
  }

  // --- Edge case: empty map -> valid 0-row file with full schema ---
  {
    ConsensusMap empty;
    std::string empty_file; NEW_TMP_FILE(empty_file)
    TEST_EQUAL(ConsensusMapArrowExport::exportToParquetStreaming(empty, empty_file, 1000), true)
    TEST_EQUAL(File::exists(empty_file), true)
    auto t = read_combined(empty_file);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), 0)
    TEST_EQUAL(t->num_columns() > 0, true)
  }
}
END_SECTION

END_TEST
