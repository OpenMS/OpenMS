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
#include <OpenMS/FORMAT/QPXFile.h>
///////////////////////////

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

using namespace OpenMS;
using namespace std;

// Helpers for the merged-run sections below.
namespace
{
  std::shared_ptr<arrow::Table> readQPXParquet(const std::string& path)
  {
    arrow::MemoryPool* pool = arrow::default_memory_pool();
    auto o = arrow::io::ReadableFile::Open(path.c_str());
    if (!o.ok()) { return nullptr; }
    auto rr = parquet::arrow::OpenFile(o.ValueOrDie(), pool);
    if (!rr.ok()) { return nullptr; }
    auto rd = std::move(rr).ValueOrDie();
    std::shared_ptr<arrow::Table> t;
    if (!rd->ReadTable(&t).ok()) { return nullptr; }
    auto c = t->CombineChunks(pool);
    return c.ok() ? c.ValueOrDie() : nullptr;
  }

  /// One ProteinIdentification whose 'spectra_data' lists @p paths - the shape a merged run has.
  std::vector<ProteinIdentification> mergedRun(const StringList& paths)
  {
    std::vector<ProteinIdentification> prot(1);
    prot[0].setIdentifier("merged");
    prot[0].setScoreType("q-value");
    prot[0].setHigherScoreBetter(false);
    prot[0].setPrimaryMSRunPath(paths);
    return prot;
  }

  /// A PSM belonging to run "merged". Pass a negative @p merge_index to omit the metavalue.
  PeptideIdentification mergedPSM(const std::string& sequence, int merge_index, bool set_index = true)
  {
    PeptideIdentification pid;
    pid.setIdentifier("merged");
    pid.setRT(42.0);
    pid.setMZ(555.5);
    pid.setScoreType("q-value");
    pid.setHigherScoreBetter(false);
    pid.setSpectrumReference("controllerType=0 controllerNumber=1 scan=17");
    if (set_index) { pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, merge_index); }
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setCharge(2);
    hit.setScore(0.001);
    pid.setHits(vector<PeptideHit>{hit});
    return pid;
  }

  std::vector<const PeptideIdentification*> ptrsOf(const PeptideIdentificationList& peps)
  {
    std::vector<const PeptideIdentification*> ptrs;
    ptrs.reserve(peps.size());
    for (const auto& p : peps) { ptrs.push_back(&p); }
    return ptrs;
  }
}

START_TEST(QPXFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(...))
{
  // Create test data
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  std::vector<std::string> pep_strs = {"PEPTIDER", "PEM(Oxidation)TIDER", "DFPIANGER"};
  for (size_t i = 0; i < pep_strs.size(); ++i)
  {
    PeptideIdentification peptide_id;
    peptide_id.setIdentifier("test_search");
    peptide_id.setRT(1234.5 + i * 100);
    peptide_id.setMZ(500.25 + i * 50);
    peptide_id.setScoreType("TestScore");
    peptide_id.setSpectrumReference("controllerType=0 controllerNumber=1 scan=" + StringUtils::toStr(1000 + i));

    PeptideHit hit;
    hit.setSequence(AASequence::fromString(pep_strs[i]));
    hit.setScore(0.95 - i * 0.1);
    hit.setCharge(2 + i);

    if (i % 2 == 0)
    {
      hit.setMetaValue("target_decoy", "target");
    }
    else
    {
      hit.setMetaValue("target_decoy", "decoy");
    }

    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + StringUtils::toStr(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});

    // Add multiple hits, but only first should be processed by default
    vector<PeptideHit> hits = {hit};
    if (i == 0)
    {
      PeptideHit second_hit;
      second_hit.setSequence(AASequence::fromString("SECOND"));
      second_hit.setScore(0.5);
      second_hit.setCharge(1);
      hits.push_back(second_hit);
    }

    peptide_id.setHits(hits);
    peptide_ids.push_back(peptide_id);
  }

  // Export to Arrow table (default: best hit only)
  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (should equal number of peptide identifications, not hits)
  TEST_EQUAL(table->num_rows(), 3)

  // Verify schema column names and count (24 columns in QPXPSMSchema)
  auto schema = table->schema();
  TEST_EQUAL(table->num_columns(), 24)

  TEST_EQUAL(schema->field(0)->name(), "sequence")
  TEST_EQUAL(schema->field(1)->name(), "peptidoform")
  TEST_EQUAL(schema->field(2)->name(), "modifications")
  TEST_EQUAL(schema->field(3)->name(), "charge")
  TEST_EQUAL(schema->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(schema->field(5)->name(), "is_decoy")
  TEST_EQUAL(schema->field(6)->name(), "calculated_mz")
  TEST_EQUAL(schema->field(7)->name(), "observed_mz")
  TEST_EQUAL(schema->field(8)->name(), "mass_error_ppm")
  TEST_EQUAL(schema->field(9)->name(), "additional_scores")
  TEST_EQUAL(schema->field(10)->name(), "predicted_rt")
  TEST_EQUAL(schema->field(11)->name(), "run_file_name")
  TEST_EQUAL(schema->field(12)->name(), "cv_params")
  TEST_EQUAL(schema->field(13)->name(), "scan")
  TEST_EQUAL(schema->field(14)->name(), "rt")
  TEST_EQUAL(schema->field(15)->name(), "ion_mobility")
  TEST_EQUAL(schema->field(16)->name(), "missed_cleavages")
  TEST_EQUAL(schema->field(17)->name(), "protein_accessions")
  TEST_EQUAL(schema->field(18)->name(), "cross_links")
  TEST_EQUAL(schema->field(19)->name(), "mz_array")
  TEST_EQUAL(schema->field(20)->name(), "intensity_array")
  TEST_EQUAL(schema->field(21)->name(), "charge_array")
  TEST_EQUAL(schema->field(22)->name(), "ion_type_array")
  TEST_EQUAL(schema->field(23)->name(), "ion_mobility_array")

  // Verify data types for key columns
  TEST_EQUAL(schema->field(3)->type()->id(), arrow::Type::INT16)   // charge is int16
  TEST_EQUAL(schema->field(4)->type()->id(), arrow::Type::DOUBLE)  // PEP is float64
  TEST_EQUAL(schema->field(6)->type()->id(), arrow::Type::FLOAT)   // calculated_mz is float32
  TEST_EQUAL(schema->field(7)->type()->id(), arrow::Type::FLOAT)   // observed_mz is float32
  TEST_EQUAL(schema->field(8)->type()->id(), arrow::Type::FLOAT)   // mass_error_ppm is float32
  TEST_EQUAL(schema->field(13)->type()->id(), arrow::Type::LIST)   // scan is list<int32>
  TEST_EQUAL(schema->field(17)->type()->id(), arrow::Type::LIST)   // protein_accessions is list
  TEST_EQUAL(schema->field(2)->type()->id(), arrow::Type::LIST)    // modifications is list

  // Verify sequence values
  auto seq_col = table->GetColumnByName("sequence");
  auto seq_arr = std::static_pointer_cast<arrow::StringArray>(seq_col->chunk(0));
  TEST_EQUAL(seq_arr->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_arr->GetString(1), "PEMTIDER")
  TEST_EQUAL(seq_arr->GetString(2), "DFPIANGER")

  // Verify peptidoform uses ProForma notation (not bracket substitution)
  auto pf_col = table->GetColumnByName("peptidoform");
  auto pf_arr = std::static_pointer_cast<arrow::StringArray>(pf_col->chunk(0));
  TEST_EQUAL(pf_arr->GetString(0), "PEPTIDER")
  // Modified peptide should have ProForma notation with bracket-enclosed modification
  std::string mod_pf = pf_arr->GetString(1);
  TEST_EQUAL(mod_pf.empty(), false)
  TEST_EQUAL(mod_pf.find('[') != std::string::npos, true) // ProForma bracket notation expected

  // Verify is_decoy non-nullable behavior
  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::BooleanArray>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->Value(0), false) // target
  TEST_EQUAL(decoy_arr->Value(1), true) // decoy
  TEST_EQUAL(decoy_arr->Value(2), false) // target

  // Verify scan extraction from native ID (now list<int32>)
  auto scan_col = table->GetColumnByName("scan");
  auto scan_list = std::static_pointer_cast<arrow::ListArray>(scan_col->chunk(0));
  // Each row should have a list with one scan number
  auto scan_values = std::static_pointer_cast<arrow::Int32Array>(scan_list->values());
  TEST_EQUAL(scan_list->value_length(0), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(0)), 1000)
  TEST_EQUAL(scan_list->value_length(1), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(1)), 1001)
  TEST_EQUAL(scan_list->value_length(2), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(2)), 1002)
}
END_SECTION

START_SECTION(static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(...) with export_all_psms)
{
  // Create test data with multiple hits per peptide identification
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");

  vector<PeptideHit> hits;
  vector<std::string> pep_strs = {"PEPTIDER", "ALTERNATIVE", "THIRDPSM"};
  for (int i = 0; i < 3; ++i)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(pep_strs[i]));
    hit.setScore(0.95 - i * 0.1);
    hit.setCharge(2 + i);
    hit.setMetaValue("target_decoy", "target");

    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + StringUtils::toStr(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});

    hits.push_back(hit);
  }

  peptide_id.setHits(hits);
  peptide_ids.push_back(peptide_id);

  // Export with all PSMs
  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids, true);
  TEST_NOT_EQUAL(table, nullptr)

  // All 3 hits should be exported
  TEST_EQUAL(table->num_rows(), 3)

  // Verify charge column (int16)
  auto charge_col = table->GetColumnByName("charge");
  auto charge_arr = std::static_pointer_cast<arrow::Int16Array>(charge_col->chunk(0));
  TEST_EQUAL(charge_arr->Value(0), 2)
  TEST_EQUAL(charge_arr->Value(1), 3)
  TEST_EQUAL(charge_arr->Value(2), 4)

  // Verify sequence values
  auto seq_col = table->GetColumnByName("sequence");
  auto seq_arr = std::static_pointer_cast<arrow::StringArray>(seq_col->chunk(0));
  TEST_EQUAL(seq_arr->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_arr->GetString(1), "ALTERNATIVE")
  TEST_EQUAL(seq_arr->GetString(2), "THIRDPSM")
}
END_SECTION

START_SECTION(static bool exportToParquet(...))
{
  // Create test data
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEM(Oxidation)TIDER"));
  hit.setScore(0.95);
  hit.setCharge(2);
  hit.setMetaValue("target_decoy", "target");

  PeptideEvidence evidence;
  evidence.setProteinAccession("TEST_PROTEIN");
  PeptideEvidence evidence2;
  evidence2.setProteinAccession("TEST_PROTEIN_2");
  hit.setPeptideEvidences(vector<PeptideEvidence>{evidence, evidence2});

  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);

  std::string output_file;
  NEW_TMP_FILE(output_file)

  // Write parquet
  bool success = QPXFile::exportToParquet(protein_ids, peptide_ids, output_file);
  TEST_EQUAL(success, true)

  // Read back and verify
  arrow::MemoryPool* pool = arrow::default_memory_pool();
  std::shared_ptr<arrow::io::ReadableFile> infile;
  auto result = arrow::io::ReadableFile::Open(output_file.c_str());
  TEST_EQUAL(result.ok(), true)
  infile = result.ValueOrDie();

  std::unique_ptr<parquet::arrow::FileReader> reader;
  PARQUET_ASSIGN_OR_THROW(reader, parquet::arrow::OpenFile(infile, pool));

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  TEST_EQUAL(read_status.ok(), true)

  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 24)

  // Verify modifications column has structured data for modified peptide
  auto mod_col = table->GetColumnByName("modifications");
  TEST_NOT_EQUAL(mod_col, nullptr)
  TEST_EQUAL(mod_col->type()->id(), arrow::Type::LIST)

  // Verify protein_accessions is a list with two entries
  auto pa_col = table->GetColumnByName("protein_accessions");
  TEST_NOT_EQUAL(pa_col, nullptr)
  TEST_EQUAL(pa_col->type()->id(), arrow::Type::LIST)

  // Verify mass_error_ppm column exists
  auto mep_col = table->GetColumnByName("mass_error_ppm");
  TEST_NOT_EQUAL(mep_col, nullptr)
  TEST_EQUAL(mep_col->type()->id(), arrow::Type::FLOAT)

  // Verify QPX file metadata is present in schema key-value metadata
  auto metadata = table->schema()->metadata();
  TEST_NOT_EQUAL(metadata, nullptr)

  auto qpx_version_idx = metadata->FindKey("qpx_version");
  TEST_EQUAL(qpx_version_idx >= 0, true)
  TEST_EQUAL(metadata->value(qpx_version_idx), "1.0")

  auto creator_idx = metadata->FindKey("creator");
  TEST_EQUAL(creator_idx >= 0, true)
  TEST_EQUAL(metadata->value(creator_idx), "OpenMS")

  auto file_type_idx = metadata->FindKey("file_type");
  TEST_EQUAL(file_type_idx >= 0, true)
  TEST_EQUAL(metadata->value(file_type_idx), "psm_file")

  // QPX serialization spec: compression_format must be one of zstd|snappy|gzip|lzo|none
  auto compression_idx = metadata->FindKey("compression_format");
  TEST_TRUE(compression_idx >= 0)
  TEST_EQUAL(metadata->value(compression_idx), "zstd")

  // scan_format is derived from the spectrum native IDs. This fixture sets none, so the
  // key must be absent rather than asserted as "scan" — the scan column cannot hold scan
  // numbers when there is no native ID to extract them from.
  TEST_TRUE(metadata->FindKey("scan_format") < 0)

  // software_provider carries the versioned library ("OpenMS <version>"); creator is the
  // bare tool/org name.
  auto software_idx = metadata->FindKey("software_provider");
  TEST_TRUE(software_idx >= 0)
  TEST_TRUE(StringUtils::hasPrefix(metadata->value(software_idx), "OpenMS "))

  // creation_date and uuid should exist (values are dynamic)
  TEST_EQUAL(metadata->FindKey("creation_date") >= 0, true)
  TEST_EQUAL(metadata->FindKey("uuid") >= 0, true)
}
END_SECTION

START_SECTION([EXTRA] QPX scan_format is derived from the spectrum native IDs)
{
  // Helper: write PSMs whose spectrum references use `ref_prefix`, read the file back and
  // return its scan_format metadata value ("" when the key is absent).
  auto scan_format_for = [](const std::string& ref) -> std::string
  {
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;

    ProteinIdentification protein_id;
    protein_id.setIdentifier("sf_search");
    protein_ids.push_back(protein_id);

    PeptideIdentification pid;
    pid.setIdentifier("sf_search");
    pid.setScoreType("TestScore");
    if (!ref.empty()) { pid.setSpectrumReference(ref); }
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDER"));
    hit.setCharge(2);
    hit.setScore(1.0);
    pid.setHits({hit});
    peptide_ids.push_back(pid);

    std::string f;
    NEW_TMP_FILE(f)
    if (!QPXFile::exportToParquet(protein_ids, peptide_ids, f)) { return "<write failed>"; }

    auto open_res = arrow::io::ReadableFile::Open(f.c_str());
    if (!open_res.ok()) { return "<read failed>"; }
    auto reader_res = parquet::arrow::OpenFile(open_res.ValueOrDie(), arrow::default_memory_pool());
    if (!reader_res.ok()) { return "<read failed>"; }
    auto reader = std::move(reader_res).ValueOrDie();
    std::shared_ptr<arrow::Table> t;
    if (!reader->ReadTable(&t).ok()) { return "<read failed>"; }

    auto md = t->schema()->metadata();
    if (!md) { return ""; }
    auto r = md->Get("scan_format");
    return r.ok() ? r.ValueOrDie() : "";
  };

  TEST_STRING_EQUAL(scan_format_for("controllerType=0 controllerNumber=1 scan=1234"), "scan")
  TEST_STRING_EQUAL(scan_format_for("scan=42"), "scan")
  TEST_STRING_EQUAL(scan_format_for("index=7"), "index")
  // Unrecognized / absent native IDs yield no evidence, so the key is omitted.
  TEST_STRING_EQUAL(scan_format_for("some_opaque_identifier"), "")
  TEST_STRING_EQUAL(scan_format_for(""), "")
}
END_SECTION

START_SECTION(Test modifications structured output)
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(100.0);
  peptide_id.setMZ(500.0);
  peptide_id.setScoreType("TestScore");

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEM(Oxidation)TIDER"));
  hit.setScore(0.9);
  hit.setCharge(2);

  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);

  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify modifications column is a list of structs (not null)
  auto mod_col = table->GetColumnByName("modifications");
  auto mod_list = std::static_pointer_cast<arrow::ListArray>(mod_col->chunk(0));
  // The modified peptide should have at least one modification entry
  TEST_EQUAL(mod_list->value_length(0) > 0, true)

  // Unmodified peptide should have empty modifications list
  PeptideIdentificationList peptide_ids2;
  PeptideIdentification peptide_id2;
  peptide_id2.setIdentifier("test_search");
  peptide_id2.setRT(100.0);
  peptide_id2.setMZ(500.0);
  peptide_id2.setScoreType("TestScore");
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setScore(0.9);
  hit2.setCharge(2);
  peptide_id2.setHits(vector<PeptideHit>{hit2});
  peptide_ids2.push_back(peptide_id2);

  auto table2 = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids2);
  auto mod_col2 = table2->GetColumnByName("modifications");
  auto mod_list2 = std::static_pointer_cast<arrow::ListArray>(mod_col2->chunk(0));
  TEST_EQUAL(mod_list2->value_length(0), 0) // empty list for unmodified
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(QPXPgSchema::schema())
{
  auto schema = QPXPgSchema::schema();
  TEST_EQUAL(schema->num_fields(), 20)

  // Required (non-nullable) fields
  TEST_EQUAL(schema->GetFieldByName("pg_accessions")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("anchor_protein")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("run_file_name")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("is_decoy")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("peptides")->nullable(), false)

  // Nullable (optional) fields
  TEST_EQUAL(schema->GetFieldByName("intensities")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("additional_intensities")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("global_qvalue")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("pg_qvalue")->nullable(), true)
}
END_SECTION

START_SECTION(ProteinGroupArrowExport::exportToArrow(vector<ProteinIdentification>, PeptideIdentificationList))
{
  // Set up minimal protein identification with one group
  ProteinIdentification prot_id;
  prot_id.setPrimaryMSRunPath({"test_run.mzML"});
  ProteinHit hit1;
  hit1.setAccession("PROT_A");
  hit1.setScore(0.99);
  hit1.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  ProteinHit hit2;
  hit2.setAccession("PROT_B");
  hit2.setScore(0.95);
  hit2.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  prot_id.setHits({hit1, hit2});
  prot_id.setScoreType("Mascot");
  prot_id.setHigherScoreBetter(true);

  ProteinIdentification::ProteinGroup group;
  group.accessions = {"PROT_A", "PROT_B"};
  group.probability = 0.01;
  prot_id.insertIndistinguishableProteins(group);

  // Set up peptide identifications
  PeptideIdentification pep_id;
  PeptideHit pep_hit;
  pep_hit.setSequence(AASequence::fromString("PEPTIDE"));
  PeptideEvidence ev;
  ev.setProteinAccession("PROT_A");
  pep_hit.setPeptideEvidences({ev});
  pep_id.setHits({pep_hit});
  PeptideIdentificationList pep_ids = {pep_id};

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, pep_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 20)

  // Verify run_file_name is derived from ProteinIdentification, without path or extension
  auto run_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(run_col->GetString(0), "test_run")

  // Verify anchor_protein
  auto anchor_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("anchor_protein")->chunk(0));
  TEST_STRING_EQUAL(anchor_col->GetString(0), "PROT_A")

  // Verify intensities is null (no quantification)
  auto int_col = table->GetColumnByName("intensities")->chunk(0);
  TEST_EQUAL(int_col->IsNull(0), true)

  // is_decoy
  auto is_decoy_col = std::static_pointer_cast<arrow::BooleanArray>(
    table->GetColumnByName("is_decoy")->chunk(0));
  TEST_EQUAL(is_decoy_col->Value(0), false)

  // global_qvalue carries group.probability
  auto qv_col = std::static_pointer_cast<arrow::DoubleArray>(
    table->GetColumnByName("global_qvalue")->chunk(0));
  TEST_EQUAL(qv_col->IsNull(0), false)
  TEST_REAL_SIMILAR(qv_col->Value(0), 0.01)

  // peptides list: one entry per group accession (2 here: PROT_A and PROT_B)
  auto peptides_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("peptides")->chunk(0));
  TEST_EQUAL(peptides_list->value_length(0), 2)
  auto peptides_struct = std::static_pointer_cast<arrow::StructArray>(peptides_list->values());
  auto pep_name_arr = std::static_pointer_cast<arrow::StringArray>(peptides_struct->field(0));
  auto pep_count_arr = std::static_pointer_cast<arrow::Int32Array>(peptides_struct->field(1));
  TEST_STRING_EQUAL(pep_name_arr->GetString(0), "PROT_A")
  TEST_EQUAL(pep_count_arr->Value(0), 1)   // one peptide evidence for PROT_A
  TEST_STRING_EQUAL(pep_name_arr->GetString(1), "PROT_B")
  TEST_EQUAL(pep_count_arr->Value(1), 0)   // no peptide evidence for PROT_B

  // peptide_counts: unique sequences = 1 (the peptide is attributed to PROT_A only);
  // total_sequences = 1 (sum of per-protein counts)
  auto pc_struct = std::static_pointer_cast<arrow::StructArray>(
    table->GetColumnByName("peptide_counts")->chunk(0));
  auto pc_unique = std::static_pointer_cast<arrow::Int32Array>(pc_struct->field(0));
  auto pc_total = std::static_pointer_cast<arrow::Int32Array>(pc_struct->field(1));
  TEST_EQUAL(pc_unique->Value(0), 1)
  TEST_EQUAL(pc_total->Value(0), 1)

  // gg_accessions / gg_names must be parallel to pg_accessions (length 2, both empty strings)
  auto gg_acc_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("gg_accessions")->chunk(0));
  auto gg_names_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("gg_names")->chunk(0));
  TEST_EQUAL(gg_acc_list->value_length(0), 2)
  TEST_EQUAL(gg_names_list->value_length(0), 2)
}
END_SECTION

START_SECTION(ProteinGroupArrowExport::exportToArrow empty groups)
{
  ProteinIdentification prot_id;
  prot_id.setPrimaryMSRunPath({"test.mzML"});
  PeptideIdentificationList pep_ids;

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, pep_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 20)
}
END_SECTION

START_SECTION(([EXTRA] ProteinGroupArrowExport writes QPX channel labels, not channel numbers))
{
  // intensities[].label is a join key (i.label = rs.label against run.samples), so it must
  // carry a canonical channel token. This covers the isobaric path, which has no TOPP test.
  //
  // Deliberately does NOT call setExperimentType: IsobaricWorkflow never does either, so the
  // mapping has to work off the headers' channel_id meta value alone.
  auto build_map = [](const std::string& method, const std::vector<std::string>& channel_names)
  {
    ConsensusMap cmap;
    for (Size c = 0; c < channel_names.size(); ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "/data/tmt_run.mzML";
      h.label = method + "_" + channel_names[c];
      h.setMetaValue("channel_name", channel_names[c]);
      h.setMetaValue("channel_id", static_cast<int>(c));   // 0-based; getLabelAsUInt adds 1
      cmap.getColumnHeaders()[c] = h;
    }

    ProteinIdentification prot_id;
    prot_id.setScoreType("q-value");
    prot_id.setPrimaryMSRunPath({"/data/tmt_run.mzML"});
    ProteinHit ph;
    ph.setAccession("PROT_A");
    prot_id.setHits({ph});

    ProteinIdentification::ProteinGroup group;
    group.accessions = {"PROT_A"};
    group.probability = 0.01;
    // Quant arrays as PeptideAndProteinQuant lays them out: abundances in float[3],
    // design filenames (already stemmed) in string[0], 1-based channels in int[0].
    group.getFloatDataArrays().resize(4);
    group.getStringDataArrays().resize(1);
    group.getIntegerDataArrays().resize(1);
    for (Size c = 0; c < channel_names.size(); ++c)
    {
      group.getFloatDataArrays()[3].push_back(1000.0f * (c + 1));
      group.getStringDataArrays()[0].push_back("tmt_run");
      group.getIntegerDataArrays()[0].push_back(static_cast<int>(c) + 1);
    }
    prot_id.insertIndistinguishableProteins(group);
    cmap.setProteinIdentifications({prot_id});
    return cmap;
  };

  auto labels_of = [](const std::shared_ptr<arrow::Table>& t)
  {
    std::vector<std::string> out;
    auto col = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("intensities")->chunk(0));
    auto st = std::static_pointer_cast<arrow::StructArray>(col->values());
    auto lab = std::static_pointer_cast<arrow::StringArray>(st->field(0));
    for (int64_t i = 0; i < lab->length(); ++i) { out.push_back(lab->GetString(i)); }
    return out;
  };

  // TMT 10-plex: channel 10 is "131" in OpenMS' naming (qpx's 11-plex-indexed map says
  // "TMT131N"); OpenMS' name is authoritative.
  {
    auto cmap = build_map("tmt10plex",
      {"126","127N","127C","128N","128C","129N","129C","130N","130C","131"});
    auto t = ProteinGroupArrowExport::exportToArrow(cmap);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), 1)
    auto labels = labels_of(t);
    TEST_EQUAL(labels.size(), 10)
    TEST_STRING_EQUAL(labels[0], "TMT126")
    TEST_STRING_EQUAL(labels[1], "TMT127N")
    TEST_STRING_EQUAL(labels[9], "TMT131")
    // run_file_name is the stemmed design filename, matching psm/feature
    auto run_col = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
    TEST_STRING_EQUAL(run_col->GetString(0), "tmt_run")
  }

  // iTRAQ 4-plex
  {
    auto cmap = build_map("itraq4plex", {"114","115","116","117"});
    auto t = ProteinGroupArrowExport::exportToArrow(cmap);
    TEST_NOT_EQUAL(t, nullptr)
    auto labels = labels_of(t);
    TEST_EQUAL(labels.size(), 4)
    TEST_STRING_EQUAL(labels[0], "ITRAQ114")
    TEST_STRING_EQUAL(labels[3], "ITRAQ117")
  }

  // An unidentifiable quantitation method must abort rather than emit a guessed join key.
  {
    auto cmap = build_map("not_a_method", {"126", "127N"});
    TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap), nullptr)
  }
}
END_SECTION

START_SECTION(([EXTRA] importFromArrow_round_trip))
{
  // Build minimal protein + peptide identifications.
  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  std::vector<ProteinIdentification> prot_ids{prot};

  // --- First PeptideIdentification ---
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

  // --- Second PeptideIdentification (same run, different spectrum) ---
  PeptideIdentification pid2;
  pid2.setIdentifier("run_1");
  pid2.setScoreType("score");
  pid2.setHigherScoreBetter(true);
  pid2.setRT(234.5);
  pid2.setMZ(678.90);
  pid2.setSpectrumReference("scan=99");
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit2.setCharge(3);
  hit2.setScore(0.80);
  hit2.setMetaValue("target_decoy", "target");
  pid2.getHits().push_back(hit2);

  PeptideIdentificationList pep_ids;
  pep_ids.push_back(pid);
  pep_ids.push_back(pid2);

  // Export to Arrow table (PSMSchema), then import back.
  // export_all_psms=true is required so that every PeptideIdentification produces a row.
  auto table = QPXFile::exportToArrow(prot_ids, pep_ids, /*export_all_psms=*/true);
  TEST_NOT_EQUAL(table.get(), nullptr);
  // Two PIDs → two rows in the table.
  TEST_EQUAL(table->num_rows(), 2);

  std::vector<ProteinIdentification> prot_ids_out = prot_ids; // shell with run identifier preserved
  PeptideIdentificationList pep_ids_out;
  TEST_TRUE(QPXFile::importFromArrow(table, prot_ids_out, pep_ids_out));

  // --- Single-group assertions (existing) ---
  TEST_EQUAL(pep_ids_out.size(), 2);
  TEST_EQUAL(pep_ids_out[0].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_out[0].getHits()[0].getSequence().toString(), "PEPTIDE");
  TEST_EQUAL(pep_ids_out[0].getHits()[0].getCharge(), 2);
  TEST_REAL_SIMILAR(pep_ids_out[0].getHits()[0].getScore(), 0.95);
  TEST_STRING_EQUAL(StringUtils::toStr(pep_ids_out[0].getHits()[0].getMetaValue("target_decoy")), "target");
  TEST_REAL_SIMILAR(double(pep_ids_out[0].getHits()[0].getMetaValue("COMET:deltaCn")), 0.5);

  // --- Round-trip detail assertions on the first PID ---
  TEST_REAL_SIMILAR(pep_ids_out[0].getRT(), 123.4);
  TEST_REAL_SIMILAR(pep_ids_out[0].getMZ(), 567.89);
  TEST_STRING_EQUAL(pep_ids_out[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(pep_ids_out[0].getSpectrumReference(), "scan=42");
  // Protein accession + positional fields round-trip on the first hit.
  TEST_EQUAL(pep_ids_out[0].getHits()[0].getPeptideEvidences().size(), 1);
  {
    const auto& ev_rt = pep_ids_out[0].getHits()[0].getPeptideEvidences()[0];
    TEST_STRING_EQUAL(ev_rt.getProteinAccession(), "sp|P12345|EXAMPLE");
    TEST_EQUAL(ev_rt.getAABefore(), 'K');
    TEST_EQUAL(ev_rt.getAAAfter(),  'R');
    TEST_EQUAL(ev_rt.getStart(), 42);
    TEST_EQUAL(ev_rt.getEnd(),   48);
  }

  // --- Multi-group dispatch assertions (second PID) ---
  TEST_EQUAL(pep_ids_out[1].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_out[1].getHits()[0].getSequence().toString(), "ACDEFGHIK");
  TEST_EQUAL(pep_ids_out[1].getHits()[0].getCharge(), 3);
  TEST_REAL_SIMILAR(pep_ids_out[1].getHits()[0].getScore(), 0.80);
  TEST_REAL_SIMILAR(pep_ids_out[1].getRT(), 234.5);
  TEST_REAL_SIMILAR(pep_ids_out[1].getMZ(), 678.90);
  TEST_STRING_EQUAL(pep_ids_out[1].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(pep_ids_out[1].getSpectrumReference(), "scan=99");
}
END_SECTION

START_SECTION(([EXTRA] importFromArrow_round_trip_peptide_evidence_unknown_sentinels))
{
  // PSMSchema protein_accessions is list<struct{accession, aa_before, aa_after, start, end}>
  // with UNKNOWN_AA / UNKNOWN_POSITION encoded as Arrow null. Verify round-trip restores
  // the C++ sentinels.
  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  std::vector<ProteinIdentification> prot_ids{prot};

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  pid.setRT(11.0);
  pid.setMZ(222.0);
  pid.setSpectrumReference("scan=1");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(0.5);
  hit.setMetaValue("target_decoy", "target");

  // Default-constructed PeptideEvidence: UNKNOWN_AA + UNKNOWN_POSITION sentinels.
  PeptideEvidence ev;
  ev.setProteinAccession("PROT_X");
  TEST_EQUAL(ev.getAABefore(), PeptideEvidence::UNKNOWN_AA);
  TEST_EQUAL(ev.getAAAfter(),  PeptideEvidence::UNKNOWN_AA);
  TEST_EQUAL(ev.getStart(), PeptideEvidence::UNKNOWN_POSITION);
  TEST_EQUAL(ev.getEnd(),   PeptideEvidence::UNKNOWN_POSITION);
  hit.addPeptideEvidence(ev);

  // Also test N/C terminal markers — must round-trip as literal '['/']'.
  PeptideEvidence ev_term;
  ev_term.setProteinAccession("PROT_TERM");
  ev_term.setAABefore(PeptideEvidence::N_TERMINAL_AA);
  ev_term.setAAAfter(PeptideEvidence::C_TERMINAL_AA);
  ev_term.setStart(0);
  ev_term.setEnd(6);
  hit.addPeptideEvidence(ev_term);

  pid.getHits().push_back(hit);
  PeptideIdentificationList pep_ids{pid};

  auto table = QPXFile::exportToArrow(prot_ids, pep_ids, /*export_all_psms=*/true);
  TEST_NOT_EQUAL(table.get(), nullptr);

  std::vector<ProteinIdentification> prot_ids_out = prot_ids;
  PeptideIdentificationList pep_ids_out;
  TEST_TRUE(QPXFile::importFromArrow(table, prot_ids_out, pep_ids_out));

  TEST_EQUAL(pep_ids_out.size(), 1);
  TEST_EQUAL(pep_ids_out[0].getHits().size(), 1);
  const auto& evs = pep_ids_out[0].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(evs.size(), 2);

  // First evidence: unknown sentinels survive null round-trip.
  TEST_STRING_EQUAL(evs[0].getProteinAccession(), "PROT_X");
  TEST_EQUAL(evs[0].getAABefore(), PeptideEvidence::UNKNOWN_AA);
  TEST_EQUAL(evs[0].getAAAfter(),  PeptideEvidence::UNKNOWN_AA);
  TEST_EQUAL(evs[0].getStart(), PeptideEvidence::UNKNOWN_POSITION);
  TEST_EQUAL(evs[0].getEnd(),   PeptideEvidence::UNKNOWN_POSITION);

  // Second evidence: literal terminal markers round-trip as length-1 strings, not as null.
  TEST_STRING_EQUAL(evs[1].getProteinAccession(), "PROT_TERM");
  TEST_EQUAL(evs[1].getAABefore(), PeptideEvidence::N_TERMINAL_AA);
  TEST_EQUAL(evs[1].getAAAfter(),  PeptideEvidence::C_TERMINAL_AA);
  TEST_EQUAL(evs[1].getStart(), 0);
  TEST_EQUAL(evs[1].getEnd(),   6);
}
END_SECTION

START_SECTION(([EXTRA] importFromArrow_appends_shell_for_unknown_run_identifier))
{
  // Build one PeptideIdentification with a run_identifier that has NO matching
  // ProteinIdentification — importFromArrow must append a shell entry.
  PeptideIdentification pid;
  pid.setIdentifier("run_unknown");
  pid.setScoreType("hyperscore");
  pid.setHigherScoreBetter(true);
  pid.setRT(99.0);
  pid.setMZ(400.0);
  pid.setSpectrumReference("scan=7");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("TESTVAK"));
  hit.setCharge(2);
  hit.setScore(0.70);
  hit.setMetaValue("target_decoy", "target");
  pid.getHits().push_back(hit);

  PeptideIdentificationList pep_ids;
  pep_ids.push_back(pid);

  // Export: no matching ProteinIdentification provided.
  std::vector<ProteinIdentification> prot_ids_empty;
  auto table = QPXFile::exportToArrow(prot_ids_empty, pep_ids, /*export_all_psms=*/true);
  TEST_NOT_EQUAL(table.get(), nullptr);

  // Import against an *empty* prot_ids_out vector — shell must be appended.
  std::vector<ProteinIdentification> prot_ids_out;
  PeptideIdentificationList pep_ids_out;
  TEST_TRUE(QPXFile::importFromArrow(table, prot_ids_out, pep_ids_out));

  // One shell ProteinIdentification must have been appended.
  TEST_EQUAL(prot_ids_out.size(), 1);
  TEST_STRING_EQUAL(prot_ids_out[0].getIdentifier(), "run_unknown");

  // The shell's score_type and higher_score_better must mirror the input PID.
  TEST_STRING_EQUAL(prot_ids_out[0].getScoreType(), "hyperscore");
  TEST_EQUAL(prot_ids_out[0].isHigherScoreBetter(), true);

  // The PeptideIdentification itself must have been populated correctly.
  TEST_EQUAL(pep_ids_out.size(), 1);
  TEST_STRING_EQUAL(pep_ids_out[0].getIdentifier(), "run_unknown");
  TEST_EQUAL(pep_ids_out[0].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_out[0].getHits()[0].getSequence().toString(), "TESTVAK");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static bool exportToParquetStreaming(const std::vector<ProteinIdentification>&, const std::vector<const PeptideIdentification*>&, const std::string&, bool, size_t, const ParquetWriteConfig&)))
{
  arrow::MemoryPool* pool = arrow::default_memory_pool();

  // Read a parquet file into a single-chunk (combined) Arrow table. Streaming writes
  // multiple row groups, so columns come back multi-chunk; combine before row indexing.
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

  // --- Build test data (mix of modified/unmodified, decoys, evidences) ---
  vector<ProteinIdentification> protein_ids;
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_id.setPrimaryMSRunPath({"/data/run1.mzML"});
  protein_ids.push_back(protein_id);

  PeptideIdentificationList peptide_ids;
  const size_t M = 2500;
  for (size_t i = 0; i < M; ++i)
  {
    PeptideIdentification pid;
    pid.setIdentifier("test_search");
    pid.setRT(100.0 + i);
    pid.setMZ(400.0 + (i % 500));
    pid.setScoreType("TestScore");
    pid.setHigherScoreBetter(true);
    pid.setSpectrumReference("controllerType=0 controllerNumber=1 scan=" + StringUtils::toStr(1000 + i));

    PeptideHit hit;
    hit.setSequence(AASequence::fromString((i % 3 == 0) ? "PEM(Oxidation)TIDER" : "PEPTIDEK"));
    hit.setCharge(2 + (i % 3));
    hit.setScore(0.99 - (i % 100) * 0.001);
    hit.setMetaValue("target_decoy", (i % 2 == 0) ? "target" : "decoy");
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_" + StringUtils::toStr(i % 50));
    hit.setPeptideEvidences(vector<PeptideEvidence>{ev});
    pid.setHits(vector<PeptideHit>{hit});
    peptide_ids.push_back(pid);
  }

  std::vector<const PeptideIdentification*> ptrs;
  ptrs.reserve(peptide_ids.size());
  for (const auto& p : peptide_ids) { ptrs.push_back(&p); }

  // --- Streaming write with a small batch (multiple row groups: 2500 / 1000) ---
  std::string stream_file;
  NEW_TMP_FILE(stream_file)
  TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, ptrs, stream_file, false, 1000), true)

  auto st_table = read_combined(stream_file);
  TEST_NOT_EQUAL(st_table, nullptr)
  TEST_EQUAL(st_table->num_rows(), (int64_t)M)
  TEST_EQUAL(st_table->num_columns(), 24)

  // --- QPX metadata must survive the streaming/metadata-Open path ---
  {
    auto md = st_table->schema()->metadata();
    TEST_NOT_EQUAL(md, nullptr)
    std::string ft, ver;
    if (md)
    {
      auto r1 = md->Get("file_type");   if (r1.ok()) { ft = r1.ValueOrDie(); }
      auto r2 = md->Get("qpx_version"); if (r2.ok()) { ver = r2.ValueOrDie(); }
    }
    TEST_STRING_EQUAL(ft, "psm_file")
    TEST_STRING_EQUAL(ver, "1.0")
  }

  // --- Equivalence vs the one-shot (non-streaming) path ---
  std::string ref_file;
  NEW_TMP_FILE(ref_file)
  TEST_EQUAL(QPXFile::exportToParquet(protein_ids, peptide_ids, ref_file), true)
  auto ref_table = read_combined(ref_file);
  TEST_NOT_EQUAL(ref_table, nullptr)
  TEST_EQUAL(ref_table->num_rows(), st_table->num_rows())

  {
    auto seq_s = std::static_pointer_cast<arrow::StringArray>(st_table->GetColumnByName("sequence")->chunk(0));
    auto seq_r = std::static_pointer_cast<arrow::StringArray>(ref_table->GetColumnByName("sequence")->chunk(0));
    auto pf_s  = std::static_pointer_cast<arrow::StringArray>(st_table->GetColumnByName("peptidoform")->chunk(0));
    auto pf_r  = std::static_pointer_cast<arrow::StringArray>(ref_table->GetColumnByName("peptidoform")->chunk(0));
    auto chg_s = std::static_pointer_cast<arrow::Int16Array>(st_table->GetColumnByName("charge")->chunk(0));
    auto chg_r = std::static_pointer_cast<arrow::Int16Array>(ref_table->GetColumnByName("charge")->chunk(0));
    auto sc_s  = std::static_pointer_cast<arrow::DoubleArray>(st_table->GetColumnByName("posterior_error_probability")->chunk(0));
    auto sc_r  = std::static_pointer_cast<arrow::DoubleArray>(ref_table->GetColumnByName("posterior_error_probability")->chunk(0));
    // rt is unique per row (100.0 + i) -> comparing it makes this an ORDER-sensitive check that
    // catches partition-reordering regressions, not just value equality.
    auto rt_s  = std::static_pointer_cast<arrow::FloatArray>(st_table->GetColumnByName("rt")->chunk(0));
    auto rt_r  = std::static_pointer_cast<arrow::FloatArray>(ref_table->GetColumnByName("rt")->chunk(0));
    bool all_eq = true;
    for (int64_t r = 0; r < st_table->num_rows(); ++r)
    {
      if (seq_s->GetString(r) != seq_r->GetString(r)) { all_eq = false; break; }
      if (pf_s->GetString(r)  != pf_r->GetString(r))  { all_eq = false; break; }
      if (chg_s->Value(r)     != chg_r->Value(r))     { all_eq = false; break; }
      if (rt_s->Value(r)      != rt_r->Value(r))      { all_eq = false; break; } // order-sensitive
      if (sc_s->IsNull(r) != sc_r->IsNull(r))         { all_eq = false; break; }
      if (!sc_s->IsNull(r) && sc_s->Value(r) != sc_r->Value(r)) { all_eq = false; break; }
    }
    TEST_EQUAL(all_eq, true)
  }

  // --- Edge case: empty input -> valid 0-row file with full schema + metadata ---
  {
    std::vector<const PeptideIdentification*> empty_ptrs;
    std::string empty_file;
    NEW_TMP_FILE(empty_file)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, empty_ptrs, empty_file), true)
    auto e_table = read_combined(empty_file);
    TEST_NOT_EQUAL(e_table, nullptr)
    TEST_EQUAL(e_table->num_rows(), 0)
    TEST_EQUAL(e_table->num_columns(), 24)
  }

  // --- Edge case: M=1 with batch_size=1 ---
  {
    std::vector<const PeptideIdentification*> one_ptr{ptrs[0]};
    std::string one_file;
    NEW_TMP_FILE(one_file)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, one_ptr, one_file, false, 1), true)
    auto o_table = read_combined(one_file);
    TEST_NOT_EQUAL(o_table, nullptr)
    TEST_EQUAL(o_table->num_rows(), 1)
  }

  // --- Edge case: batch_size=0 must not hang (guarded to default) and write all rows ---
  {
    std::string zero_file;
    NEW_TMP_FILE(zero_file)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, ptrs, zero_file, false, 0), true)
    auto z_table = read_combined(zero_file);
    TEST_NOT_EQUAL(z_table, nullptr)
    TEST_EQUAL(z_table->num_rows(), (int64_t)M)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] exportToParquetStreaming parallel build (n_threads)))
{
  arrow::MemoryPool* pool = arrow::default_memory_pool();

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

  // Compare two combined tables row-by-row on stable columns. Returns true if equal.
  auto tables_equal = [](const std::shared_ptr<arrow::Table>& a, const std::shared_ptr<arrow::Table>& b) -> bool
  {
    if (!a || !b) return false;
    if (a->num_rows() != b->num_rows()) return false;
    if (a->num_rows() == 0) return true;
    auto sa = std::static_pointer_cast<arrow::StringArray>(a->GetColumnByName("sequence")->chunk(0));
    auto sb = std::static_pointer_cast<arrow::StringArray>(b->GetColumnByName("sequence")->chunk(0));
    auto pa = std::static_pointer_cast<arrow::StringArray>(a->GetColumnByName("peptidoform")->chunk(0));
    auto pb = std::static_pointer_cast<arrow::StringArray>(b->GetColumnByName("peptidoform")->chunk(0));
    auto ca = std::static_pointer_cast<arrow::Int16Array>(a->GetColumnByName("charge")->chunk(0));
    auto cb = std::static_pointer_cast<arrow::Int16Array>(b->GetColumnByName("charge")->chunk(0));
    // rt is unique per row (100.0 + i) -> order-sensitive check (catches partition reordering).
    auto ra = std::static_pointer_cast<arrow::FloatArray>(a->GetColumnByName("rt")->chunk(0));
    auto rb = std::static_pointer_cast<arrow::FloatArray>(b->GetColumnByName("rt")->chunk(0));
    for (int64_t r = 0; r < a->num_rows(); ++r)
    {
      if (sa->GetString(r) != sb->GetString(r)) return false;
      if (pa->GetString(r) != pb->GetString(r)) return false;
      if (ca->Value(r) != cb->Value(r)) return false;
      if (ra->Value(r) != rb->Value(r)) return false;
    }
    return true;
  };

  vector<ProteinIdentification> protein_ids;
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_id.setPrimaryMSRunPath({"/data/run1.mzML"});
  protein_ids.push_back(protein_id);

  // --- Single-hit dataset (best-hit export) ---
  PeptideIdentificationList peptide_ids;
  const size_t M = 2500;
  for (size_t i = 0; i < M; ++i)
  {
    PeptideIdentification pid;
    pid.setIdentifier("test_search");
    pid.setRT(100.0 + i);
    pid.setMZ(400.0 + (i % 500));
    pid.setScoreType("TestScore");
    pid.setHigherScoreBetter(true);
    pid.setSpectrumReference("controllerType=0 controllerNumber=1 scan=" + StringUtils::toStr(1000 + i));
    PeptideHit hit;
    hit.setSequence(AASequence::fromString((i % 3 == 0) ? "PEM(Oxidation)TIDEK" : "PEPTIDEK"));
    hit.setCharge(2 + (i % 3));
    hit.setScore(0.99 - (i % 100) * 0.001);
    hit.setMetaValue("target_decoy", (i % 2 == 0) ? "target" : "decoy");
    PeptideEvidence ev;
    ev.setProteinAccession("PROT_" + StringUtils::toStr(i % 50));
    hit.setPeptideEvidences(vector<PeptideEvidence>{ev});
    pid.setHits(vector<PeptideHit>{hit});
    peptide_ids.push_back(pid);
  }
  std::vector<const PeptideIdentification*> ptrs;
  ptrs.reserve(peptide_ids.size());
  for (const auto& p : peptide_ids) { ptrs.push_back(&p); }

  // Reference: one-shot (non-streaming) output.
  std::string ref_file; NEW_TMP_FILE(ref_file)
  TEST_EQUAL(QPXFile::exportToParquet(protein_ids, peptide_ids, ref_file), true)
  auto ref_table = read_combined(ref_file);
  TEST_NOT_EQUAL(ref_table, nullptr)

  // Determinism/equivalence across thread counts (batch_size small => many batches x partitions).
  for (int nthreads : {1, 2, 8})
  {
    std::string f; NEW_TMP_FILE(f)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, ptrs, f, false, 1000, ParquetWriteConfig{}, nthreads), true)
    auto tbl = read_combined(f);
    TEST_NOT_EQUAL(tbl, nullptr)
    TEST_EQUAL(tbl->num_rows(), (int64_t)M)
    TEST_EQUAL(tables_equal(tbl, ref_table), true)
  }

  // --- export_all_psms=true with variable hit counts across partition boundaries ---
  PeptideIdentificationList multi_ids;
  size_t expected_rows = 0;
  for (size_t i = 0; i < 777; ++i)
  {
    PeptideIdentification pid;
    pid.setIdentifier("test_search");
    pid.setRT(10.0 + i);
    pid.setMZ(300.0 + (i % 100));
    pid.setScoreType("TestScore");
    pid.setHigherScoreBetter(true);
    pid.setSpectrumReference("controllerType=0 controllerNumber=1 scan=" + StringUtils::toStr(i));
    size_t n_hits = 1 + (i % 3); // 1..3 hits
    vector<PeptideHit> hits;
    for (size_t h = 0; h < n_hits; ++h)
    {
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(h == 0 ? "PEPTIDEK" : "DFPIANGER"));
      hit.setCharge(2 + (int)h);
      hit.setScore(0.9 - 0.01 * h);
      hits.push_back(hit);
    }
    expected_rows += n_hits;
    pid.setHits(hits);
    multi_ids.push_back(pid);
  }
  std::vector<const PeptideIdentification*> multi_ptrs;
  for (const auto& p : multi_ids) { multi_ptrs.push_back(&p); }

  std::string all_ref; NEW_TMP_FILE(all_ref)
  TEST_EQUAL(QPXFile::exportToParquet(protein_ids, multi_ids, all_ref, /*export_all_psms=*/true), true)
  auto all_ref_tbl = read_combined(all_ref);
  TEST_NOT_EQUAL(all_ref_tbl, nullptr)
  TEST_EQUAL(all_ref_tbl->num_rows(), (int64_t)expected_rows)

  std::string all_par; NEW_TMP_FILE(all_par)
  TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, multi_ptrs, all_par, /*export_all_psms=*/true, 50, ParquetWriteConfig{}, 8), true)
  auto all_par_tbl = read_combined(all_par);
  TEST_NOT_EQUAL(all_par_tbl, nullptr)
  TEST_EQUAL(all_par_tbl->num_rows(), (int64_t)expected_rows)
  TEST_EQUAL(tables_equal(all_par_tbl, all_ref_tbl), true)

  // --- Non-empty input that yields zero PSM rows (all hits empty) + 8 threads ---
  {
    PeptideIdentificationList empty_hits(50); // 50 default-constructed PeptideIdentifications, no hits
    std::vector<const PeptideIdentification*> eh_ptrs;
    for (const auto& p : empty_hits) { eh_ptrs.push_back(&p); }
    std::string f; NEW_TMP_FILE(f)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, eh_ptrs, f, false, 10, ParquetWriteConfig{}, 8), true)
    auto tbl = read_combined(f);
    TEST_NOT_EQUAL(tbl, nullptr)
    TEST_EQUAL(tbl->num_rows(), 0)
    TEST_EQUAL(tbl->num_columns(), 24)
  }

  // --- Edge cases with parallelism: M=0, M=1, rows < threads ---
  {
    std::vector<const PeptideIdentification*> empty_ptrs;
    std::string f0; NEW_TMP_FILE(f0)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, empty_ptrs, f0, false, 1000, ParquetWriteConfig{}, 8), true)
    auto t0 = read_combined(f0);
    TEST_NOT_EQUAL(t0, nullptr)
    TEST_EQUAL(t0->num_rows(), 0)

    std::vector<const PeptideIdentification*> one_ptr{ptrs[0]};
    std::string f1; NEW_TMP_FILE(f1)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, one_ptr, f1, false, 1000, ParquetWriteConfig{}, 8), true)
    auto t1 = read_combined(f1);
    TEST_NOT_EQUAL(t1, nullptr)
    TEST_EQUAL(t1->num_rows(), 1)

    // rows (3) < threads (8): W clamps to rows
    std::vector<const PeptideIdentification*> few(ptrs.begin(), ptrs.begin() + 3);
    std::string f3; NEW_TMP_FILE(f3)
    TEST_EQUAL(QPXFile::exportToParquetStreaming(protein_ids, few, f3, false, 1000, ParquetWriteConfig{}, 8), true)
    auto t3 = read_combined(f3);
    TEST_NOT_EQUAL(t3, nullptr)
    TEST_EQUAL(t3->num_rows(), 3)
  }
}
END_SECTION

START_SECTION(([EXTRA] exportToParquetStreaming dedicated-column values (index path + precedence)))
{
  // Locks the Approach-B index-based dedicated-column extraction: exact values + alias precedence,
  // identical between the parallel streaming path and the one-shot path.
  arrow::MemoryPool* pool = arrow::default_memory_pool();
  auto read1 = [pool](const std::string& path) -> std::shared_ptr<arrow::Table>
  {
    auto o = arrow::io::ReadableFile::Open(path.c_str());
    if (!o.ok()) { return nullptr; }
    auto rr = parquet::arrow::OpenFile(o.ValueOrDie(), pool);
    if (!rr.ok()) { return nullptr; }
    auto rd = std::move(rr).ValueOrDie();
    std::shared_ptr<arrow::Table> t;
    if (!rd->ReadTable(&t).ok()) { return nullptr; }
    auto c = t->CombineChunks(pool);
    return c.ok() ? c.ValueOrDie() : nullptr;
  };

  vector<ProteinIdentification> prot(1);
  prot[0].setIdentifier("run1");
  prot[0].setScoreType("q-value");
  prot[0].setHigherScoreBetter(false);
  prot[0].setPrimaryMSRunPath({"/data/fallback.mzML"});

  PeptideIdentificationList peps;
  PeptideIdentification pid;
  pid.setIdentifier("run1");
  pid.setRT(42.0);
  pid.setMZ(555.5);
  pid.setScoreType("q-value"); // NOT a PEP type -> PEP is looked up in hit metavalues
  pid.setHigherScoreBetter(false);
  pid.setSpectrumReference("controllerType=0 controllerNumber=1 scan=4242");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(3);
  hit.setScore(0.001);
  hit.setMetaValue("target_decoy", "decoy");
  hit.setMetaValue("Posterior Error Probability", 0.25); // -> posterior_error_probability column
  hit.setMetaValue("predicted_RT", 40.0);
  hit.setMetaValue("predicted_rt", 999.0);               // precedence: predicted_RT must win
  hit.setMetaValue("ion_mobility", 1.23);
  hit.setMetaValue("IM", 9.99);                          // precedence: ion_mobility must win
  hit.setMetaValue("missed_cleavages", 2);
  hit.setMetaValue("reference_file_name", "/data/real.mzML"); // must override primary MS run path
  pid.setHits(vector<PeptideHit>{hit});

  // Replicate so the parallel path actually partitions (n_threads=4, batch_size=3 => multiple
  // batches AND multiple partitions per batch), then check the first and last row.
  const size_t NREP = 10;
  for (size_t i = 0; i < NREP; ++i) { peps.push_back(pid); }
  std::vector<const PeptideIdentification*> ptrs;
  ptrs.reserve(peps.size());
  for (const auto& p : peps) { ptrs.push_back(&p); }

  std::string sf; NEW_TMP_FILE(sf)
  TEST_EQUAL(QPXFile::exportToParquetStreaming(prot, ptrs, sf, false, 3, ParquetWriteConfig{}, 4), true)
  std::string of; NEW_TMP_FILE(of)
  TEST_EQUAL(QPXFile::exportToParquet(prot, peps, of), true)

  for (const std::string& f : {sf, of})
  {
    auto t = read1(f);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), (int64_t)NREP)
    auto is_decoy = std::static_pointer_cast<arrow::BooleanArray>(t->GetColumnByName("is_decoy")->chunk(0));
    auto prt      = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("predicted_rt")->chunk(0));
    auto im       = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("ion_mobility")->chunk(0));
    auto mc       = std::static_pointer_cast<arrow::Int16Array>(t->GetColumnByName("missed_cleavages")->chunk(0));
    auto rfn      = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
    auto pep      = std::static_pointer_cast<arrow::DoubleArray>(t->GetColumnByName("posterior_error_probability")->chunk(0));
    for (int64_t r : {(int64_t)0, (int64_t)(NREP - 1)})
    {
      TEST_EQUAL(is_decoy->Value(r), true)
      TEST_REAL_SIMILAR(prt->Value(r), 40.0)  // predicted_RT wins over predicted_rt
      TEST_REAL_SIMILAR(im->Value(r), 1.23)   // ion_mobility wins over IM
      TEST_EQUAL(mc->Value(r), 2)
      TEST_EQUAL(rfn->GetString(r), "real") // reference_file_name overrides fallback, and is stemmed
      TEST_REAL_SIMILAR(pep->Value(r), 0.25)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] run_file_name resolves per PSM on a merged run))
{
  // ProteomicsLFQ merges all inputs into one ProteinIdentification whose 'spectra_data' is the
  // ordered concatenation of every input file; each PeptideIdentification points at one of them via
  // 'id_merge_index'. Before this was honored, every PSM was labelled with the run's first file.
  auto prot = mergedRun({"/data/runA.mzML", "/data/runB.mzML.gz"});
  PeptideIdentificationList peps;
  peps.push_back(mergedPSM("PEPTIDEK", 0));
  peps.push_back(mergedPSM("TESTPEPTIDER", 1));
  peps.push_back(mergedPSM("ANOTHERPEPK", 1));
  auto ptrs = ptrsOf(peps);

  std::string one_shot; NEW_TMP_FILE(one_shot)
  TEST_TRUE(QPXFile::exportToParquet(prot, peps, one_shot))
  std::string streamed; NEW_TMP_FILE(streamed)
  TEST_TRUE(QPXFile::exportToParquetStreaming(prot, ptrs, streamed, false, 2, ParquetWriteConfig{}, 4))

  for (const std::string& f : {one_shot, streamed})
  {
    auto t = readQPXParquet(f);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), 3)
    auto rfn = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
    // Bare stems, and the compound '.mzML.gz' extension is stripped as a whole
    TEST_STRING_EQUAL(rfn->GetString(0), "runA")
    TEST_STRING_EQUAL(rfn->GetString(1), "runB")
    TEST_STRING_EQUAL(rfn->GetString(2), "runB")
  }
}
END_SECTION

START_SECTION(([EXTRA] merged run without a usable id_merge_index is refused))
{
  const StringList paths = {"/data/runA.mzML", "/data/runB.mzML"};

  // (a) metavalue missing entirely
  {
    auto prot = mergedRun(paths);
    PeptideIdentificationList peps;
    peps.push_back(mergedPSM("PEPTIDEK", 0));
    peps.push_back(mergedPSM("TESTPEPTIDER", 0, /*set_index=*/false));
    auto ptrs = ptrsOf(peps);

    std::string f1; NEW_TMP_FILE(f1)
    TEST_EXCEPTION(Exception::MissingInformation, QPXFile::exportToParquet(prot, peps, f1))
    std::string f2; NEW_TMP_FILE(f2)
    TEST_EXCEPTION(Exception::MissingInformation,
                   QPXFile::exportToParquetStreaming(prot, ptrs, f2, false, 2, ParquetWriteConfig{}, 4))
    // The check is a serial pre-pass, so the streaming path must not have created (and truncated)
    // the output file before failing.
    TEST_FALSE(File::exists(f2))
  }

  // (b) index past the end of 'spectra_data'
  {
    auto prot = mergedRun(paths);
    PeptideIdentificationList peps;
    peps.push_back(mergedPSM("PEPTIDEK", 2));
    auto ptrs = ptrsOf(peps);
    std::string f1; NEW_TMP_FILE(f1)
    TEST_EXCEPTION(Exception::MissingInformation, QPXFile::exportToParquet(prot, peps, f1))
    std::string f2; NEW_TMP_FILE(f2)
    TEST_EXCEPTION(Exception::MissingInformation,
                   QPXFile::exportToParquetStreaming(prot, ptrs, f2, false, 2, ParquetWriteConfig{}, 4))
  }

  // (c) negative index
  {
    auto prot = mergedRun(paths);
    PeptideIdentificationList peps;
    peps.push_back(mergedPSM("PEPTIDEK", -1));
    std::string f1; NEW_TMP_FILE(f1)
    TEST_EXCEPTION(Exception::MissingInformation, QPXFile::exportToParquet(prot, peps, f1))
  }

  // (d) non-integer index
  {
    auto prot = mergedRun(paths);
    PeptideIdentificationList peps;
    PeptideIdentification pid = mergedPSM("PEPTIDEK", 0);
    pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, "not_an_index");
    peps.push_back(pid);
    std::string f1; NEW_TMP_FILE(f1)
    TEST_EXCEPTION(Exception::MissingInformation, QPXFile::exportToParquet(prot, peps, f1))
  }

  // (e) unmerged input (a single file in the run) needs no index and is unaffected
  {
    auto prot = mergedRun({"/data/only.mzML"});
    PeptideIdentificationList peps;
    peps.push_back(mergedPSM("PEPTIDEK", 0, /*set_index=*/false));
    std::string f1; NEW_TMP_FILE(f1)
    TEST_TRUE(QPXFile::exportToParquet(prot, peps, f1))
    auto t = readQPXParquet(f1);
    TEST_NOT_EQUAL(t, nullptr)
    auto rfn = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
    TEST_STRING_EQUAL(rfn->GetString(0), "only") // only change vs. before: path and extension gone
  }
}
END_SECTION

START_SECTION(([EXTRA] run_file_name stem collisions warn but do not fail the export))
{
  // Same file name in two directories is a legitimate layout. QPX cannot represent the difference,
  // but the origin is known - so warn (not visible here) and export anyway.
  auto prot = mergedRun({"/a/run1.mzML", "/b/run1.mzML"});
  PeptideIdentificationList peps;
  peps.push_back(mergedPSM("PEPTIDEK", 0));
  peps.push_back(mergedPSM("TESTPEPTIDER", 1));

  std::string f; NEW_TMP_FILE(f)
  TEST_TRUE(QPXFile::exportToParquet(prot, peps, f))
  auto t = readQPXParquet(f);
  TEST_NOT_EQUAL(t, nullptr)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "run1")
  TEST_STRING_EQUAL(rfn->GetString(1), "run1")
}
END_SECTION

START_SECTION(([EXTRA] per-PSM origin survives an internal-format round trip))
{
  // The internal PSMSchema stores the origin as a full path in 'reference_file_name', which
  // importFromArrow restores as a per-hit metavalue - and the QPX writer gives that metavalue
  // precedence over the run-derived fallback. If the internal export collapsed to the run's first
  // file, the QPX table would be wrong again after any round trip, even with the merged-run fix.
  auto prot = mergedRun({"/data/runA.mzML", "/data/runB.mzML"});
  PeptideIdentificationList peps;
  peps.push_back(mergedPSM("PEPTIDEK", 0));
  peps.push_back(mergedPSM("TESTPEPTIDER", 1));

  auto internal = QPXFile::exportToArrow(prot, peps, false);
  TEST_NOT_EQUAL(internal, nullptr)
  auto ref = std::static_pointer_cast<arrow::StringArray>(internal->GetColumnByName("reference_file_name")->chunk(0));
  TEST_STRING_EQUAL(ref->GetString(0), "/data/runA.mzML") // full path: the internal format round-trips it
  TEST_STRING_EQUAL(ref->GetString(1), "/data/runB.mzML")

  vector<ProteinIdentification> prot_in;
  PeptideIdentificationList peps_in;
  TEST_TRUE(QPXFile::importFromArrow(internal, prot_in, peps_in))

  auto qpx = QPXFile::exportPSMsToQPXArrow(prot_in, peps_in, false);
  TEST_NOT_EQUAL(qpx, nullptr)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(qpx->GetColumnByName("run_file_name")->chunk(0));
  TEST_EQUAL(qpx->num_rows(), 2)
  TEST_STRING_EQUAL(rfn->GetString(0), "runA")
  TEST_STRING_EQUAL(rfn->GetString(1), "runB")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
