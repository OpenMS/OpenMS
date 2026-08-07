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
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/QPXValueValidation.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

#include <cmath>
#include <set>
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

std::shared_ptr<arrow::Table> replaceColumn(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& name,
  const std::shared_ptr<arrow::Array>& values)
{
  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns = table->columns();
  columns.at(static_cast<size_t>(table->schema()->GetFieldIndex(name))) =
    std::make_shared<arrow::ChunkedArray>(values);
  return arrow::Table::Make(table->schema(), std::move(columns));
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
  TEST_TRUE(table->num_columns() > 0)
}
END_SECTION

START_SECTION(exportToArrow - basic export)
{
  ConsensusMap cmap = createTestConsensusMap();

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);

  TEST_NOT_EQUAL(table, nullptr)
  // The QPX feature view is long: one row per (ConsensusFeature, run). The fixture's cf1
  // spans both runs (handles on maps 0 and 1) while cf2 and cf3 have one handle each,
  // so 3 consensus features melt into 2 + 1 + 1 = 4 rows.
  TEST_EQUAL(table->num_rows(), 4)
  {
    auto run_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
    std::multiset<std::string> runs;
    for (int64_t r = 0; r < run_col->length(); ++r) { runs.insert(run_col->GetString(r)); }
    // stemmed, and cf1 contributes one row to each run
    TEST_EQUAL(runs.count("sample1"), 3)
    TEST_EQUAL(runs.count("sample2"), 1)
  }

  // Check schema has expected columns
  auto schema = table->schema();
  TEST_TRUE(schema->GetFieldIndex("sequence") >= 0)
  TEST_TRUE(schema->GetFieldIndex("peptidoform") >= 0)
  TEST_TRUE(schema->GetFieldIndex("observed_mz") >= 0)
  TEST_TRUE(schema->GetFieldIndex("rt") >= 0)
  TEST_TRUE(schema->GetFieldIndex("charge") >= 0)
  TEST_TRUE(schema->GetFieldIndex("intensities") >= 0)
  TEST_TRUE(schema->GetFieldIndex("pg_accessions") >= 0)
}
END_SECTION

START_SECTION(([EXTRA] feature value validation rejects duplicate primary keys before writing))
{
  ConsensusMap cmap = createTestConsensusMap();
  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)

  QPXValueValidation validator(QPXValueValidation::View::FEATURE);
  TEST_TRUE(validator.validate(table).valid)

  auto duplicate = arrow::ConcatenateTables({table->Slice(0, 1), table->Slice(0, 1)}).ValueOrDie();
  QPXValueValidation duplicate_validator(QPXValueValidation::View::FEATURE);
  auto duplicate_result = duplicate_validator.validate(duplicate);
  TEST_FALSE(duplicate_result.valid)
  TEST_TRUE(duplicate_result.toString().find("primary key") != std::string::npos)

  auto invalid_keys = table->Slice(0, 1);
  arrow::StringBuilder null_peptidoform_builder;
  TEST_TRUE(null_peptidoform_builder.AppendNull().ok())
  invalid_keys = replaceColumn(invalid_keys, QPXFeatureSchema::PEPTIDOFORM,
                               null_peptidoform_builder.Finish().ValueOrDie());
  arrow::Int16Builder null_charge_builder;
  TEST_TRUE(null_charge_builder.AppendNull().ok())
  invalid_keys = replaceColumn(invalid_keys, QPXFeatureSchema::CHARGE,
                               null_charge_builder.Finish().ValueOrDie());
  QPXValueValidation invalid_key_validator(QPXValueValidation::View::FEATURE);
  const auto invalid_key_result = invalid_key_validator.validate(invalid_keys);
  TEST_FALSE(invalid_key_result.valid)
  TEST_TRUE(invalid_key_result.toString().find("null 'peptidoform'") != std::string::npos)
  TEST_TRUE(invalid_key_result.toString().find("null 'charge'") != std::string::npos)

  // Duplicating the source feature exercises the one-shot writer's validation placement: the
  // output path remains absent because validation runs before FileOutputStream::Open().
  cmap.push_back(cmap.front());
  std::string output;
  NEW_TMP_FILE(output);
  TEST_FALSE(ConsensusMapArrowExport::exportToParquet(cmap, output))
  TEST_FALSE(File::exists(output))

  std::string stream_output;
  NEW_TMP_FILE(stream_output);
  TEST_FALSE(ConsensusMapArrowExport::exportToParquetStreaming(cmap, stream_output, 1))
  TEST_FALSE(File::exists(stream_output))
}
END_SECTION

START_SECTION(([EXTRA] unmapped features are separated by observed_mz, not by rt alone))
{
  // The feature primary key is (peptidoform, charge, run_file_name, rt, observed_mz). For an
  // unmapped feature the peptidoform is empty, so without observed_mz the key would collapse to
  // (charge, run, rt) -- and 'rt' is float32 on write, one ULP being 244 us at 3000 s. Two
  // co-eluting unmapped features of the same charge in one run would then collide and the whole
  // export would be refused. ProteomicsLFQ with seeds leaves most of its features unmapped, so
  // this is the common shape, not an exotic one.
  ConsensusMap cmap = createTestConsensusMap();
  while (!cmap.empty())
  {
    cmap.pop_back();
  }

  const auto make_unmapped = [](double rt, double mz, double intensity)
  {
    ConsensusFeature cf;
    cf.setRT(rt);
    cf.setMZ(mz);
    cf.setIntensity(intensity);
    cf.setCharge(2);
    BaseFeature bf;
    bf.setRT(rt);
    bf.setMZ(mz);
    bf.setIntensity(static_cast<float>(intensity));
    bf.setCharge(2);
    cf.insert(0, bf);
    return cf; // no PeptideIdentification: an unmapped feature
  };

  // Same run, same charge, retention times that differ but round to the SAME float32.
  // Their m/z differ by 100 Th, so they are plainly different features.
  const double rt_a = 3000.0;
  const double rt_b = 3000.0001;
  TEST_EQUAL(static_cast<float>(rt_a), static_cast<float>(rt_b)) // the collision this guards
  cmap.push_back(make_unmapped(rt_a, 400.0, 5000.0));
  cmap.push_back(make_unmapped(rt_b, 500.0, 6000.0));

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  QPXValueValidation validator(QPXValueValidation::View::FEATURE);
  TEST_TRUE(validator.validate(table).valid)

  std::string output;
  NEW_TMP_FILE(output);
  TEST_TRUE(ConsensusMapArrowExport::exportToParquet(cmap, output))
  TEST_TRUE(File::exists(output))

  // The check must still fire when the m/z matches too: then the two rows really are
  // indistinguishable and one of them is a duplicate.
  ConsensusMap same_mz = cmap;
  same_mz.pop_back();
  same_mz.push_back(make_unmapped(rt_b, 400.0, 6000.0));
  auto same_mz_table = ConsensusMapArrowExport::exportToArrow(same_mz);
  TEST_NOT_EQUAL(same_mz_table, nullptr)
  QPXValueValidation same_mz_validator(QPXValueValidation::View::FEATURE);
  const auto same_mz_result = same_mz_validator.validate(same_mz_table);
  TEST_FALSE(same_mz_result.valid)
  TEST_TRUE(same_mz_result.toString().find("primary key") != std::string::npos)
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

  // A row in the QPX feature view is keyed on (feature, run), so the feature needs a handle
  // and a describing column header to produce one at all.
  ConsensusMap::ColumnHeader ch;
  ch.filename = "mods_run.mzML";
  cmap.getColumnHeaders()[0] = ch;

  ConsensusFeature cf;
  cf.setMZ(500.0);
  cf.setRT(100.0);
  cf.setCharge(2);

  BaseFeature bf;
  bf.setIntensity(1234.0f);
  bf.setMZ(500.0);
  bf.setRT(100.0);
  bf.setCharge(2);
  cf.insert(0, bf);

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

START_SECTION(([EXTRA] feature view melts to one row per (feature, run)))
{
  // Label-free: each row carries its own run's coordinates and a single "LFQ" intensity,
  // matching docs/spec/feature.md ("a quantified peptidoform in a specific run file").
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < 2; ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = "/data/run" + StringUtils::toStr(m + 1) + ".mzML";
    ch.label = "label-free";
    cmap.getColumnHeaders()[m] = ch;
  }

  ConsensusFeature cf;
  cf.setMZ(500.0);
  cf.setRT(100.0);
  cf.setCharge(2);
  // Per-run coordinates deliberately differ from the consensus centroid.
  BaseFeature b0; b0.setIntensity(10.0f); b0.setMZ(500.1); b0.setRT(101.0); b0.setCharge(2);
  BaseFeature b1; b1.setIntensity(20.0f); b1.setMZ(500.2); b1.setRT(102.0); b1.setCharge(2);
  cf.insert(0, b0);
  cf.insert(1, b1);
  cmap.push_back(cf);

  // A feature with no handles has no run and no quantification -> no row.
  ConsensusFeature orphan;
  orphan.setMZ(600.0);
  orphan.setRT(200.0);
  cmap.push_back(orphan);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 2)   // one per run; the handleless feature is dropped

  auto run_col = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("run_file_name")->chunk(0));
  auto rt_col  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("rt")->chunk(0));
  auto mz_col  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("observed_mz")->chunk(0));
  TEST_STRING_EQUAL(run_col->GetString(0), "run1")
  TEST_STRING_EQUAL(run_col->GetString(1), "run2")
  // Label-free rows take RT/mz from the run's own handle, not the consensus centroid (100.0/500.0)
  TEST_REAL_SIMILAR(rt_col->Value(0), 101.0)
  TEST_REAL_SIMILAR(rt_col->Value(1), 102.0)
  TEST_REAL_SIMILAR(mz_col->Value(0), 500.1)
  TEST_REAL_SIMILAR(mz_col->Value(1), 500.2)

  // Exactly one intensity per label-free row, labelled "LFQ" and holding that run's value.
  auto int_col = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("intensities")->chunk(0));
  auto st  = std::static_pointer_cast<arrow::StructArray>(int_col->values());
  auto lab = std::static_pointer_cast<arrow::StringArray>(st->field(0));
  auto val = std::static_pointer_cast<arrow::FloatArray>(st->field(1));
  TEST_EQUAL(int_col->value_length(0), 1)
  TEST_EQUAL(int_col->value_length(1), 1)
  TEST_STRING_EQUAL(lab->GetString(0), "LFQ")
  TEST_STRING_EQUAL(lab->GetString(1), "LFQ")
  TEST_REAL_SIMILAR(val->Value(0), 10.0)
  TEST_REAL_SIMILAR(val->Value(1), 20.0)
}
END_SECTION

START_SECTION(([EXTRA] every identification-level column comes from the winning identification))
{
  // A ConsensusFeature can carry several PSMs. If sequence comes from the best-scoring one
  // while scan/ion-mobility/PEP come from pep_ids[0], a row describes one peptide's identity
  // next to another peptide's acquisition metadata.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);

  // First identification: WORSE score (q-value, lower is better), distinct metadata.
  PeptideIdentification worse;
  worse.setScoreType("q-value");
  worse.setHigherScoreBetter(false);
  worse.setSpectrumReference("controllerType=0 controllerNumber=1 scan=111");
  worse.setMetaValue("ion_mobility", 1.11);
  PeptideHit wh;
  wh.setSequence(AASequence::fromString("SHAREDPEPTIDEK"));
  wh.setCharge(2);
  wh.setScore(0.50);
  worse.setHits({wh});

  // Second identification: BETTER score.
  PeptideIdentification better;
  better.setScoreType("q-value");
  better.setHigherScoreBetter(false);
  better.setSpectrumReference("controllerType=0 controllerNumber=1 scan=222");
  better.setMetaValue("ion_mobility", 2.22);
  PeptideHit bh;
  bh.setSequence(AASequence::fromString("SHAREDPEPTIDEK"));
  bh.setCharge(2);
  bh.setScore(0.01);
  better.setHits({bh});

  cf.setPeptideIdentifications({worse, better});
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)

  auto seq  = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("sequence")->chunk(0));
  auto scan = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("scan")->chunk(0));
  auto im   = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("ion_mobility")->chunk(0));

  // Both identifications name the same peptide, so the annotation is unambiguous and the
  // strict identity contract accepts it -- but they carry DIFFERENT acquisition metadata.
  // Every identification-level column must therefore describe the same one of them; reading
  // scan from pep_ids[0] while sequence/PEP come from the winner produced a mixed row.
  TEST_STRING_EQUAL(seq->GetString(0), "SHAREDPEPTIDEK")
  auto scan_vals = std::static_pointer_cast<arrow::Int32Array>(scan->values());
  TEST_EQUAL(scan->value_length(0), 1)
  const int got_scan = scan_vals->Value(scan->value_offset(0));
  const double got_im = im->Value(0);
  // scan and ion mobility must come from the SAME identification, not one from each.
  TEST_TRUE((got_scan == 222 && std::fabs(got_im - 2.22) < 1e-4)
         || (got_scan == 111 && std::fabs(got_im - 1.11) < 1e-4))
}
END_SECTION

START_SECTION(([EXTRA] label-free row coordinates are mutually consistent))
{
  // observed_mz, calculated_mz, mass_error_ppm and the rt bounds must all derive from one
  // source. Computing ppm error against the consensus centroid while reporting the handle's
  // m/z yields an error measured against a value the row does not contain.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ConsensusFeature cf;
  cf.setMZ(999.0);        // consensus centroid, deliberately far from the handle
  cf.setRT(500.0);
  cf.setCharge(2);
  cf.setWidth(0.0);
  BaseFeature bf;
  bf.setIntensity(10.0f);
  bf.setMZ(738.8628);     // ~ the 2+ m/z of PEPTIDEKPEPTIDER
  bf.setRT(101.0);
  bf.setCharge(2);
  bf.setWidth(10.0);
  cf.insert(0, bf);

  PeptideIdentification pid;
  pid.setScoreType("q-value");
  pid.setHigherScoreBetter(false);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEKPEPTIDER"));
  hit.setCharge(2);
  hit.setScore(0.01);
  pid.setHits({hit});
  cf.setPeptideIdentifications({pid});
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)

  auto obs  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("observed_mz")->chunk(0));
  auto calc = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("calculated_mz")->chunk(0));
  auto ppm  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("mass_error_ppm")->chunk(0));
  auto rt   = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("rt")->chunk(0));
  auto rt0  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("rt_start")->chunk(0));
  auto rt1  = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("rt_stop")->chunk(0));

  // Coordinates come from the handle, not the consensus centroid (999.0 / 500.0).
  TEST_REAL_SIMILAR(obs->Value(0), 738.8628)
  TEST_REAL_SIMILAR(rt->Value(0), 101.0)
  // rt bounds derive from the SAME rt and the handle's width (10.0), not the feature's (0.0)
  TEST_REAL_SIMILAR(rt0->Value(0), 96.0)
  TEST_REAL_SIMILAR(rt1->Value(0), 106.0)
  // ppm error must be measured against the m/z the row actually reports ...
  const double expected_ppm = (obs->Value(0) - calc->Value(0)) / calc->Value(0) * 1e6;
  TEST_REAL_SIMILAR(ppm->Value(0), expected_ppm)
  // ... and NOT against the consensus centroid the row does not contain. Stated as a
  // comparison rather than an absolute bound so the test does not depend on the fixture's
  // m/z happening to match the peptide's theoretical mass.
  const double wrong_ppm = (999.0 - calc->Value(0)) / calc->Value(0) * 1e6;
  TEST_TRUE(std::fabs(ppm->Value(0) - wrong_ppm) > 1.0)
}
END_SECTION

// NOTE: a test asserting that identity selection uses ONE score orientation regardless of
// identification order was removed when the strict identity contract landed. It required a
// feature carrying two DIFFERENT peptides to observe which one won, which
// requireUnambiguousIdentities() now rejects outright. Same-peptide identifications differ
// only in score, which the feature view does not export per hit, so the scenario is no longer
// observable through the exporter. The orientation rule itself is still exercised by the
// leading-hitless case below and by the strictness tests in ProteinGroupArrowExport_test.
START_SECTION(([EXTRA] a leading hitless identification does not discard the feature identity))
{
  // The positional pick (pep_ids[0].getHits()[0]) lost a feature's identity whenever its first
  // identification carried no hits, even though a later one did.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);

  PeptideIdentification hitless;
  hitless.setScoreType("q-value");
  hitless.setHigherScoreBetter(false);

  PeptideIdentification real_id;
  real_id.setScoreType("q-value");
  real_id.setHigherScoreBetter(false);
  PeptideHit h; h.setSequence(AASequence::fromString("FOUNDPEPTIDEK")); h.setCharge(2); h.setScore(0.02);
  real_id.setHits({h});

  cf.setPeptideIdentifications({hitless, real_id});
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  auto seq = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("sequence")->chunk(0));
  TEST_STRING_EQUAL(seq->GetString(0), "FOUNDPEPTIDEK")
}
END_SECTION

START_SECTION(([EXTRA] a feature whose identifications are all hitless does not crash))
{
  // pep_ids is non-empty but every identification is hitless, so there is no winning hit.
  // Guarding the scan column on !pep_ids.empty() instead of on the selected identification
  // dereferenced a null pointer here.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/run1.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);

  PeptideIdentification empty_pid;                 // no hits at all
  empty_pid.setScoreType("q-value");
  empty_pid.setHigherScoreBetter(false);
  PeptideIdentification another_empty_pid;
  another_empty_pid.setScoreType("q-value");
  another_empty_pid.setHigherScoreBetter(false);
  cf.setPeptideIdentifications({empty_pid, another_empty_pid});
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)   // quantified, just unidentified
  auto seq = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("sequence")->chunk(0));
  TEST_STRING_EQUAL(seq->GetString(0), "")
}
END_SECTION

START_SECTION(([EXTRA] two runs sharing a stem stay separate rows))
{
  // Fractionated layouts repeat basenames across directories. Grouping melted rows by the
  // stemmed name merged them into one row carrying both runs' intensities -- worse than the
  // ambiguous-but-separate rows the psm exporter warns about.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  for (Size m = 0; m < 2; ++m)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = (m == 0) ? "/frac_a/run1.mzML" : "/frac_b/run1.mzML";  // same stem!
    ch.label = "label-free";
    cmap.getColumnHeaders()[m] = ch;
  }

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature b0; b0.setIntensity(10.0f); b0.setMZ(500.0); b0.setRT(100.0); b0.setCharge(2);
  BaseFeature b1; b1.setIntensity(20.0f); b1.setMZ(500.0); b1.setRT(100.0); b1.setCharge(2);
  cf.insert(0, b0);
  cf.insert(1, b1);
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 2)   // one per source run, not one merged row
  auto ints = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("intensities")->chunk(0));
  TEST_EQUAL(ints->value_length(0), 1)   // 2 would mean the runs were fused
  TEST_EQUAL(ints->value_length(1), 1)
}
END_SECTION

START_SECTION(([EXTRA] isobaric feature rows keep consensus coordinates and list all channels))
{
  // Isobaric handles are reporter ions: their m/z is the reporter mass and their RT the
  // quantification scan, so the row must keep the ConsensusFeature's precursor coordinates
  // and use the handles only for channel intensities.
  ConsensusMap cmap;
  const std::vector<std::string> channels = {"126", "127N", "127C"};
  for (Size c = 0; c < channels.size(); ++c)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = "/data/tmt_run.mzML";      // all channels share one run
    ch.label = "tmt10plex_" + channels[c];
    ch.setMetaValue("channel_name", channels[c]);
    ch.setMetaValue("channel_id", static_cast<int>(c));
    cmap.getColumnHeaders()[c] = ch;
  }

  ConsensusFeature cf;
  cf.setMZ(700.5);   // precursor
  cf.setRT(300.0);
  cf.setCharge(2);
  for (Size c = 0; c < channels.size(); ++c)
  {
    BaseFeature b;
    b.setIntensity(100.0f * (c + 1));
    b.setMZ(126.0 + c);      // reporter mass, must NOT become observed_mz
    b.setRT(301.0);          // quant scan RT, must NOT become rt
    cf.insert(c, b);
  }
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  TEST_EQUAL(t->num_rows(), 1)   // one run -> one row, with all channels in it

  auto rt_col = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("rt")->chunk(0));
  auto mz_col = std::static_pointer_cast<arrow::FloatArray>(t->GetColumnByName("observed_mz")->chunk(0));
  TEST_REAL_SIMILAR(rt_col->Value(0), 300.0)   // consensus, not 301.0
  TEST_REAL_SIMILAR(mz_col->Value(0), 700.5)   // precursor, not a reporter mass

  auto int_col = std::static_pointer_cast<arrow::ListArray>(t->GetColumnByName("intensities")->chunk(0));
  auto st  = std::static_pointer_cast<arrow::StructArray>(int_col->values());
  auto lab = std::static_pointer_cast<arrow::StringArray>(st->field(0));
  TEST_EQUAL(int_col->value_length(0), 3)
  TEST_STRING_EQUAL(lab->GetString(0), "TMT126")
  TEST_STRING_EQUAL(lab->GetString(1), "TMT127N")
  TEST_STRING_EQUAL(lab->GetString(2), "TMT127C")
}
END_SECTION

START_SECTION(([EXTRA] SILAC feature labels use the canonical QPX two-plex vocabulary))
{
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS1");
  const std::vector<std::string> raw_labels = {"no_label", "Arg6"};
  for (Size channel = 0; channel < raw_labels.size(); ++channel)
  {
    ConsensusMap::ColumnHeader header;
    header.filename = "/data/silac_run.mzML";
    header.label = raw_labels[channel];
    header.setMetaValue("channel_id", static_cast<int>(channel));
    cmap.getColumnHeaders()[channel] = header;
  }

  ConsensusFeature feature;
  feature.setMZ(600.0);
  feature.setRT(120.0);
  feature.setCharge(2);
  for (Size channel = 0; channel < raw_labels.size(); ++channel)
  {
    BaseFeature handle;
    handle.setIntensity(100.0f * (channel + 1));
    handle.setMZ(600.0);
    handle.setRT(120.0);
    handle.setCharge(2);
    feature.insert(channel, handle);
  }
  cmap.push_back(feature);

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  auto intensities = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("intensities")->chunk(0));
  auto entries = std::static_pointer_cast<arrow::StructArray>(intensities->values());
  auto labels = std::static_pointer_cast<arrow::StringArray>(entries->field(0));
  TEST_EQUAL(intensities->value_length(0), 2)
  TEST_STRING_EQUAL(labels->GetString(0), "SILAC light")
  // Arg6 is the second/heavy-role channel in a two-plex, not "SILAC medium".
  TEST_STRING_EQUAL(labels->GetString(1), "SILAC heavy")

  QPXValueValidation validator(QPXValueValidation::View::FEATURE);
  TEST_TRUE(validator.validate(table).valid)
}
END_SECTION

START_SECTION(exportToParquet - basic export)
{
  ConsensusMap cmap = createTestConsensusMap();

  std::string filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename);
  TEST_TRUE(success)

  // Verify file was created
  TEST_TRUE(File::exists(filename))
  TEST_FALSE(File::empty(filename))
}
END_SECTION

START_SECTION(exportToParquet - empty consensus map)
{
  ConsensusMap cmap;

  std::string filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = ConsensusMapArrowExport::exportToParquet(cmap, filename);
  TEST_TRUE(success)

  // Verify file was created (even for empty data)
  TEST_TRUE(File::exists(filename))
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
  TEST_TRUE(success)

  // Verify file was created
  TEST_TRUE(File::exists(filename))
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
  TEST_TRUE(success)

  TEST_TRUE(File::exists(filename))
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
  TEST_TRUE(success)

  TEST_TRUE(File::exists(filename))
}
END_SECTION

START_SECTION(([EXTRA] exportToParquet leaves no file behind when the write fails))
{
  // The value validation runs before the file is opened, so a refused table never creates one.
  // A write that fails afterwards is a different matter: FileOutputStream::Open has created and
  // truncated the file and the Parquet writer has already put its magic bytes there, so what
  // survives is a headerless fragment that reports as corrupt rather than as the smaller table
  // it looks like. The streaming writer handles this; the one-shot writer must too.
  // A row group size of 0 is refused by Parquet for a non-empty table, which reaches that
  // failure deterministically on every platform.
  ConsensusMap cmap = createTestConsensusMap();

  std::string filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  ParquetWriteConfig pq_config;
  pq_config.row_group_size = 0;

  TEST_FALSE(ConsensusMapArrowExport::exportToParquet(cmap, filename, pq_config))
  TEST_FALSE(File::exists(filename))
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

  // --- One-shot reference ---
  std::string ref_file; NEW_TMP_FILE(ref_file)
  TEST_TRUE(ConsensusMapArrowExport::exportToParquet(cmap, ref_file))
  auto ref_table = read_combined(ref_file);
  TEST_NOT_EQUAL(ref_table, nullptr)
  TEST_EQUAL(ref_table->num_rows(), (int64_t)N)

  // --- Streaming must equal the one-shot output for ANY thread count (deterministic build).
  //     Small batch (1000) forces multiple row groups; n_threads exercises serial (1), fixed
  //     parallel (2, 8), and the production auto path (0 = omp_get_max_threads, used by IsobaricWorkflow). ---
  for (int nt : {0, 1, 2, 8})
  {
    std::string stream_file; NEW_TMP_FILE(stream_file)
    TEST_TRUE(ConsensusMapArrowExport::exportToParquetStreaming(cmap, stream_file, 1000, ParquetWriteConfig{}, nt))
    auto st_table = read_combined(stream_file);
    TEST_NOT_EQUAL(st_table, nullptr)
    TEST_EQUAL(st_table->num_rows(), (int64_t)N)
    TEST_EQUAL(st_table->num_columns(), ref_table->num_columns())

    // Order-sensitive scalar check (unique rt) — catches any partition reordering under threads.
    auto seq_s = std::static_pointer_cast<arrow::StringArray>(st_table->GetColumnByName("sequence")->chunk(0));
    auto seq_r = std::static_pointer_cast<arrow::StringArray>(ref_table->GetColumnByName("sequence")->chunk(0));
    auto rt_s  = std::static_pointer_cast<arrow::FloatArray>(st_table->GetColumnByName("rt")->chunk(0));
    auto rt_r  = std::static_pointer_cast<arrow::FloatArray>(ref_table->GetColumnByName("rt")->chunk(0));
    auto dec_s = std::static_pointer_cast<arrow::BooleanArray>(st_table->GetColumnByName("is_decoy")->chunk(0));
    auto dec_r = std::static_pointer_cast<arrow::BooleanArray>(ref_table->GetColumnByName("is_decoy")->chunk(0));
    bool all_eq = true;
    for (int64_t r = 0; r < st_table->num_rows(); ++r)
    {
      if (seq_s->GetString(r) != seq_r->GetString(r)) { all_eq = false; break; }
      if (rt_s->Value(r)      != rt_r->Value(r))      { all_eq = false; break; } // order-sensitive
      if (dec_s->Value(r)     != dec_r->Value(r))     { all_eq = false; break; }
    }
    TEST_TRUE(all_eq)

    // Full logical equivalence: also covers nested/list columns (modifications, intensities,
    // pg_accessions, gg_*) and lookup-derived columns (pg_global_qvalue). Both tables were
    // CombineChunks()ed and round-tripped through parquet the same way, so schema (incl. nested
    // types) and values must match; ignore schema metadata.
    TEST_TRUE(st_table->Equals(*ref_table))
  }

  // --- batch_size >= N: single batch built across all threads, still correct ---
  {
    std::string one_batch; NEW_TMP_FILE(one_batch)
    TEST_TRUE(ConsensusMapArrowExport::exportToParquetStreaming(cmap, one_batch, N + 1000, ParquetWriteConfig{}, 8))
    auto t = read_combined(one_batch);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), (int64_t)N)
    TEST_TRUE(t->Equals(*ref_table))
  }

  // --- batch_size == 0 guard: treated as default, valid file ---
  {
    std::string zero_batch; NEW_TMP_FILE(zero_batch)
    TEST_TRUE(ConsensusMapArrowExport::exportToParquetStreaming(cmap, zero_batch, 0))
    auto t = read_combined(zero_batch);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), (int64_t)N)
  }

  // --- Edge case: empty map -> valid 0-row file with full schema ---
  {
    ConsensusMap empty;
    std::string empty_file; NEW_TMP_FILE(empty_file)
    TEST_TRUE(ConsensusMapArrowExport::exportToParquetStreaming(empty, empty_file, 1000))
    TEST_TRUE(File::exists(empty_file))
    auto t = read_combined(empty_file);
    TEST_NOT_EQUAL(t, nullptr)
    TEST_EQUAL(t->num_rows(), 0)
    TEST_TRUE(t->num_columns() > 0)
  }
}
END_SECTION

START_SECTION([EXTRA] exportToArrow - dedicated metavalue columns resolve to correct values)
{
  // Guards the index-based (lock-free) metavalue reads in buildFeatureTableRange: one identified
  // feature with every dedicated metavalue set to a known value; assert each exported column. The
  // streaming<->one-shot equivalence test cannot catch an index-name typo here, since both paths
  // now share buildFeatureTableRange.
  ConsensusMap cmap;
  ConsensusMap::ColumnHeaders headers;
  ConsensusMap::ColumnHeader ch0; ch0.filename = "run.mzML"; headers[0] = ch0;
  cmap.setColumnHeaders(headers);

  ProteinIdentification prot_id; prot_id.setIdentifier("PI_0");
  ProteinHit ph; ph.setAccession("PROT_A"); ph.setScore(0.9);
  prot_id.setHits({ph});
  ProteinIdentification::ProteinGroup pgrp; pgrp.probability = 0.02; pgrp.accessions = {"PROT_A"};
  prot_id.insertProteinGroup(pgrp);
  cmap.setProteinIdentifications({prot_id});

  ConsensusFeature cf;
  cf.setMZ(555.5); cf.setRT(321.0); cf.setCharge(2);
  cf.setMetaValue("start_ion_mobility", 0.80);
  cf.setMetaValue("stop_ion_mobility", 0.90);
  cf.setMetaValue("rt_start", 300.0);
  cf.setMetaValue("rt_stop", 340.0);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(555.5); bf.setRT(321.0);
  cf.insert(0, bf);

  PeptideIdentification pid;
  pid.setScoreType("Posterior Error Probability"); pid.setHigherScoreBetter(false);
  pid.setMetaValue("ion_mobility", 0.85);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  hit.setScore(0.02);                       // main score -> posterior_error_probability
  hit.setMetaValue("target_decoy", "decoy");
  hit.setMetaValue("predicted_RT", 319.0);
  hit.setMetaValue("missed_cleavages", 1);
  hit.setMetaValue("q-value", 0.05);        // known score (QVAL, lower-better) -> additional_scores
  hit.setMetaValue("zzz_custom_flag", 7);   // int, but NOT a known score -> filtered out
  hit.setMetaValue("zzz_custom_note", "x"); // string -> filtered out (not INT/DOUBLE)
  PeptideEvidence ev; ev.setProteinAccession("PROT_A");
  hit.setPeptideEvidences({ev});
  pid.setHits({hit});
  cf.setPeptideIdentifications({pid});
  cmap.push_back(cf);

  // Second feature exercises the FALLBACK precedence branches: predicted_rt (lowercase, only reached
  // when predicted_RT is absent) and IM (only reached when ion_mobility is absent). A predicted_RT<->
  // predicted_rt or ion_mobility<->IM index-name swap would still pass row 0 but fails here.
  ConsensusFeature cf2;
  cf2.setMZ(600.0); cf2.setRT(400.0); cf2.setCharge(2);
  BaseFeature bf2; bf2.setIntensity(5.0f); bf2.setMZ(600.0); bf2.setRT(400.0);
  cf2.insert(0, bf2);
  PeptideIdentification pid2;
  pid2.setScoreType("Posterior Error Probability"); pid2.setHigherScoreBetter(false);
  pid2.setMetaValue("IM", 1.10);            // fallback path (no "ion_mobility")
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("PEPTIDEK"));
  hit2.setCharge(2); hit2.setScore(0.03);
  hit2.setMetaValue("predicted_rt", 405.0); // fallback path (lowercase, no "predicted_RT")
  PeptideEvidence ev2; ev2.setProteinAccession("PROT_A");
  hit2.setPeptideEvidences({ev2});
  pid2.setHits({hit2});
  cf2.setPeptideIdentifications({pid2});
  cmap.push_back(cf2);

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  auto f32 = [&](const char* n){ return std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName(n)->chunk(0)); };
  auto is_decoy = std::static_pointer_cast<arrow::BooleanArray>(table->GetColumnByName("is_decoy")->chunk(0));
  auto pep      = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("posterior_error_probability")->chunk(0));
  auto mc       = std::static_pointer_cast<arrow::Int16Array>(table->GetColumnByName("missed_cleavages")->chunk(0));
  auto pgq      = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("pg_global_qvalue")->chunk(0));
  auto as_col   = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("additional_scores")->chunk(0));

  TEST_TRUE(is_decoy->Value(0))
  TEST_FALSE(pep->IsNull(0))
  TEST_REAL_SIMILAR(pep->Value(0), 0.02)
  TEST_FALSE(mc->IsNull(0))
  TEST_EQUAL(mc->Value(0), 1)
  TEST_FALSE(pgq->IsNull(0))
  TEST_REAL_SIMILAR(pgq->Value(0), 0.02)

  TEST_FALSE(f32("predicted_rt")->IsNull(0))
  TEST_REAL_SIMILAR(f32("predicted_rt")->Value(0), 319.0)
  TEST_FALSE(f32("ion_mobility")->IsNull(0))
  TEST_REAL_SIMILAR(f32("ion_mobility")->Value(0), 0.85)
  TEST_FALSE(f32("ion_mobility_start")->IsNull(0))
  TEST_REAL_SIMILAR(f32("ion_mobility_start")->Value(0), 0.80)
  TEST_FALSE(f32("ion_mobility_stop")->IsNull(0))
  TEST_REAL_SIMILAR(f32("ion_mobility_stop")->Value(0), 0.90)
  TEST_FALSE(f32("rt_start")->IsNull(0))
  TEST_REAL_SIMILAR(f32("rt_start")->Value(0), 300.0)
  TEST_FALSE(f32("rt_stop")->IsNull(0))
  TEST_REAL_SIMILAR(f32("rt_stop")->Value(0), 340.0)

  // Row 1: fallback precedence (predicted_rt lowercase, IM) must resolve to the right column/value.
  TEST_FALSE(f32("predicted_rt")->IsNull(1))
  TEST_REAL_SIMILAR(f32("predicted_rt")->Value(1), 405.0)
  TEST_FALSE(f32("ion_mobility")->IsNull(1))
  TEST_REAL_SIMILAR(f32("ion_mobility")->Value(1), 1.10)

  // additional_scores: only the known "q-value" is emitted via the index-based metaBegin pass; the
  // non-score int and the string metavalue are filtered out. Pin the exact struct fields (name, value,
  // and higher_better = false for a q-value, which the score-type switcher maps to QVAL).
  TEST_EQUAL(as_col->value_length(0), 1)
  {
    auto structs = std::static_pointer_cast<arrow::StructArray>(as_col->values());
    auto sname = std::static_pointer_cast<arrow::StringArray>(structs->GetFieldByName("score_name"));
    auto sval  = std::static_pointer_cast<arrow::DoubleArray>(structs->GetFieldByName("score_value"));
    auto shb   = std::static_pointer_cast<arrow::BooleanArray>(structs->GetFieldByName("higher_better"));
    const int64_t off = as_col->value_offset(0);
    TEST_EQUAL(sname->GetString(off), "q-value")
    TEST_REAL_SIMILAR(sval->Value(off), 0.05)
    TEST_FALSE(shb->IsNull(off))
    TEST_FALSE(shb->Value(off))
  }
}
END_SECTION

START_SECTION([EXTRA] exportToArrow - run_file_name is emitted without path or extension)
{
  // QPX defines run_file_name as the spectrum file name without path or extension. The feature,
  // psm and pg tables each derive it from a different source, so all three must stem it for the
  // cross-table join (and the Hive partitioning) to work.
  ConsensusMap cmap = createTestConsensusMap();
  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "sample1") // ColumnHeader::filename is "sample1.mzML"
}
END_SECTION

START_SECTION([EXTRA] exportToArrow - isobaric channels of one file share one run_file_name)
{
  // For isobaric input every channel of a file gets its own ColumnHeader with the same `filename`;
  // the channel identity lives in `label`. So whichever channel column the first FeatureHandle
  // happens to point at, the row must resolve to that file's single stem.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");

  ConsensusMap::ColumnHeaders headers;
  const StringList labels = {"tmt6plex_126", "tmt6plex_127", "tmt6plex_128"};
  for (Size i = 0; i < labels.size(); ++i)
  {
    ConsensusMap::ColumnHeader ch;
    ch.filename = "/raw/20150820_Haura-Pilot-TMT1-bRPLC01-1.mzML";
    ch.label = labels[i];
    headers[i] = ch;
  }
  cmap.setColumnHeaders(headers);

  ProteinIdentification prot_id;
  prot_id.setIdentifier("PI_0");
  cmap.setProteinIdentifications({prot_id});

  // Two features anchored at different channel columns (2 and 0)
  for (UInt64 first_map : {UInt64(2), UInt64(0)})
  {
    ConsensusFeature cf;
    cf.setMZ(500.25);
    cf.setRT(100.5);
    cf.setCharge(2);
    for (UInt64 m : {first_map, UInt64(1)})
    {
      BaseFeature bf;
      bf.setIntensity(1000.0f);
      bf.setMZ(500.25);
      bf.setRT(100.5);
      cf.insert(m, bf);
    }
    cmap.push_back(cf);
  }

  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)
  auto rfn = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(rfn->GetString(0), "20150820_Haura-Pilot-TMT1-bRPLC01-1")
  TEST_STRING_EQUAL(rfn->GetString(1), "20150820_Haura-Pilot-TMT1-bRPLC01-1")
}
END_SECTION

START_SECTION(([EXTRA] a merged run without id_merge_index is refused, not resolved to file 0))
{
  // id_run_file_name is resolved through IdentifierMSRunMapper, which falls back to the run's
  // FIRST file when the index is missing -- so the column would name a run the best PSM did not
  // come from. The psm and pg views refuse this; the feature view must too, and it matters that
  // it does: ProteomicsLFQ writes the feature view FIRST and pg third, so accepting here and
  // refusing there leaves a wrong feature.parquet on disk beside a half-written collection.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/runA.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath({"/data/runA.mzML", "/data/runB.mzML"}); // merged run
  cmap.setProteinIdentifications({prot});

  PeptideIdentification pid;                 // no id_merge_index: origin undetermined
  pid.setIdentifier("merged");
  pid.setScoreType("q-value");
  pid.setHigherScoreBetter(false);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  pid.setHits({hit});

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);
  cf.setPeptideIdentifications({pid});
  cmap.push_back(cf);

  TEST_EXCEPTION(Exception::MissingInformation, ConsensusMapArrowExport::exportToArrow(cmap))

  // With the index present the same input exports and names the right run.
  cmap[0].getPeptideIdentifications()[0].setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 1);
  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  auto idrun = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("id_run_file_name")->chunk(0));
  TEST_STRING_EQUAL(idrun->GetString(0), "runB")
}
END_SECTION

START_SECTION(([EXTRA] only the selected identification must be resolvable))
{
  // The row draws every identification-derived value from the winning hit, so a sibling
  // identification that carries no hit cannot affect the output and must not be able to refuse
  // the export. Validating every attached identification rejected input that exports correctly.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/runA.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setHigherScoreBetter(false);
  prot.setPrimaryMSRunPath({"/data/runA.mzML", "/data/runB.mzML"});
  cmap.setProteinIdentifications({prot});

  PeptideIdentification hitless;              // no hits, no id_merge_index: never selected
  hitless.setIdentifier("merged");
  hitless.setScoreType("q-value");
  hitless.setHigherScoreBetter(false);

  PeptideIdentification winner;               // carries the hit, and a usable index
  winner.setIdentifier("merged");
  winner.setScoreType("q-value");
  winner.setHigherScoreBetter(false);
  winner.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 1);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  winner.setHits({hit});

  ConsensusFeature cf;
  cf.setMZ(500.0); cf.setRT(100.0); cf.setCharge(2);
  BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0); bf.setCharge(2);
  cf.insert(0, bf);
  cf.setPeptideIdentifications({hitless, winner});
  cmap.push_back(cf);

  auto t = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(t, nullptr)
  auto idrun = std::static_pointer_cast<arrow::StringArray>(t->GetColumnByName("id_run_file_name")->chunk(0));
  TEST_STRING_EQUAL(idrun->GetString(0), "runB")
}
END_SECTION

START_SECTION(([EXTRA] duplicate identification-run identifiers are refused))
{
  // IdentifierMSRunMapper assigns its forward map with operator[], so two runs sharing an
  // identifier leave only the last one's paths -- and a feature of the first would then be
  // stamped with a file it never came from, with one row and no warning to give it away.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/runA.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification a;
  a.setIdentifier("dup");
  a.setPrimaryMSRunPath({"/data/runA.mzML"});
  ProteinIdentification b;
  b.setIdentifier("dup");                     // same identifier, different file
  b.setPrimaryMSRunPath({"/data/runB.mzML"});
  cmap.setProteinIdentifications({a, b});

  TEST_EXCEPTION(Exception::InvalidValue, ConsensusMapArrowExport::exportToArrow(cmap))
}
END_SECTION

START_SECTION(([EXTRA] feature_id is derived from the row this exporter wrote, in the declared order))
{
  // Re-derives every row's feature_id from the columns AS WRITTEN, passing them in the order
  // QPXIdentity::FEATURE_COMPOSITE declares. QPXIdentity_test pins the derivation itself against
  // the reference implementation's vectors; what this pins is the exporter's use of it -- an
  // exporter that hashed the run name where the peptidoform belongs, or an rt it had not yet
  // narrowed to float32, would produce ids no reader could reproduce, and neither the schema
  // check nor a value reference would notice.
  ConsensusMap cmap = createTestConsensusMap();
  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)

  auto ids  = std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName("feature_id")->chunk(0));
  auto run  = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  auto pf   = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("peptidoform")->chunk(0));
  auto chg  = std::static_pointer_cast<arrow::Int16Array>(table->GetColumnByName("charge")->chunk(0));
  auto rt   = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("rt")->chunk(0));
  auto mz   = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("observed_mz")->chunk(0));
  auto scan = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("scan")->chunk(0));
  auto scan_values = std::static_pointer_cast<arrow::Int32Array>(scan->values());

  std::set<Int64> seen;
  for (int64_t r = 0; r < table->num_rows(); ++r)
  {
    std::vector<Int32> components;
    for (int64_t i = 0; i < scan->value_length(r); ++i)
    {
      components.push_back(scan_values->Value(scan->value_offset(r) + i));
    }
    const std::optional<float> row_rt =
      rt->IsNull(r) ? std::optional<float>{} : std::optional<float>{rt->Value(r)};

    TEST_FALSE(ids->IsNull(r))
    TEST_EQUAL(ids->Value(r), QPXIdentity::featureId(run->GetString(r), pf->GetString(r),
                                                     chg->Value(r), row_rt, components, mz->Value(r)))
    // The identity is also the view's primary key.
    TEST_TRUE(seen.insert(ids->Value(r)).second)
  }
}
END_SECTION

START_SECTION(([EXTRA] the feature<->PSM cross-references resolve and agree in both directions))
{
  // qpx validates three things about these columns: that psm.feature_id resolves to a feature,
  // that feature.psm_ids resolve to PSMs, and that the two directions describe the SAME edge
  // (its "reciprocal desync" check). Producing both from one pass is what makes that hold, so
  // this asserts the property rather than the implementation detail.
  ConsensusMap cmap = createTestConsensusMap();
  // The shared fixture leaves the identification run's MS run path unset and its PSMs' run
  // identifier unlinked. Both are legitimate, but together they mean the psm view cannot resolve
  // run_file_name -- and a PSM with no run matches no feature row, so the cross-references would
  // come out empty and every assertion below would hold vacuously. The `total_listed > 0` guard
  // further down is what catches that if these two lines are ever lost.
  cmap.getProteinIdentifications()[0].setPrimaryMSRunPath({"sample1.mzML"});
  for (auto& cf : cmap)
  {
    for (auto& pid : cf.getPeptideIdentifications()) { pid.setIdentifier("PI_0"); }
  }

  QPXIdentity::FeatureLinks links;
  auto features = ConsensusMapArrowExport::exportToArrow(cmap, &links);
  TEST_NOT_EQUAL(features, nullptr)

  PeptideIdentificationList all_pepids;
  for (const auto& cf : cmap)
  {
    for (const auto& pid : cf.getPeptideIdentifications()) { all_pepids.push_back(pid); }
  }
  for (const auto& pid : cmap.getUnassignedPeptideIdentifications()) { all_pepids.push_back(pid); }
  auto psms = QPXFile::exportPSMsToQPXArrow(cmap.getProteinIdentifications(), all_pepids,
                                            /*export_all_psms=*/false, &links);
  TEST_NOT_EQUAL(psms, nullptr)

  auto feature_ids = std::static_pointer_cast<arrow::Int64Array>(features->GetColumnByName("feature_id")->chunk(0));
  auto psm_ids_col = std::static_pointer_cast<arrow::ListArray>(features->GetColumnByName("psm_ids")->chunk(0));
  auto psm_ids_values = std::static_pointer_cast<arrow::Int64Array>(psm_ids_col->values());
  auto psm_id  = std::static_pointer_cast<arrow::Int64Array>(psms->GetColumnByName("psm_id")->chunk(0));
  auto psm_feature_id = std::static_pointer_cast<arrow::Int64Array>(psms->GetColumnByName("feature_id")->chunk(0));

  std::set<Int64> known_features;
  for (int64_t r = 0; r < features->num_rows(); ++r) { known_features.insert(feature_ids->Value(r)); }
  std::set<Int64> known_psms;
  for (int64_t r = 0; r < psms->num_rows(); ++r) { known_psms.insert(psm_id->Value(r)); }

  // feature.psm_ids -> psm.psm_id, and the reverse edge each one implies
  std::map<Int64, Int64> listed_by;    // psm_id -> the feature that lists it
  Size total_listed = 0;
  for (int64_t r = 0; r < features->num_rows(); ++r)
  {
    for (int64_t i = 0; i < psm_ids_col->value_length(r); ++i)
    {
      const Int64 referenced = psm_ids_values->Value(psm_ids_col->value_offset(r) + i);
      TEST_EQUAL(known_psms.count(referenced), 1)   // no dangling psm_id
      listed_by[referenced] = feature_ids->Value(r);
      ++total_listed;
    }
  }
  TEST_TRUE(total_listed > 0)   // an assertion over an empty set proves nothing

  // psm.feature_id -> feature.feature_id, and it must be the feature that lists it back
  Size total_linked = 0;
  for (int64_t r = 0; r < psms->num_rows(); ++r)
  {
    if (psm_feature_id->IsNull(r)) { continue; }   // unassigned identification: correctly unlinked
    TEST_EQUAL(known_features.count(psm_feature_id->Value(r)), 1)  // no dangling feature_id
    const auto back = listed_by.find(psm_id->Value(r));
    TEST_TRUE(back != listed_by.end())
    if (back != listed_by.end()) { TEST_EQUAL(back->second, psm_feature_id->Value(r)) }
    ++total_linked;
  }
  TEST_EQUAL(total_linked, total_listed)   // the two directions describe exactly one edge set
}
END_SECTION

START_SECTION(([EXTRA] a psm view exported without a feature view leaves feature_id null))
{
  // The spec's answer for an unlinked PSM is null, not a fabricated reference. Every psm-only
  // producer (ProSE, any search engine) takes this path.
  ConsensusMap cmap = createTestConsensusMap();
  PeptideIdentificationList all_pepids;
  for (const auto& cf : cmap)
  {
    for (const auto& pid : cf.getPeptideIdentifications()) { all_pepids.push_back(pid); }
  }
  auto psms = QPXFile::exportPSMsToQPXArrow(cmap.getProteinIdentifications(), all_pepids);
  TEST_NOT_EQUAL(psms, nullptr)
  TEST_TRUE(psms->num_rows() > 0)
  auto psm_feature_id = std::static_pointer_cast<arrow::Int64Array>(psms->GetColumnByName("feature_id")->chunk(0));
  TEST_EQUAL(psm_feature_id->null_count(), psms->num_rows())
}
END_SECTION

START_SECTION(([EXTRA] two feature rows cannot deliver a PSM claimed by both))
{
  // Both directions of the cross-reference are built from one loop, so they agree per row. The
  // remaining way they could disagree is ACROSS rows: two consensus features in one run whose
  // identifications share all four psm-key columns derive the same psm_id, so both rows list it
  // while psm.feature_id can name only one of them -- the "reciprocal desync" qpx checks for.
  //
  // That state cannot reach disk, and this pins the reason: the same two identifications are
  // exactly a repeated psm primary key, which the psm view refuses. The feature view is written
  // first, so the refusal arrives after feature.parquet exists -- which is what QPXCollectionExport's
  // transaction is for, and why this is a refusal rather than a silently inconsistent collection.
  ConsensusMap cmap;
  cmap.setExperimentType("label-free");
  ConsensusMap::ColumnHeader ch;
  ch.filename = "/data/runA.mzML";
  ch.label = "label-free";
  cmap.getColumnHeaders()[0] = ch;

  ProteinIdentification prot;
  prot.setIdentifier("PI_0");
  prot.setPrimaryMSRunPath({"/data/runA.mzML"});
  ProteinIdentification::ProteinGroup group;
  group.probability = 0.01;
  group.accessions = {"PROT_A"};
  prot.insertProteinGroup(group);
  ProteinHit ph;
  ph.setAccession("PROT_A");
  prot.setHits({ph});
  cmap.setProteinIdentifications({prot});

  // One identification, two consensus features. Same run, scan, peptidoform and charge -> same
  // psm_id; different rt and m/z -> two distinct feature rows, both listing it.
  const auto make_id = []()
  {
    PeptideIdentification pid;
    pid.setIdentifier("PI_0");
    pid.setSpectrumReference("scan=1234");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    hit.setScore(0.01);
    pid.setHits({hit});
    return pid;
  };

  PeptideIdentificationList all_pepids;
  for (int i = 0; i < 2; ++i)
  {
    ConsensusFeature cf;
    cf.setMZ(500.0 + i);
    cf.setRT(100.0 + i);
    cf.setCharge(2);
    BaseFeature bf;
    bf.setIntensity(1000.0f);
    bf.setMZ(500.0 + i);
    bf.setRT(100.0 + i);
    bf.setCharge(2);
    cf.insert(0, bf);
    cf.setPeptideIdentifications({make_id()});
    cmap.push_back(cf);
    all_pepids.push_back(make_id());
  }

  QPXIdentity::FeatureLinks links;
  auto features = ConsensusMapArrowExport::exportToArrow(cmap, &links);
  TEST_NOT_EQUAL(features, nullptr)
  TEST_EQUAL(features->num_rows(), 2)

  // The premise: the two rows really do claim one psm_id, so this is not vacuous.
  auto psm_ids_col = std::static_pointer_cast<arrow::ListArray>(features->GetColumnByName("psm_ids")->chunk(0));
  auto psm_ids_values = std::static_pointer_cast<arrow::Int64Array>(psm_ids_col->values());
  TEST_EQUAL(psm_ids_col->value_length(0), 1)
  TEST_EQUAL(psm_ids_col->value_length(1), 1)
  TEST_EQUAL(psm_ids_values->Value(psm_ids_col->value_offset(0)),
             psm_ids_values->Value(psm_ids_col->value_offset(1)))
  TEST_EQUAL(links.size(), 1)   // one psm_id, so one entry however the two writes are resolved

  // And the refusal: the psm view sees the same two identifications as one repeated key.
  auto psms = QPXFile::exportPSMsToQPXArrow(cmap.getProteinIdentifications(), all_pepids,
                                            /*export_all_psms=*/false, &links);
  TEST_NOT_EQUAL(psms, nullptr)
  TEST_EQUAL(psms->num_rows(), 2)
  QPXValueValidation psm_validator(QPXValueValidation::View::PSM);
  const auto verdict = psm_validator.validate(psms);
  TEST_FALSE(verdict.valid)
  TEST_TRUE(verdict.toString().find("repeats the QPX psm identity") != std::string::npos)
}
END_SECTION

END_TEST
