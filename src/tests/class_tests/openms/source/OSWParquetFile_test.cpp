// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/OSWParquetFile.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <fstream>
#include <iterator>
#include <optional>
#include <utility>

using namespace OpenMS;

namespace
{
template<typename Builder, typename Value>
std::shared_ptr<arrow::Array> array_(std::initializer_list<Value> values)
{
  Builder builder;
  for (const auto& value : values)
  {
    ParquetFile::appendOrThrow(builder.Append(value), "fixture");
  }
  return ParquetFile::finishArray(builder, "fixture");
}

std::shared_ptr<arrow::Array> doubles_(std::initializer_list<std::optional<double>> values)
{
  arrow::DoubleBuilder builder;
  for (const auto& value : values)
  {
    ParquetFile::appendOrThrow(value ? builder.Append(*value) : builder.AppendNull(), "fixture");
  }
  return ParquetFile::finishArray(builder, "fixture");
}

void table_(const std::string& filename, const std::vector<std::pair<std::string, std::shared_ptr<arrow::Array>>>& columns)
{
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;
  for (const auto& [name, values] : columns)
  {
    fields.push_back(arrow::field(name, values->type()));
    arrays.push_back(values);
  }
  ParquetFile::writeTable(arrow::Table::Make(arrow::schema(fields), arrays), filename);
}

OpenSwath::LightTargetedExperiment library_()
{
  OpenSwath::LightTargetedExperiment library;
  for (int i = 0; i < 2; ++i)
  {
    OpenSwath::LightProtein protein;
    protein.id = i == 0 ? "P_TARGET" : "P_DECOY";
    library.proteins.push_back(protein);
    OpenSwath::LightCompound compound;
    compound.id = i == 0 ? "0" : "7";
    compound.sequence = i == 0 ? "PEPTIDEA" : "PEPTIDEC";
    compound.charge = 2;
    compound.rt = 100.0;
    compound.protein_refs = {protein.id};
    compound.gene_name = i == 0 ? "GENE_TARGET" : "GENE_DECOY";
    library.compounds.push_back(compound);
    OpenSwath::LightTransition transition;
    transition.transition_name = i == 0 ? "3" : "100";
    transition.peptide_ref = compound.id;
    transition.precursor_mz = 500.0 + i;
    transition.product_mz = 200.0 + i;
    transition.fragment_charge = 1;
    transition.fragment_nr = 3;
    transition.library_intensity = 1000.0;
    transition.setDetectingTransition(true);
    transition.setDecoy(i == 1);
    transition.setFragmentType(i == 0 ? "y" : "b");
    library.transitions.push_back(transition);
  }
  return library;
}

void fixture_(const std::string& directory, bool scored = true)
{
  File::makeDir(directory + "/runs");
  table_(directory + "/runs/runs.parquet",
         {{"run_id", array_<arrow::Int64Builder>({11, 22})}, {"filename", array_<arrow::StringBuilder>({"sample one.mzML.gz", "sample two.mzML"})}});
  for (const Int64 run : {11, 22})
  {
    const std::string base = directory + "/runs/run_id=" + StringUtils::toStr(run);
    File::makeDir(base);
    std::vector<std::pair<std::string, std::shared_ptr<arrow::Array>>> columns
      = {{"feature_id", array_<arrow::Int64Builder>({run * 10, run * 10 + 1, run * 10 + 2})},
         {"precursor_id", array_<arrow::Int64Builder>({0, 0, 7})},
         {"exp_rt", doubles_({99.0, 101.0, 100.0})},
         {"ms2_area_intensity", doubles_({100.0, 200.0, 300.0})},
         {"score_ms2_qvalue", run == 11 ? doubles_({0.01, 0.02, std::nullopt}) : doubles_({0.01, 0.02, 0.03})},
         {"score_ms2_pep", doubles_({0.02, std::nullopt, 0.05})},
         {"unrelated_column", array_<arrow::StringBuilder>({"keep A", "keep B", "keep C"})}};
    if (scored) { columns.emplace_back("score_ms2_score", run == 11 ? doubles_({2.0, 5.0, std::nullopt}) : doubles_({6.0, 4.0, 3.0})); }
    table_(base + "/features.parquet", columns);
    table_(base + "/feature_transition.parquet", {{"feature_id", array_<arrow::Int64Builder>({run * 10, run * 10 + 1, run * 10 + 2})},
                                                  {"transition_id", array_<arrow::Int64Builder>({3, 3, 100})},
                                                  {"area_intensity", doubles_({10.0, 20.0, 30.0})},
                                                  {"apex_intensity", doubles_({1.0, 2.0, 3.0})},
                                                  {"score_transition_pep", doubles_({0.01, 0.9, std::nullopt})}});
  }
}

OpenSwathExportFilterConfig filters_()
{
  OpenSwathExportFilterConfig config;
  config.ipf_mode = OpenSwathIPFExportMode::Disable;
  config.peptide = false;
  config.protein = false;
  config.exclude_decoys = false;
  return config;
}

OSWParquetFile::InferenceResults results_(InferenceContext context = InferenceContext::Global)
{
  OSWParquetFile::InferenceResults result;
  result.level = InferenceLevel::Peptide;
  result.context = context;
  LevelContextResultRow row;
  row.entity_id = 0;
  row.context = context;
  row.score = 6.0;
  row.pvalue = 0.01;
  row.qvalue = 0.02;
  row.pep = 0.03;
  if (context != InferenceContext::Global) row.run_id = 11;
  result.rows.push_back(row);
  return result;
}

std::string contents_(const std::string& filename)
{
  std::ifstream input(filename, std::ios::binary);
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}
} // namespace

START_TEST(OSWParquetFile, "$Id$")

START_SECTION((readLevelContextData preserves canonical IDs, contexts, decoys, and null scores))
{
  TempDir temporary;
  const std::string directory = temporary.getPath() + "/workflow with spaces";
  fixture_(directory);
  OSWParquetFile file(directory, library_());
  TEST_EQUAL(file.readRunBasenames().at(11), "sample one")
  TEST_EQUAL(file.readRunBasenames().at(22), "sample two")
  const auto global = file.readLevelContextData(InferenceLevel::Peptide, InferenceContext::Global);
  TEST_EQUAL(global.size(), 2)
  TEST_FALSE(global[0].run_id.has_value())
  TEST_EQUAL(global[0].entity_id, 0)
  TEST_EQUAL(global[0].score, 6.0)
  TEST_FALSE(global[0].decoy)
  TEST_EQUAL(global[1].entity_id, 1)
  TEST_EQUAL(global[1].score, 3.0)
  TEST_TRUE(global[1].decoy)
  for (const auto context : {InferenceContext::ExperimentWide, InferenceContext::RunSpecific})
  {
    const auto rows = file.readLevelContextData(InferenceLevel::Peptide, context);
    TEST_EQUAL(rows.size(), 3)
    TEST_EQUAL(*rows[0].run_id, 11)
    TEST_EQUAL(rows[0].score, 5.0)
    TEST_EQUAL(*rows[1].run_id, 22)
    TEST_EQUAL(rows[1].score, 6.0)
  }
  TEST_EQUAL(file.readLevelContextData(InferenceLevel::Protein, InferenceContext::Global).size(), 2)
  TEST_EQUAL(file.readLevelContextData(InferenceLevel::Gene, InferenceContext::Global).size(), 2)
  TEST_EXCEPTION(Exception::Precondition, file.readLevelContextData(InferenceLevel::Peptidoform, InferenceContext::Global))
  OSWParquetFile moved(std::move(file));
  TEST_EQUAL(moved.readRunBasenames().size(), 2)
}
END_SECTION

START_SECTION((writeLevelContextResults updates all requested contexts and preserves unrelated columns))
{
  TempDir temporary;
  fixture_(temporary.getPath());
  OSWParquetFile file(temporary.getPath(), library_());
  file.writeLevelContextResults({results_(), results_(InferenceContext::RunSpecific)});
  OpenSwathParquetExportConfig config;
  config.filters = filters_();
  const auto table = file.readOpenSwathFeatureScoreTable(config);
  TEST_EQUAL(table.rows.size(), 6)
  Size target_rows = 0;
  for (const auto& row : table.rows)
  {
    if (row.precursor_id == 0)
    {
      ++target_rows;
      TEST_REAL_SIMILAR(*row.score_peptide_global_qvalue, 0.02)
      TEST_EQUAL(row.score_peptide_run_specific_qvalue.has_value(), row.run_id == 11)
    }
    else
    {
      TEST_EQUAL(row.precursor_id, 7)
      TEST_FALSE(row.score_peptide_global_qvalue.has_value())
    }
  }
  TEST_EQUAL(target_rows, 4)
  const auto features = ParquetFile::readTable(temporary.getPath() + "/runs/run_id=11/features.parquet");
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(features, "unrelated_column"), 1), "keep B")
  TEST_TRUE(ParquetFile::getColumn(features, "score_ms2_score")->IsNull(2))
  TEST_TRUE(ParquetFile::getColumn(features, "score_ms2_pep")->IsNull(1))

  // Empty results must clear their requested columns without touching other contexts.
  auto empty = results_();
  empty.rows.clear();
  file.writeLevelContextResults({empty});
  const auto cleared = file.readOpenSwathFeatureScoreTable(config);
  for (const auto& row : cleared.rows)
  {
    TEST_FALSE(row.score_peptide_global_qvalue.has_value())
    if (row.precursor_id == 0 && row.run_id == 11) { TEST_REAL_SIMILAR(*row.score_peptide_run_specific_qvalue, 0.02) }
  }
  auto filtered = filters_();
  filtered.peptide = true;
  TEST_EQUAL(file.readOpenSwathExportRows(filtered).size(), 0)
  const std::string path = temporary.getPath() + "/runs/run_id=11/features.parquet";
  const std::string before = contents_(path);
  auto invalid = results_();
  invalid.rows[0].context = InferenceContext::RunSpecific;
  TEST_EXCEPTION(Exception::InvalidValue, file.writeLevelContextResults({invalid}))
  TEST_EXCEPTION(Exception::InvalidValue, file.writeLevelContextResults({results_(), results_()}))
  invalid = results_();
  invalid.level = InferenceLevel::Peptidoform;
  TEST_EXCEPTION(Exception::Precondition, file.writeLevelContextResults({invalid}))
  TEST_EQUAL(contents_(path), before)
}
END_SECTION

START_SECTION((exports retain canonical precursor and transition IDs and optional scores))
{
  TempDir temporary;
  fixture_(temporary.getPath());
  OSWParquetFile file(temporary.getPath(), library_());
  const auto rows = file.readOpenSwathExportRows(filters_());
  // One unscored feature has no MS2 q-value and is excluded.
  TEST_EQUAL(rows.size(), 5)
  Size target_rows = 0;
  for (const auto& row : rows)
  {
    if (row.precursor_id == 0)
    {
      ++target_rows;
      TEST_EQUAL(row.sequence, "PEPTIDEA")
      TEST_EQUAL(row.protein_name, "P_TARGET")
      TEST_FALSE(row.decoy)
      if (row.feature_id == 110 || row.feature_id == 220)
      {
        TEST_EQUAL(row.aggr_peak_area, "10.0")
        TEST_EQUAL(row.aggr_fragment_annotation, "3_y3_1")
      }
      else
      {
        TEST_EQUAL(row.aggr_peak_area, "") // high transition PEP
      }
    }
    else
    {
      TEST_EQUAL(row.precursor_id, 7)
      TEST_EQUAL(row.sequence, "PEPTIDEC")
      TEST_TRUE(row.decoy)
    }
  }
  TEST_EQUAL(target_rows, 4)
  OpenSwathParquetExportConfig config;
  config.filters = filters_();
  const auto transitions = file.readOpenSwathTransitionScoreTable(config);
  TEST_EQUAL(transitions.rows.size(), 6)
  for (const auto& row : transitions.rows)
  {
    TEST_EQUAL(row.transition_id, row.precursor_id == 0 ? 3 : 100)
  }
  auto invalid = filters_();
  invalid.use_alignment = true;
  TEST_EXCEPTION(Exception::Precondition, file.readOpenSwathExportRows(invalid))
}
END_SECTION

START_SECTION((archives commit explicitly and preserve the original across failure))
{
  TempDir temporary;
  const std::string directory = temporary.getPath() + "/source";
  const std::string archive = temporary.getPath() + "/workflow with spaces.oswpq";
  fixture_(directory);
  ZipArchiveFile::zipDirectory(directory, archive);
  ZipArchiveFile::writeSidecarIndex(archive);
  const std::string original = contents_(archive);
  {
    OSWParquetFile file(archive, library_());
    file.commit();
    TEST_EQUAL(contents_(archive), original)
    file.writeLevelContextResults({results_()});
    TEST_EQUAL(contents_(archive), original)
    // A failing later operation must not prevent explicit preservation of inference.
    auto invalid = filters_();
    invalid.use_alignment = true;
    TEST_EXCEPTION(Exception::Precondition, file.readOpenSwathExportRows(invalid))
    file.commit();
  }
  const std::string committed = contents_(archive);
  TEST_NOT_EQUAL(committed, original)
  OpenSwathParquetExportConfig config;
  config.filters = filters_();
  {
    OSWParquetFile reopened(archive, library_());
    const auto rows = reopened.readOpenSwathFeatureScoreTable(config).rows;
    TEST_REAL_SIMILAR(*rows.front().score_peptide_global_qvalue, 0.02)
    auto empty = results_();
    empty.rows.clear();
    reopened.writeLevelContextResults({empty});
    // Destruction does not persist changes in a disposable archive workspace.
  }
  TEST_EQUAL(contents_(archive), committed)
  {
    OSWParquetFile file(archive, library_());
    file.writeLevelContextResults({results_()});
    TEST_TRUE(File::rename(archive, archive + ".saved", false))
    File::makeDir(archive);
    const std::string sentinel = archive + "/unrelated.txt";
    std::ofstream(sentinel) << "must survive";
    TEST_EXCEPTION(Exception::FileNotWritable, file.commit())
    TEST_EQUAL(contents_(sentinel), "must survive")
    TEST_EQUAL(contents_(archive + ".saved"), committed)
    File::removeDirRecursively(archive);
    TEST_TRUE(File::rename(archive + ".saved", archive, false))
    file.commit(); // a failed commit remains retryable
    OSWParquetFile reopened(archive, library_());
    TEST_EQUAL(reopened.readRunBasenames().size(), 2)
  }
}
END_SECTION

START_SECTION((missing required scores and unknown canonical IDs fail explicitly))
{
  TempDir temporary;
  fixture_(temporary.getPath(), false);
  OSWParquetFile file(temporary.getPath(), library_());
  TEST_EXCEPTION(Exception::Precondition, file.readLevelContextData(InferenceLevel::Peptide, InferenceContext::Global))
  fixture_(temporary.getPath());
  auto library = library_();
  library.compounds[1].id = "42";
  library.transitions[1].peptide_ref = "42";
  OSWParquetFile wrong_library(temporary.getPath(), library);
  TEST_EXCEPTION(Exception::MissingInformation, wrong_library.readLevelContextData(InferenceLevel::Peptide, InferenceContext::Global))
}
END_SECTION

END_TEST
