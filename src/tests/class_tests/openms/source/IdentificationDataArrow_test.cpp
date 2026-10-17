// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdentificationDataArrow.h>
#include <arrow/api.h>
#include <arrow/io/file.h>
#include <filesystem>
#include <parquet/arrow/reader.h>

using namespace OpenMS;
using ID = IdentificationData;
using Arrow = IdentificationDataArrow;

namespace
{
/// Two searches with the same score schema and a scoreless sequence catalog.
struct Fixture
{
  ID data;
  std::vector<ID::MatchId> a, b;
  Fixture()
  {
    ID::ScoreDefinition score;
    score.name = "score";
    ID::ScoreDefinition qvalue;
    qvalue.name = "q-value";
    qvalue.higher_better = false;
    for (const std::string name : {"A", "catalog", "B"})
    {
      auto& run = data.addRun(name);
      const bool catalog = name == "catalog";
      if (catalog)
      {
        ID::RunSettings settings;
        settings.setMetaValue("identification:catalog", "true");
        run.setSettings(settings);
      }
      else
      {
        run.setPrimaryScore(run.addScore(score));
        run.addScore(qvalue);
      }
      const auto source = run.addSource({});
      ID::Database database;
      database.path = "db.fasta";
      const auto db = run.addDatabase(database);
      for (Size i = 0; i < 3; ++i)
      {
        ID::Observation observation;
        observation.data_id = "scan=" + std::to_string(i);
        observation.rt = 10.0 * double(i);
        const auto query = run.addIdentification(source, observation);
        ID::MatchData match;
        match.representation = i == 1 ? "PEPTIDEK" : "PEPTIDER";
        match.charge = 2;
        match.target_decoy = i == 2 ? ID::TargetDecoy::DECOY : ID::TargetDecoy::TARGET;
        match.sequence_evidence = {{db, "P" + std::to_string(i), 3, 10, 'K', 'A'}};
        match.setMetaValue("rank", int(i));
        if (catalog) run.addMatch(query, match);
        else
        {
          const auto id = run.addMatch(query, match, {double(i), i == 0 ? std::nullopt : std::optional<double>(0.01 * double(i))});
          (name == "A" ? a : b).push_back(id);
        }
      }
    }
  }
  const ID::Run& run(const std::string& name) const
  { return data.getRun(name); }
};

std::shared_ptr<arrow::Array> finish(arrow::ArrayBuilder& builder)
{
  std::shared_ptr<arrow::Array> array;
  if (! builder.Finish(&array).ok()) throw std::runtime_error("cannot build test column");
  return array;
}
std::shared_ptr<arrow::Array> strings(const std::vector<std::string>& values)
{
  arrow::StringBuilder builder;
  for (const auto& value : values)
    (void) builder.Append(value);
  return finish(builder);
}
template<class Builder, class T>
std::shared_ptr<arrow::Array> numbers(const std::vector<std::optional<T>>& values)
{
  Builder builder;
  for (const auto& value : values)
    (void) (value ? builder.Append(*value) : builder.AppendNull());
  return finish(builder);
}
std::shared_ptr<arrow::Array> doubles(const std::vector<std::optional<double>>& values)
{ return numbers<arrow::DoubleBuilder>(values); }

/// A patch table: run_uuid and match_id keys, then the given columns.
std::shared_ptr<arrow::Table> patch(const std::vector<std::string>& runs, const std::vector<ID::MatchId>& matches,
                                    const std::vector<std::pair<std::string, std::shared_ptr<arrow::Array>>>& columns = {})
{
  std::vector<std::optional<UInt64>> ids;
  for (const auto& match : matches)
    ids.push_back(match.value);
  std::vector<std::shared_ptr<arrow::Field>> fields {arrow::field("run_uuid", arrow::utf8()), arrow::field("match_id", arrow::uint64())};
  std::vector<std::shared_ptr<arrow::Array>> arrays {strings(runs), numbers<arrow::UInt64Builder>(ids)};
  for (const auto& [name, array] : columns)
  {
    fields.push_back(arrow::field(name, array->type()));
    arrays.push_back(array);
  }
  return arrow::Table::Make(arrow::schema(fields), arrays);
}

/// The run UUIDs of a dictionary-encoded key column.
std::vector<std::string> keys(const arrow::ChunkedArray& column)
{
  std::vector<std::string> values;
  for (const auto& chunk : column.chunks())
  {
    const auto& encoded = static_cast<const arrow::DictionaryArray&>(*chunk);
    for (int64_t row = 0; row < encoded.length(); ++row)
      values.push_back(static_cast<const arrow::StringArray&>(*encoded.dictionary()).GetString(encoded.GetValueIndex(row)));
  }
  return values;
}
Arrow::PatchOptions adding(const ID::ScoreDefinition& definition)
{
  Arrow::PatchOptions options;
  options.add_scores = {definition};
  return options;
}
Arrow::PatchOptions expecting(const std::string& run, UInt64 revision)
{
  Arrow::PatchOptions options;
  options.expected_revisions[run] = revision;
  return options;
}
std::optional<double> score(const ID::Run& run, ID::MatchId match, UInt32 index)
{ return run.getScore(match, run.getScoreId(index)); }
} // namespace

START_TEST(IdentificationDataArrow, "$Id$")

START_SECTION((static std::shared_ptr<arrow::Table> matchTable(const IdentificationData&)))
{
  Fixture fixture;
  const auto table = Arrow::matchTable(fixture.data);
  TEST_EQUAL(table->num_rows(), 9)
  // The columns of matches.parquet in a bundle, with the same values (apart from the metadata column, whose
  // descriptor indices refer to the descriptors of the table instead of the bundle manifest).
  std::string path;
  NEW_TMP_FILE(path)
  IdentificationDataFile::store(path, fixture.data);
  std::shared_ptr<arrow::Table> stored;
  {
    auto source = arrow::io::ReadableFile::Open(path + "/matches.parquet").ValueOrDie();
    stored = parquet::arrow::OpenFile(source, arrow::default_memory_pool()).ValueOrDie()->ReadTable().ValueOrDie();
  }
  std::filesystem::remove_all(path);
  TEST_EQUAL(table->num_columns(), stored->num_columns())
  for (int i = 0; i < table->num_columns(); ++i)
  {
    const auto& name = table->schema()->field(i)->name();
    TEST_EQUAL(name, stored->schema()->field(i)->name())
    TEST_TRUE(table->schema()->field(i)->type()->Equals(*stored->schema()->field(i)->type()))
    if (name == "run_uuid") TEST_TRUE(keys(*table->column(i)) == keys(*stored->column(i)))
    else if (name != "metadata")
      TEST_TRUE(table->column(i)->Equals(*stored->column(i)))
  }
  const std::vector<std::string> uuids {fixture.run("A").getUuid(), fixture.run("A").getUuid(), fixture.run("A").getUuid(),
                                        fixture.run("catalog").getUuid(), fixture.run("catalog").getUuid(), fixture.run("catalog").getUuid(),
                                        fixture.run("B").getUuid(), fixture.run("B").getUuid(), fixture.run("B").getUuid()};
  TEST_TRUE(keys(*table->GetColumnByName("run_uuid")) == uuids)
  const auto qvalues = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("score_q_value")->chunk(0));
  TEST_TRUE(qvalues->IsNull(0))
  TEST_REAL_SIMILAR(qvalues->Value(2), 0.02)
  TEST_TRUE(qvalues->IsNull(4)) // the catalog has no scores

  // The schema metadata describes the score columns, the metadata descriptors and the revisions of the runs.
  // (JSON, as in the manifest of a bundle.)
  const auto& metadata = *table->schema()->metadata();
  const auto contains = [&](const std::string& key, const std::string& text) { return metadata.Get(key).ValueOrDie().find(text) != std::string::npos; };
  TEST_TRUE(contains("openms:score_definitions", "{\"name\":\"score\""))
  TEST_TRUE(contains("openms:score_definitions", "{\"name\":\"q-value\""))
  TEST_EQUAL(metadata.Get("openms:primary_score").ValueOrDie(), "0")
  TEST_TRUE(contains("openms:metadata_descriptors", "[{\"name\":\"rank\""))
  for (const auto& run : fixture.data.getRuns())
    TEST_TRUE(contains("openms:revisions", "\"" + run.getUuid() + "\":" + std::to_string(run.getRevision())))

  // An empty dataset has the fixed columns and no rows.
  const auto empty = Arrow::matchTable(ID {});
  TEST_EQUAL(empty->num_rows(), 0)
  TEST_EQUAL(empty->num_columns(), 15)
  TEST_EQUAL(empty->schema()->metadata()->Get("openms:primary_score").ValueOrDie(), "null")
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> matchTable(const IdentificationData&, const IdentificationDataFile::Projection&)))
{
  Fixture fixture;
  IdentificationDataFile::Projection projection;
  projection.molecule = projection.evidence = projection.annotations = projection.metadata = projection.all_scores = false;
  projection.score_ids = {1};
  const auto table = Arrow::matchTable(fixture.data, projection);
  TEST_TRUE(table->schema()->field_names() == (std::vector<std::string> {"match_id", "query_id", "score_q_value", "run_uuid"}))
  TEST_EQUAL(table->num_rows(), 9)
  TEST_TRUE(table->schema()->metadata()->Contains("openms:revisions"))
  projection.score_ids = {2};
  TEST_EXCEPTION(Exception::InvalidValue, Arrow::matchTable(fixture.data, projection))
}
END_SECTION

START_SECTION((static PatchResult applyPatch(IdentificationData&, const arrow::Table&)))
{
  Fixture fixture;
  const auto a = fixture.run("A").getUuid(), b = fixture.run("B").getUuid();
  auto& run_b = fixture.data.getRun("B");
  const auto view = run_b.bindScore(run_b.getScoreId(1));
  const auto revision = run_b.getRevision();
  const auto result = Arrow::applyPatch(fixture.data, *patch({a, b, b}, {fixture.a[0], fixture.b[1], fixture.b[2]},
                                                              {{"score_q_value", doubles({0.5, std::nullopt, 0.25})},
                                                               {"target_decoy", numbers<arrow::Int32Builder, int>({2, 0, 1})}}));
  TEST_EQUAL(result.rows, 3)
  TEST_TRUE(result.columns == (std::vector<std::string> {"score_q_value", "target_decoy"}))
  TEST_TRUE(result.added.empty())
  const auto& run_a = fixture.run("A");
  TEST_REAL_SIMILAR(*score(run_a, fixture.a[0], 1), 0.5)
  TEST_REAL_SIMILAR(*score(run_a, fixture.a[0], 0), 0.0) // the primary score stays
  TEST_TRUE(run_a.getMatch(fixture.a[0]).target_decoy == ID::TargetDecoy::DECOY)
  TEST_TRUE(! score(run_b, fixture.b[1], 1)) // null removes a supplementary value
  TEST_REAL_SIMILAR(*score(run_b, fixture.b[2], 1), 0.25)
  TEST_TRUE(run_b.getMatch(fixture.b[2]).target_decoy == ID::TargetDecoy::TARGET)
  TEST_TRUE(run_b.getMatch(fixture.b[1]).target_decoy == ID::TargetDecoy::UNKNOWN)
  TEST_TRUE(! score(run_b, fixture.b[0], 1)) // untouched rows keep their values
  TEST_REAL_SIMILAR(*score(run_a, fixture.a[2], 1), 0.02)
  TEST_TRUE(run_b.getRevision() > revision)
  // The patched run is a new state: score views bound before are rejected, new views read the patched values.
  TEST_EXCEPTION(Exception::InvalidValue, view(run_b.getMatch(fixture.b[2])))
  TEST_REAL_SIMILAR(*run_b.bindScore(run_b.getScoreId(1))(run_b.getMatch(fixture.b[2])), 0.25)
  fixture.data.validate();
  // The catalog was not touched.
  TEST_EQUAL(fixture.run("catalog").getScoreDefinitions().size(), 0)

  // A table of the dataset is a patch of its own score columns that changes no value; dictionary keys are accepted.
  const auto before = fixture.data;
  const auto table = Arrow::matchTable(fixture.data);
  const auto own = table->SelectColumns({table->num_columns() - 1, 0, table->num_columns() - 3, table->num_columns() - 2}).ValueOrDie();
  TEST_TRUE(own->schema()->field_names() == (std::vector<std::string> {"run_uuid", "match_id", "score_score", "score_q_value"}))
  // The catalog rows have no scores to patch.
  TEST_EXCEPTION(Exception::InvalidParameter, Arrow::applyPatch(fixture.data, *own))
  const auto searches = arrow::ConcatenateTables({own->Slice(0, 3), own->Slice(6, 3)}).ValueOrDie();
  TEST_EQUAL(Arrow::applyPatch(fixture.data, *searches).rows, 6)
  TEST_TRUE(fixture.data == before)
}
END_SECTION

START_SECTION((static PatchResult applyPatch(IdentificationData&, const arrow::Table&, const PatchOptions&)))
{
  Fixture fixture;
  const auto a = fixture.run("A").getUuid(), b = fixture.run("B").getUuid();
  Arrow::PatchOptions options;
  ID::ScoreDefinition rescored;
  rescored.name = "PEP";
  rescored.software = "rescorer";
  rescored.higher_better = false;
  options.add_scores = {rescored, rescored};
  options.expected_revisions[b] = fixture.run("B").getRevision();
  const auto result = Arrow::applyPatch(fixture.data, *patch({b, b}, {fixture.b[0], fixture.b[2]}, {{"score_pep", doubles({0.1, std::nullopt})}}), options);
  TEST_EQUAL(result.rows, 2)
  TEST_TRUE(result.added == std::vector<std::string> {"score_pep"})
  TEST_TRUE(result.columns == std::vector<std::string> {"score_pep"})
  // The definition is added to every run of the score schema (once), so the dataset keeps one schema.
  TEST_EQUAL(fixture.data.getScoreDefinitions().size(), 3)
  TEST_TRUE(fixture.data.getScoreDefinitions()[2] == rescored)
  TEST_EQUAL(fixture.run("A").getScoreDefinitions().size(), 3)
  TEST_EQUAL(fixture.run("catalog").getScoreDefinitions().size(), 0)
  TEST_TRUE(! score(fixture.run("A"), fixture.a[0], 2))
  TEST_REAL_SIMILAR(*score(fixture.run("B"), fixture.b[0], 2), 0.1)
  TEST_TRUE(! score(fixture.run("B"), fixture.b[1], 2))
  fixture.data.validate();
  // An existing definition is reused.
  TEST_TRUE(Arrow::applyPatch(fixture.data, *patch({a}, {fixture.a[1]}, {{"score_pep", doubles({0.2})}}), adding(rescored)).added.empty())
  TEST_REAL_SIMILAR(*score(fixture.run("A"), fixture.a[1], 2), 0.2)

  // Invalid patches change nothing.
  const auto before = fixture.data;
  const auto rejected = [&](const std::shared_ptr<arrow::Table>& table, const Arrow::PatchOptions& with = {}) {
    try
    {
      Arrow::applyPatch(fixture.data, *table, with);
      return false;
    }
    catch (const Exception::InvalidParameter&)
    {
      return fixture.data == before;
    }
  };
  const auto valid = doubles({0.3, 0.4});
  // Every row is checked before anything changes: the first row is valid, the second not.
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], ID::MatchId {999}}, {{"score_q_value", valid}})))
  TEST_TRUE(rejected(patch({a, "unknown"}, {fixture.a[0], fixture.a[1]}, {{"score_q_value", valid}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[0]}, {{"score_q_value", valid}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"score_unknown", valid}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"charge", numbers<arrow::Int32Builder, int>({1, 2})}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"score_score", doubles({1.0, std::nullopt})}}))) // the primary score
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"score_q_value", doubles({0.1, std::numeric_limits<double>::infinity()})}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"score_q_value", strings({"0.1", "0.2"})}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"target_decoy", numbers<arrow::Int32Builder, int>({1, 4})}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"target_decoy", numbers<arrow::Int32Builder, int>({1, -1})}})))
  TEST_TRUE(rejected(patch({a, a}, {fixture.a[0], fixture.a[1]}, {{"score_q_value", valid}, {"score_q_value", valid}})))
  TEST_TRUE(rejected(arrow::Table::Make(arrow::schema({arrow::field("match_id", arrow::uint64())}), {numbers<arrow::UInt64Builder, UInt64>({1})})))
  // A score definition added by a rejected patch is not added either.
  ID::ScoreDefinition other;
  other.name = "other";
  TEST_TRUE(rejected(patch({a}, {ID::MatchId {999}}, {{"score_other", doubles({0.5})}}), adding(other)))
  // A run that changed since the patch was computed.
  TEST_TRUE(rejected(patch({a}, {fixture.a[0]}, {{"score_q_value", doubles({0.5})}}),
                     expecting(a, fixture.run("A").getRevision() + 1)))
  TEST_TRUE(rejected(patch({a}, {fixture.a[0]}), expecting("unknown", 0)))
  // A patch of keys only checks them, and edits nothing: score views stay valid.
  const auto& run_a = fixture.run("A");
  const auto view = run_a.bindScore(run_a.getScoreId(1));
  TEST_EQUAL(Arrow::applyPatch(fixture.data, *patch({a}, {fixture.a[0]})).rows, 1)
  TEST_TRUE(fixture.data == before)
  TEST_REAL_SIMILAR(*view(run_a.getMatch(fixture.a[2])), 0.02)
}
END_SECTION

START_SECTION((static PatchResult applyPatch(IdentificationData&, arrow::RecordBatchReader&, const PatchOptions&)))
{
  Fixture fixture;
  const auto a = fixture.run("A").getUuid();
  const auto table = patch({a, a}, {fixture.a[0], fixture.a[2]}, {{"score_q_value", doubles({0.125, 0.5})}});
  arrow::TableBatchReader reader(*table);
  reader.set_chunksize(1);
  TEST_EQUAL(Arrow::applyPatch(fixture.data, reader, {}).rows, 2)
  TEST_REAL_SIMILAR(*score(fixture.run("A"), fixture.a[0], 1), 0.125)
  TEST_REAL_SIMILAR(*score(fixture.run("A"), fixture.a[2], 1), 0.5)
  // Keys as other libraries write them: large strings and string views (polars), signed integers (pandas).
  for (const auto& type : {arrow::large_utf8(), arrow::utf8_view()})
  {
    arrow::StringViewBuilder views;
    arrow::LargeStringBuilder large;
    arrow::ArrayBuilder& builder = type->id() == arrow::Type::STRING_VIEW ? static_cast<arrow::ArrayBuilder&>(views) : large;
    TEST_TRUE((type->id() == arrow::Type::STRING_VIEW ? views.Append(a) : large.Append(a)).ok())
    const auto keys = arrow::Table::Make(arrow::schema({arrow::field("run_uuid", type), arrow::field("match_id", arrow::int64()),
                                                        arrow::field("score_q_value", arrow::float32())}),
                                         {finish(builder), numbers<arrow::Int64Builder, Int64>({Int64(fixture.a[1].value)}), numbers<arrow::FloatBuilder, float>({0.25f})});
    TEST_EQUAL(Arrow::applyPatch(fixture.data, *keys).rows, 1)
    TEST_REAL_SIMILAR(*score(fixture.run("A"), fixture.a[1], 1), 0.25)
  }
}
END_SECTION

END_TEST
