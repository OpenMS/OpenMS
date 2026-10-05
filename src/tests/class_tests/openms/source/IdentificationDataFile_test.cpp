// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <bit>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>

using namespace OpenMS;
using ID = IdentificationData;
using Native = IdentificationDataFile;
namespace fs = std::filesystem;

namespace
{
struct Fixture
{
  ID data;
  std::string uuid;
  ID::QueryId query, empty;
  ID::MatchId first, selected, last;
  Fixture()
  {
    auto& run = data.addRun("search-one");
    uuid = run.getUuid();
    ID::ScoreDefinition raw;
    raw.name = "engine score";
    raw.software = "engine";
    raw.aggregation = "PSM";
    raw.parameters.setMetaValue("threshold", -0.0);
    auto score = run.addScore(raw);
    ID::ScoreDefinition qvalue;
    qvalue.name = "q-value";
    qvalue.higher_better = false;
    run.addScore(qvalue);
    ID::SourceFile source;
    source.identifier = "raw-source";
    source.path = "/exact path/ä/sample.raw";
    source.primary_files = {"first.raw", "second.raw", "first.raw"};
    source.setMetaValue("run", IntList {1, 2});
    auto source_id = run.addSource(source);
    ID::Observation observation;
    observation.data_id = "controllerType=0 scan=1";
    observation.rt = 1.25;
    observation.mz = 345.67;
    observation.setMetaValue("zero", 0);
    observation.setMetaValue("empty", DataValue {});
    query = run.addIdentification(source_id, observation);
    observation.data_id = "empty-observation";
    empty = run.addIdentification(source_id, observation);
    ID::MatchData match;
    match.representation = "PEP(Custom definition deliberately not registered)TIDE";
    match.charge = 2;
    match.calculated_mz = 345.66;
    match.target_decoy = ID::TargetDecoy::BOTH;
    match.name = "candidate";
    match.formula = "C2H4";
    match.identifiers = {{"db", "molecule"}, {"other", "molecule"}};
    match.adduct.emplace("original adduct name", EmpiricalFormula("H2"), 2, 1);
    match.parent_evidence = {{{"db", "P1"}, std::nullopt, 7, "-", "K"}, {{"other", "P1"}, 2, std::nullopt, "R", "-"}};
    PeptideHit::PeakAnnotation annotation;
    annotation.annotation = "y7-H2O";
    annotation.charge = 1;
    annotation.mz = 700.5;
    annotation.intensity = 0.25;
    match.peak_annotations.push_back(annotation);
    DataValue units(DoubleList {0.0, -0.0, std::bit_cast<double>(UInt64 {0x7ff8000000000042}), std::numeric_limits<double>::infinity()});
    units.setUnitType(DataValue::MS_ONTOLOGY);
    units.setUnit(1000001);
    match.setMetaValue("double-list", units);
    match.setMetaValue("empty-string", "");
    match.setMetaValue("empty-list", StringList {});
    match.setMetaValue("empty", DataValue {});
    match.setMetaValue("integer", static_cast<long long>(9223372036854775807LL));
    match.setMetaValue("integers", IntList {-2, 0, 3});
    match.setMetaValue("strings", StringList {"", "α", "comma,value"});
    match.setMetaValue("nan", std::bit_cast<double>(UInt64 {0x7ff8000000000031}));
    first = run.addMatch(query, match, {10.0, std::nullopt});
    selected = run.addMatch(query, match, {20.0, 0.0});
    match.adduct.reset();
    match.formula.reset();
    last = run.addMatch(query, match, {30.0, 0.01});
    run.setSelectedMatch(query, selected);
    run.setPrimaryScore(score);
    ID::ParentRecord parent;
    parent.identity = {"db", "P1"};
    parent.sequence = "PEPTIDE";
    parent.description = "parent description";
    parent.target_decoy = ID::TargetDecoy::TARGET;
    parent.setMetaValue("parent-int", IntList {});
    run.setParents(std::vector<ID::ParentRecord> {parent});
    ID::InferenceResult inference;
    inference.identifier = "pooled";
    ProteinHit protein(0.99, 1, "P1", "PEPTIDE");
    protein.setMetaValue("integer", 7);
    inference.proteins.setHits({protein});
    inference.qualified_accessions["P1"] = {"db", "P1"};
    ID::InferenceInput input;
    input.run_identifier = run.getIdentifier();
    input.run_uuid = uuid;
    input.score = raw;
    input.matches = {selected, first, selected};
    input.selection = "all candidates";
    inference.inputs.push_back(input);
    ID::InferenceInput absent;
    absent.run_identifier = "not-in-export";
    absent.run_uuid = "11111111-1111-4111-8111-111111111111";
    absent.matches = {ID::MatchId {selected.value}};
    inference.inputs.push_back(absent);
    inference.assignments.push_back({run.getIdentifier(), uuid, selected, UInt64 {0}, {{"db", "P1"}}});
    inference.assignments.push_back({run.getIdentifier(), uuid, first, UInt64 {0}, {}});
    data.addInferenceResult(inference);
    auto& other = data.addRun("empty-compounds", ID::MoleculeKind::COMPOUND);
    other.setParents(std::vector<ID::ParentRecord> {});
  }
};
Native::Options tiny()
{
  Native::Options options;
  options.batch_rows = 1;
  options.row_group_rows = 1;
  options.batch_bytes = 64;
  options.row_group_bytes = 128;
  return options;
}
void replaceText(const fs::path& path, const std::string& from, const std::string& to)
{
  std::ifstream in(path, std::ios::binary);
  std::string content((std::istreambuf_iterator<char>(in)), {});
  in.close();
  auto offset = content.find(from);
  if (offset == std::string::npos) throw std::runtime_error("test replacement not found");
  content.replace(offset, from.size(), to);
  std::ofstream out(path, std::ios::binary);
  out << content;
}
} // namespace

START_TEST(IdentificationDataFile, "$Id$")

START_SECTION((shared row groups preserve run boundaries, metadata and supplementary score layouts))
{
  ID data;
  ID::ScoreDefinition primary;
  primary.name = "common";
  for (Size index = 0; index < 4; ++index)
  {
    auto& run = data.addRun("run" + std::to_string(index));
    auto score = run.addScore(primary);
    if (index == 3)
    {
      auto additional = primary;
      additional.name = "extra";
      run.addScore(additional);
    }
    run.setPrimaryScore(score);
    auto source = run.addSource({});
    for (Size i = 0; i < 3; ++i)
    {
      ID::Observation observation;
      observation.data_id = "scan=" + std::to_string(i);
      auto query = run.addIdentification(source, observation);
      ID::MatchData match;
      match.representation = "PEPTIDE";
      match.setMetaValue("run" + std::to_string(index), static_cast<int>(index));
      std::vector<std::optional<double>> values {double(index * 10 + i)};
      if (index == 3) values.push_back(std::nullopt);
      run.addMatch(query, match, values);
    }
  }
  Native::Options options;
  options.row_group_rows = 5; // boundaries cross runs, and runs cross row groups
  options.batch_rows = 2;
  std::string path, filtered;
  NEW_TMP_FILE(path)
  NEW_TMP_FILE(filtered)
  Native::store(path, data, options);
  TEST_TRUE(fs::exists(fs::path(path) / "queries-0.parquet"))
  TEST_TRUE(fs::exists(fs::path(path) / "matches-0.parquet"))
  TEST_TRUE(fs::exists(fs::path(path) / "matches-1.parquet"))
  TEST_FALSE(fs::exists(fs::path(path) / "runs"))
  ID loaded;
  Native::load(path, loaded, options);
  TEST_EQUAL(loaded.getRuns().size(), 4)
  for (Size index = 0; index < 4; ++index)
  {
    const auto& run = loaded.getRuns()[index];
    TEST_EQUAL(run.getNumberOfMatches(), 3)
    const auto& match = run.getSourceBlocks()[0].identifications[0].getMatches()[0];
    TEST_EQUAL(match.getMetaValue("run" + std::to_string(index)), static_cast<int>(index))
    TEST_REAL_SIMILAR(*run.getScore(match.getId(), *run.getPrimaryScore()), index * 10)
  }
  auto selected = Native::loadRun(path, "run1", options);
  TEST_EQUAL(selected.getNumberOfMatches(), 3)
  Native::filter(path, filtered, [](const auto&, const Native::MatchRecord& match) { return *match.scores[0] >= 20; },
                 ID::InferencePolicy::DISCARD, false, options);
  Native::load(filtered, loaded, options);
  TEST_EQUAL(loaded.getRuns()[0].getNumberOfMatches(), 0)
  TEST_EQUAL(loaded.getRuns()[2].getNumberOfMatches(), 3)
  TEST_EQUAL(loaded.getRuns()[3].getNumberOfMatches(), 3)
  replaceText(fs::path(path) / "manifest.json", "\"partition\": 0", "\"partition\": 99");
  TEST_EXCEPTION(Exception::InvalidValue, Native::scan(path, {}, {}, {}))
  fs::remove_all(path);
  fs::remove_all(filtered);
}
END_SECTION

START_SECTION((common primary score is checked before writing and scanning))
{
  ID data;
  ID::ScoreDefinition score;
  score.name = "common";
  auto& a = data.addRun("A");
  a.setPrimaryScore(a.addScore(score));
  auto& b = data.addRun("B");
  b.setPrimaryScore(b.addScore(score));
  std::string directory;
  NEW_TMP_FILE(directory)
  Native::store(directory, data);
  // Individually valid but incompatible definitions fail before streaming any table.
  replaceText(fs::path(directory) / "manifest.json", "\"name\": \"common\"", "\"name\": \"different\"");
  TEST_EXCEPTION(Exception::InvalidValue, Native::scan(directory, {}, {}, {}))
  fs::remove_all(directory);
  score.name = "other";
  score.higher_better = false;
  b.setPrimaryScore(b.addScore(score));
  TEST_EXCEPTION(Exception::InvalidValue, Native::store(directory, data))
  TEST_FALSE(fs::exists(directory))
}
END_SECTION

START_SECTION((static void store(const std::string&, const IdentificationData&, const Options&)))
{
  Fixture fixture;
  std::string directory;
  NEW_TMP_FILE(directory);
  fs::remove_all(directory);
  Native::store(directory, fixture.data, tiny());
  TEST_TRUE(Native::isNativeFile(directory))
  TEST_TRUE(fs::exists(fs::path(directory) / "queries-0.parquet"))
  TEST_TRUE(fs::exists(fs::path(directory) / "matches-0.parquet"))
  TEST_TRUE(fs::exists(fs::path(directory) / "input_members-0.parquet"))
  TEST_FALSE(fs::exists(fs::path(directory) / "inference.json"))
  ID loaded;
  Native::load(directory, loaded, tiny());
  const auto& run = loaded.getRun("search-one");
  TEST_EQUAL(run.getUuid(), fixture.uuid)
  TEST_EQUAL(run.getNextMatchId(), fixture.data.getRun("search-one").getNextMatchId())
  TEST_EQUAL(run.getNumberOfIdentifications(), 2)
  TEST_EQUAL(run.getNumberOfMatches(), 3)
  TEST_EQUAL(run.getIdentification(fixture.empty).getMatches().size(), 0)
  TEST_EQUAL(run.getIdentification(fixture.query).getSelectedMatch()->value, fixture.selected.value)
  TEST_EQUAL(run.getSourceBlocks()[0].source.path, "/exact path/ä/sample.raw")
  TEST_EQUAL(run.getSourceBlocks()[0].source.primary_files.size(), 3)
  TEST_TRUE(run.getScoreDefinitions() == fixture.data.getRun("search-one").getScoreDefinitions())
  const auto& match = run.getMatch(fixture.first);
  TEST_EQUAL(match.representation, fixture.data.getRun("search-one").getMatch(fixture.first).representation)
  TEST_TRUE(match.adduct == fixture.data.getRun("search-one").getMatch(fixture.first).adduct)
  TEST_FALSE(run.getMatch(fixture.last).adduct.has_value())
  TEST_FALSE(run.getMatch(fixture.last).formula.has_value())
  TEST_TRUE(match.parent_evidence == fixture.data.getRun("search-one").getMatch(fixture.first).parent_evidence)
  TEST_EQUAL(match.peak_annotations[0].annotation, "y7-H2O")
  TEST_FALSE(match.getScores()[1].has_value())
  TEST_EQUAL(*run.getMatch(fixture.selected).getScores()[1], 0.0)
  TEST_TRUE(match.metaValueExists("empty"))
  TEST_EQUAL(match.getMetaValue("empty").valueType(), DataValue::EMPTY_VALUE)
  TEST_EQUAL(match.getMetaValue("empty-list").valueType(), DataValue::STRING_LIST)
  TEST_EQUAL(match.getMetaValue("empty-list").toStringList().size(), 0)
  TEST_EQUAL(static_cast<long long>(match.getMetaValue("integer")), 9223372036854775807LL)
  TEST_EQUAL(std::bit_cast<UInt64>(static_cast<double>(match.getMetaValue("nan"))), UInt64 {0x7ff8000000000031})
  const auto values = match.getMetaValue("double-list").toDoubleList();
  TEST_EQUAL(std::bit_cast<UInt64>(values[0]), UInt64 {0})
  TEST_EQUAL(std::bit_cast<UInt64>(values[1]), UInt64 {0x8000000000000000})
  TEST_EQUAL(std::bit_cast<UInt64>(values[2]), UInt64 {0x7ff8000000000042})
  TEST_EQUAL(match.getMetaValue("double-list").getUnit(), 1000001)
  TEST_EQUAL(match.getMetaValue("double-list").getUnitType(), DataValue::MS_ONTOLOGY)
  TEST_EQUAL(run.getParents()->at(0).description, "parent description")
  TEST_TRUE(loaded.getRun("empty-compounds").getParents().has_value())
  TEST_EQUAL(loaded.getRun("empty-compounds").getParents()->size(), 0)
  TEST_TRUE(loaded.getInferenceResults()[0].inputs[0].matches == fixture.data.getInferenceResults()[0].inputs[0].matches)
  TEST_EQUAL(loaded.getInferenceResults()[0].assignments[1].parents.size(), 0)
  TEST_EQUAL(loaded.getInferenceResults()[0].inputs[1].run_uuid, "11111111-1111-4111-8111-111111111111")
  TEST_EXCEPTION(Exception::InvalidValue, Native::store(directory, fixture.data))
  fs::remove_all(directory);
}
END_SECTION

START_SECTION((static ScanStatistics scan(const std::string&, const ScanOptions&, const QueryCallback&, const MatchCallback&)))
{
  Fixture fixture;
  std::string directory;
  NEW_TMP_FILE(directory);
  fs::remove_all(directory);
  Native::store(directory, fixture.data, tiny());
  const auto headers = Native::inspect(directory);
  TEST_EQUAL(headers.size(), 2)
  TEST_EQUAL(headers[0].query_count, 2)
  TEST_EQUAL(headers[0].match_count, 3)
  auto single = Native::loadRun(directory, fixture.uuid, tiny());
  TEST_EQUAL(single.getNumberOfMatches(), 3)
  Native::ScanOptions scan;
  scan.buffering = tiny();
  scan.runs = {fixture.uuid};
  scan.validate_unique_ids = true;
  scan.projection.molecule = false;
  scan.projection.evidence = false;
  scan.projection.annotations = false;
  scan.projection.metadata = false;
  scan.projection.all_scores = false;
  scan.projection.score_ids = {1};
  Size query_count = 0, match_count = 0;
  double sum = 0;
  auto statistics = Native::scan(
    directory, scan,
    [&](const std::string& uuid, const std::vector<Native::QueryRecord>& records) {
      TEST_EQUAL(uuid, fixture.uuid) TEST_EQUAL(records.size(), 1) query_count += records.size();
    },
    [&](const std::string& uuid, const std::vector<Native::MatchRecord>& records) {
      TEST_EQUAL(uuid, fixture.uuid) TEST_EQUAL(records.size(), 1) for (const auto& m : records)
      {
        ++match_count;
        TEST_TRUE(m.data.representation.empty())
        TEST_TRUE(m.data.parent_evidence.empty()) TEST_FALSE(m.scores[0].has_value()) if (m.scores[1]) sum += *m.scores[1];
      }
    });
  TEST_EQUAL(query_count, 2)
  TEST_EQUAL(match_count, 3)
  TEST_REAL_SIMILAR(sum, 0.01)
  TEST_EQUAL(statistics.queries, 2)
  TEST_EQUAL(statistics.matches, 3)
  TEST_TRUE(statistics.descriptor_bytes > 0)
  // Run selection must not touch another run's missing tables.
  fs::remove(fs::path(directory) / "matches-1.parquet");
  single = Native::loadRun(directory, "search-one", tiny());
  TEST_EQUAL(single.getNumberOfMatches(), 3)
  ID old;
  old.addRun("destination");
  TEST_EXCEPTION(Exception::InvalidValue, Native::load(directory, old, tiny()))
  TEST_EQUAL(old.getRuns().front().getIdentifier(), "destination")
  fs::remove_all(directory);
}
END_SECTION

START_SECTION(
  (static void filter(const std::string&, const std::string&, const MatchPredicate&, IdentificationData::InferencePolicy, bool, const Options&)))
{
  Fixture fixture;
  std::string input, preserved, discarded, selected;
  NEW_TMP_FILE(input);
  fs::remove_all(input);
  NEW_TMP_FILE(preserved);
  fs::remove_all(preserved);
  NEW_TMP_FILE(discarded);
  fs::remove_all(discarded);
  NEW_TMP_FILE(selected);
  fs::remove_all(selected);
  Native::store(input, fixture.data, tiny());
  Native::MatchPredicate keep_first = [&](const std::string&, const Native::MatchRecord& m) { return m.match_id == fixture.first.value; };
  Native::filter(input, preserved, keep_first, ID::InferencePolicy::PRESERVE, false, tiny());
  ID reduced;
  Native::load(preserved, reduced, tiny());
  auto& run = reduced.getRun("search-one");
  TEST_EQUAL(run.getUuid(), fixture.uuid)
  TEST_EQUAL(run.getNumberOfIdentifications(), 1)
  TEST_EQUAL(run.getNumberOfMatches(), 1)
  TEST_EQUAL(run.getMatch(fixture.first).getId().value, fixture.first.value)
  TEST_FALSE(run.getIdentification(fixture.query).getSelectedMatch().has_value())
  TEST_TRUE(run.findMatch(fixture.selected) == nullptr)
  TEST_TRUE(reduced.getInferenceResults()[0].inputs[0].matches == fixture.data.getInferenceResults()[0].inputs[0].matches)
  TEST_EQUAL(reduced.getInferenceResults()[0].assignments[0].match.value, fixture.selected.value)
  auto appended = run.addMatch(fixture.query, run.getMatch(fixture.first).getData(), {40.0, 0.02});
  TEST_TRUE(appended.value > fixture.last.value)
  Native::filter(input, discarded, keep_first, ID::InferencePolicy::DISCARD, true, tiny());
  Native::load(discarded, reduced, tiny());
  TEST_EQUAL(reduced.getInferenceResults().size(), 0)
  TEST_EQUAL(reduced.getRun("search-one").getNumberOfIdentifications(), 2)
  TEST_FALSE(fs::exists(fs::path(discarded) / "inference"))
  Native::MatchPredicate keep_selected = [&](const std::string&, const Native::MatchRecord& m) { return m.match_id == fixture.selected.value; };
  Native::filter(input, selected, keep_selected, ID::InferencePolicy::PRESERVE, true, tiny());
  Native::load(selected, reduced, tiny());
  TEST_EQUAL(reduced.getRun("search-one").getIdentification(fixture.query).getSelectedMatch()->value, fixture.selected.value)
  for (const auto& path : {input, preserved, discarded, selected})
    fs::remove_all(path);
}
END_SECTION

START_SECTION((write failures and callback exceptions never publish a dataset))
{
  Fixture fixture;
  std::string input, callback_output, text_output, size_output;
  NEW_TMP_FILE(input);
  fs::remove_all(input);
  NEW_TMP_FILE(callback_output);
  fs::remove_all(callback_output);
  NEW_TMP_FILE(text_output);
  fs::remove_all(text_output);
  NEW_TMP_FILE(size_output);
  fs::remove_all(size_output);
  Native::store(input, fixture.data, tiny());
  Native::MatchPredicate throwing = [](const std::string&, const Native::MatchRecord&) -> bool { throw std::runtime_error("callback"); };
  TEST_EXCEPTION(std::runtime_error, Native::filter(input, callback_output, throwing, ID::InferencePolicy::PRESERVE, false, tiny()))
  TEST_FALSE(fs::exists(callback_output))
  Native::ScanOptions options;
  options.buffering = tiny();
  Native::MatchCallback throwing_scan = [](const std::string&, const std::vector<Native::MatchRecord>&) { throw std::runtime_error("callback"); };
  TEST_EXCEPTION(std::runtime_error, Native::scan(input, options, {}, throwing_scan))
  ID bad;
  auto& run = bad.addRun("bad");
  auto source = run.addSource(ID::SourceFile {});
  auto query = run.addIdentification(source, ID::Observation {});
  ID::MatchData bytes;
  bytes.representation = std::string("invalid\xc0\x80", 9);
  run.addMatch(query, bytes);
  TEST_EXCEPTION(Exception::InvalidValue, Native::store(text_output, bad, tiny()))
  TEST_FALSE(fs::exists(text_output))
  auto small = tiny();
  small.max_record_bytes = 8;
  TEST_EXCEPTION(Exception::InvalidValue, Native::store(size_output, fixture.data, small))
  TEST_FALSE(fs::exists(size_output))
  fs::remove_all(input);
}
END_SECTION

START_SECTION((invalid descriptors fail transactionally))
{
  Fixture fixture;
  std::string directory;
  NEW_TMP_FILE(directory);
  fs::remove_all(directory);
  Native::store(directory, fixture.data, tiny());
  ID destination;
  destination.addRun("original");
  const auto original_uuid = destination.getRuns().front().getUuid();
  replaceText(fs::path(directory) / "manifest.json", "\"next_match_id\": 4", "\"next_match_id\": 2");
  TEST_EXCEPTION(Exception::InvalidValue, Native::load(directory, destination, tiny()))
  TEST_EQUAL(destination.getRuns().front().getUuid(), original_uuid)
  replaceText(fs::path(directory) / "manifest.json", "\"next_match_id\": 2", "\"next_match_id\": 4");
  replaceText(fs::path(directory) / "manifest.json", "\"unit_type\": 2", "\"unit_type\": 4294967298");
  TEST_EXCEPTION(Exception::InvalidValue, Native::load(directory, destination, tiny()))
  TEST_EQUAL(destination.getRuns().front().getUuid(), original_uuid)
  fs::remove_all(directory);
}
END_SECTION

START_SECTION((empty datasets and absent parent catalogues remain distinct))
{
  ID empty;
  std::string directory;
  NEW_TMP_FILE(directory);
  fs::remove_all(directory);
  Native::store(directory, empty, tiny());
  ID loaded;
  loaded.addRun("replace");
  Native::load(directory, loaded);
  TEST_TRUE(loaded.getRuns().empty())
  TEST_TRUE(loaded.getInferenceResults().empty())
  TEST_TRUE(Native::inspect(directory).empty())
  Native::ScanOptions scan;
  const auto statistics = Native::scan(directory, scan, {}, {});
  TEST_EQUAL(statistics.queries, 0)
  fs::remove_all(directory);
}
END_SECTION

END_TEST
