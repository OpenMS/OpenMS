// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <bit>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>

using namespace OpenMS;
using ID = IdentificationData;
namespace fs = std::filesystem;

namespace
{
struct RemoveDirectory
{
  fs::path path;
  ~RemoveDirectory()
  {
    std::error_code ec;
    fs::remove_all(path, ec);
  }
};

ProteinIdentification processing()
{
  ProteinIdentification result;
  result.setIdentifier("processing id");
  result.setSearchEngine("custom search");
  result.setSearchEngineVersion("2.3");
  result.setInferenceEngine("pooled inference");
  result.setInferenceEngineVersion("4.5");
  DateTime date;
  date.set("2026-10-05T14:31:22.137");
  result.setDateTime(date);
  result.setScoreType("protein posterior");
  result.setHigherScoreBetter(true);
  result.setSignificanceThreshold(0.99);
  result.setPrimaryMSRunPath({"C:\\input files\\one.mzML", "/second/input.mzML"});
  result.setPrimaryMSRunPath({"C:\\input files\\one.raw"}, true);
  result.setMetaValue("empty", DataValue::EMPTY);
  result.setMetaValue("empty strings", StringList {});
  auto& search = result.getSearchParameters();
  search.db = "custom database";
  search.db_version = "2026-10";
  search.taxonomy = "all";
  search.charges = "2:6";
  search.mass_type = ProteinIdentification::PeakMassType::AVERAGE;
  search.fixed_modifications = {"fixed custom", "fixed custom"};
  search.variable_modifications = {"variable custom"};
  search.missed_cleavages = 2;
  search.fragment_mass_tolerance = 0.01;
  search.fragment_mass_tolerance_ppm = false;
  search.precursor_mass_tolerance = 7.5;
  search.precursor_mass_tolerance_ppm = true;
  search.enzyme_term_specificity = EnzymaticDigestion::SPEC_NOCTERM;
  EmpiricalFormula n_term("H");
  n_term.setCharge(1);
  search.digestion_enzyme = Protease("custom enzyme", "(?<=[KR])", {"synonym A", "synonym B"}, "custom rule", n_term, EmpiricalFormula("OH"),
                                     "MS:1000001", "[KR]|{P}", 12, 13, 14);
  search.setMetaValue("custom_modification_definitions", "1|custom|custom (M)|custom|M|Anywhere||12.34|12.34|");
  DataValue tolerance(4.0);
  tolerance.setUnitType(DataValue::MS_ONTOLOGY);
  tolerance.setUnit(1000040);
  search.setMetaValue("unit value", tolerance);
  return result;
}

ResidueModification modification()
{
  ResidueModification mod;
  mod.setId("custom");
  mod.setFullId("custom (M)");
  mod.setFullName("custom protein modification");
  mod.setName("short custom");
  mod.setPSIMODAccession("MOD:00777");
  mod.setUniModRecordId(777);
  mod.setTermSpecificity(ResidueModification::PROTEIN_N_TERM);
  mod.setOrigin('M');
  mod.setSourceClassification(ResidueModification::CHEMICAL_DERIVATIVE);
  mod.setProvenance(ResidueModification::DEFINED);
  mod.setAverageMass(100.25);
  mod.setMonoMass(99.75);
  mod.setDiffAverageMass(12.5);
  mod.setDiffMonoMass(12.34);
  mod.setFormula("C2 H4");
  EmpiricalFormula diff("C2H4");
  diff.setCharge(-1);
  mod.setDiffFormula(diff);
  mod.setSynonyms({"synonym 1", "synonym 2"});
  EmpiricalFormula loss("H2O");
  loss.setCharge(1);
  mod.setNeutralLossDiffFormulas({loss, EmpiricalFormula("NH3")});
  mod.setNeutralLossMonoMasses({18.01, 17.02});
  mod.setNeutralLossAverageMasses({18.1, 17.1});
  return mod;
}

ID::ScoreDefinition definition()
{
  ID::ScoreDefinition result;
  result.name = "posterior";
  result.accession = "MS:1000000";
  result.scope = ID::ScoreScope::MATCH;
  result.higher_better = true;
  result.software = "inference";
  result.software_version = "3";
  result.calibration = "joint calibration";
  result.aggregation = "pooled";
  result.parameters.setMetaValue("integer list", IntList {1, 2});
  return result;
}

ID::InferenceResult inference(const ID::Run& run)
{
  ID::InferenceResult result;
  result.identifier = "pooled result";
  result.proteins = processing();
  result.parent_score = definition();
  result.parent_score->scope = ID::ScoreScope::PROTEIN;
  result.group_score = definition();
  result.group_score->scope = ID::ScoreScope::PROTEIN_GROUP;
  ProteinHit a(0.9, 2, "A", "PEPTIDE");
  a.setCoverage(42.5);
  a.setDescription(" original description ");
  a.setMetaValue("target_decoy", "target");
  a.setMetaValue("empty list", DoubleList {});
  DataValue unit(12);
  unit.setUnitType(DataValue::UNIT_ONTOLOGY);
  unit.setUnit(10);
  a.setMetaValue("integer unit", unit);
  std::set<std::pair<Size, ResidueModification>> mods {{0, modification()}};
  a.setModifications(mods);
  result.proteins.insertHit(a);
  ProteinHit b(0.8, 3, "B", "OTHER");
  ResidueModification incomplete;
  incomplete.setId("name without full id");
  mods = {{1, incomplete}};
  b.setModifications(mods);
  result.proteins.insertHit(b);
  ProteinIdentification::ProteinGroup group;
  group.probability = 0.95;
  group.accessions = {"B", "only in group", "A", "A"};
  DataArrays::FloatDataArray floats;
  floats.setName("float values");
  floats.insert(floats.end(), {1.25f, -2.5f});
  floats.setMetaValue("unit", unit);
  auto processed = std::make_shared<DataProcessing>();
  processed->setCompletionTime(processing().getDateTime());
  processed->setProcessingActions({DataProcessing::IDENTIFICATION, DataProcessing::FILTERING});
  processed->getSoftware().setName("array producer");
  processed->getSoftware().setVersion("1.2");
  processed->getSoftware().setMetaValue("flag", DataValue::EMPTY);
  processed->setMetaValue("weights", DoubleList {0.25, 0.75});
  floats.setDataProcessing({processed});
  group.getFloatDataArrays().push_back(floats);
  DataArrays::IntegerDataArray integers;
  integers.setName("integer values");
  integers.insert(integers.end(), {2, -4});
  integers.setMetaValue("empty integer list", IntList {});
  group.getIntegerDataArrays().push_back(integers);
  DataArrays::StringDataArray strings;
  strings.setName("string values");
  strings.insert(strings.end(), {"a", "", "b"});
  strings.setMetaValue("string-list", StringList {"x", "y"});
  group.getStringDataArrays().push_back(strings);
  result.proteins.insertProteinGroup(group);
  group = {};
  group.probability = 0.5;
  result.proteins.insertIndistinguishableProteins(group);
  result.qualified_accessions = {{"A", {"db", "parent A"}}, {"B", {"db", "parent B"}}, {"only in group", {"other db", "group member"}}};
  ID::InferenceInput input;
  input.run_identifier = run.getIdentifier();
  input.run_uuid = run.getUuid();
  input.score = definition();
  input.selection = "top candidates, before filtering";
  result.inputs.push_back(input);
  ID::Run absent("absent run");
  input.run_identifier = absent.getIdentifier();
  input.run_uuid = absent.getUuid();
  input.score.reset();
  input.selection = "No candidates passed the input filter";
  result.inputs.push_back(input);
  input.selection = "Imported legacy run-level provenance";
  result.inputs.push_back(input);
  return result;
}

void replaceJsonNumber(const fs::path& path, const std::string& key, const std::string& replacement)
{
  std::ifstream in(path, std::ios::binary);
  std::string content {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
  in.close();
  auto position = content.find('"' + key + '"');
  if (position == std::string::npos) throw std::runtime_error("Missing JSON fixture key");
  position = content.find(':', position);
  position = content.find_first_not_of(" \t\r\n", position + 1);
  const auto end = content.find_first_of(",}\r\n \t", position);
  content.replace(position, end - position, replacement);
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << content;
}
} // namespace

START_TEST(IdentificationDataFileInference, "$Id$")

START_SECTION((inline group membership obeys the whole-record size limit))
{
  ID data;
  ID::InferenceResult result;
  result.identifier = "large group";
  ProteinIdentification::ProteinGroup group;
  group.accessions.assign(20, std::string(64, 'A'));
  result.proteins.insertProteinGroup(group);
  data.addInferenceResult(result);
  std::string path, rejected, filtered;
  NEW_TMP_FILE(path)
  NEW_TMP_FILE(rejected)
  NEW_TMP_FILE(filtered)
  RemoveDirectory cleanup {path}, cleanup_rejected {rejected}, cleanup_filtered {filtered};
  IdentificationDataFile::Options small;
  small.max_record_bytes = 256;
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::store(rejected, data, small))
  TEST_FALSE(fs::exists(rejected))
  IdentificationDataFile::store(path, data);
  ID destination;
  destination.addRun("unchanged");
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::load(path, destination, small))
  TEST_EQUAL(destination.getRuns().size(), 1)
  TEST_EQUAL(destination.getRuns().front().getIdentifier(), "unchanged")
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::filter(path, filtered,
    [](const auto&, const auto&) { return true; }, ID::InferencePolicy::PRESERVE, false, small))
  TEST_FALSE(fs::exists(filtered))
}
END_SECTION

START_SECTION((multiple inference results share files without mixing rows or metadata))
{
  ID data;
  for (Size i = 0; i < 5; ++i)
  {
    ID::InferenceResult result;
    result.identifier = "result" + std::to_string(i);
    ProteinHit hit(double(i), 1, "P" + std::to_string(i), "PEPTIDE");
    hit.setMetaValue("key" + std::to_string(i), static_cast<int>(i));
    result.proteins.insertHit(hit);
    ProteinIdentification::ProteinGroup group;
    group.probability = i * 0.1;
    group.accessions = {hit.getAccession()};
    result.proteins.getProteinGroups().push_back(group);
    data.addInferenceResult(result);
  }
  std::string path, filtered;
  NEW_TMP_FILE(path)
  NEW_TMP_FILE(filtered)
  IdentificationDataFile::Options options;
  options.row_group_rows = 2;
  IdentificationDataFile::store(path, data, options);
  TEST_TRUE(fs::exists(fs::path(path) / "proteins.parquet"))
  TEST_FALSE(fs::exists(fs::path(path) / "proteins-1.parquet"))
  ID loaded;
  IdentificationDataFile::load(path, loaded, options);
  TEST_EQUAL(loaded.getInferenceResults().size(), 5)
  for (Size i = 0; i < 5; ++i)
    TEST_TRUE(loaded.getInferenceResults()[i].proteins == data.getInferenceResults()[i].proteins)
  IdentificationDataFile::filter(path, filtered, [](const auto&, const auto&) { return true; }, ID::InferencePolicy::PRESERVE);
  IdentificationDataFile::load(filtered, loaded);
  for (Size i = 0; i < 5; ++i)
    TEST_TRUE(loaded.getInferenceResults()[i].proteins == data.getInferenceResults()[i].proteins)
  fs::remove_all(path);
  fs::remove_all(filtered);
}
END_SECTION

START_SECTION((typed pooled inference preserves complete protein values, group arrays and ordered provenance))
{
  ID data;
  auto& run = data.addRun("run with processing");
  run.setProcessingMetadata(processing());
  ID::ParentRecord parent;
  parent.identity = {"db", "parent A"};
  parent.target_decoy = ID::TargetDecoy::BOTH;
  parent.sequence = "PEPTIDE";
  parent.description = " original parent ";
  parent.setMetaValue("empty", DataValue::EMPTY);
  run.setParents(std::vector<ID::ParentRecord> {parent});
  const auto input = inference(run);
  data.addInferenceResult(input);
  std::string path;
  NEW_TMP_FILE(path)
  RemoveDirectory cleanup {path};
  IdentificationDataFile::Options options;
  options.batch_rows = 1;
  options.row_group_rows = 2;
  IdentificationDataFile::store(path, data, options);
  ID loaded;
  IdentificationDataFile::load(path, loaded, options);
  TEST_EQUAL(loaded.getInferenceResults().size(), 1)
  const auto& result = loaded.getInferenceResults().at(0);
  TEST_EQUAL(result.identifier, input.identifier)
  TEST_TRUE(result.proteins == input.proteins)
  TEST_TRUE(result.parent_score == input.parent_score)
  TEST_TRUE(result.group_score == input.group_score)
  TEST_TRUE(result.qualified_accessions == input.qualified_accessions)
  TEST_EQUAL(result.inputs.size(), 3)
  for (Size i = 0; i < result.inputs.size(); ++i)
  {
    TEST_EQUAL(result.inputs[i].run_uuid, input.inputs[i].run_uuid)
    TEST_EQUAL(result.inputs[i].run_identifier, input.inputs[i].run_identifier)
    TEST_TRUE(result.inputs[i].score == input.inputs[i].score)
    TEST_EQUAL(result.inputs[i].selection, input.inputs[i].selection)
  }

  TEST_TRUE(loaded.getRun(run.getIdentifier()).getProcessingMetadata() == run.getProcessingMetadata())
  const auto& parents = *loaded.getRun(run.getIdentifier()).getParents();
  TEST_EQUAL(parents.size(), 1)
  TEST_TRUE(parents[0].identity == parent.identity)
  TEST_TRUE(parents[0].target_decoy == parent.target_decoy)
  TEST_EQUAL(parents[0].sequence, parent.sequence)
  TEST_EQUAL(parents[0].description, parent.description)
  TEST_TRUE(static_cast<const MetaInfoInterface&>(parents[0]) == static_cast<const MetaInfoInterface&>(parent))
  TEST_EQUAL(result.proteins.getHits()[1].getModifications().begin()->second.getFullId(), "")
  TEST_EQUAL(result.proteins.getHits()[0].getModifications().begin()->second.getProvenance(), ResidueModification::DEFINED)
  TEST_EQUAL(loaded.getRun(run.getIdentifier()).getNextMatchId(), 1)
  TEST_TRUE(fs::exists(fs::path(path) / "proteins.parquet"))
  TEST_FALSE(fs::exists(fs::path(path) / "inference.json"))
}
END_SECTION

START_SECTION((inference floating point bit patterns survive dictionary encoding and row groups))
{
  const std::vector<UInt64> bits {0, 0x8000000000000000ULL, 0x7ff8000000000001ULL, 0x7ff8000000000002ULL};
  const std::vector<UInt32> float_bits {0, 0x80000000U, 0x7fc00001U, 0x7fc00002U};
  ID data;
  ID::InferenceResult result;
  result.identifier = "float bits";
  for (Size i = 0; i < bits.size(); ++i)
  {
    ProteinHit protein(std::bit_cast<double>(bits[i]), 0, "P" + std::to_string(i), "PEPTIDE");
    protein.setCoverage(std::bit_cast<double>(bits[(i + 1) % bits.size()]));
    ResidueModification mod;
    mod.setDiffMonoMass(std::bit_cast<double>(bits[i]));
    mod.setNeutralLossMonoMasses({std::bit_cast<double>(bits[i])});
    std::set<std::pair<Size, ResidueModification>> modifications {{0, mod}};
    protein.setModifications(modifications);
    result.proteins.insertHit(protein);
    ProteinIdentification::ProteinGroup group;
    group.probability = std::bit_cast<double>(bits[i]);
    DataArrays::FloatDataArray values;
    for (auto raw : float_bits)
      values.push_back(std::bit_cast<float>(raw));
    group.getFloatDataArrays().push_back(values);
    result.proteins.insertProteinGroup(group);
  }
  data.addInferenceResult(result);
  std::string path;
  NEW_TMP_FILE(path)
  RemoveDirectory cleanup {path};
  IdentificationDataFile::Options options;
  options.batch_rows = 2;
  options.row_group_rows = 3;
  IdentificationDataFile::store(path, data, options);
  ID loaded;
  IdentificationDataFile::load(path, loaded, options);
  const auto& output = loaded.getInferenceResults().at(0).proteins;
  TEST_EQUAL(output.getHits().size(), bits.size())
  TEST_EQUAL(output.getProteinGroups().size(), bits.size())
  for (Size i = 0; i < bits.size(); ++i)
  {
    TEST_EQUAL(std::bit_cast<UInt64>(output.getHits()[i].getScore()), bits[i])
    TEST_EQUAL(std::bit_cast<UInt64>(output.getHits()[i].getCoverage()), bits[(i + 1) % bits.size()])
    const auto& mod = output.getHits()[i].getModifications().begin()->second;
    TEST_EQUAL(std::bit_cast<UInt64>(mod.getDiffMonoMass()), bits[i])
    TEST_EQUAL(std::bit_cast<UInt64>(mod.getNeutralLossMonoMasses()[0]), bits[i])
    TEST_EQUAL(std::bit_cast<UInt64>(output.getProteinGroups()[i].probability), bits[i])
    const auto& values = output.getProteinGroups()[i].getFloatDataArrays()[0];
    TEST_EQUAL(values.size(), float_bits.size())
    for (Size j = 0; j < values.size(); ++j)
      TEST_EQUAL(std::bit_cast<UInt32>(values[j]), float_bits[j])
  }
}
END_SECTION

START_SECTION((invalid configuration integers fail transactionally and detached aliases are explicit errors))
{
  ID data;
  data.addRun("run").setProcessingMetadata(processing());
  std::string path;
  NEW_TMP_FILE(path)
  RemoveDirectory cleanup {path};
  IdentificationDataFile::store(path, data);
  replaceJsonNumber(fs::path(path) / "manifest.json", "missed_cleavages", "-1");
  ID destination;
  destination.addRun("sentinel");
  const auto sentinel_uuid = destination.getRun("sentinel").getUuid();
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::load(path, destination))
  TEST_EQUAL(destination.getRuns().size(), 1)
  TEST_EQUAL(destination.getRun("sentinel").getUuid(), sentinel_uuid)
  replaceJsonNumber(fs::path(path) / "manifest.json", "missed_cleavages", "4294967296");
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::load(path, destination))
  replaceJsonNumber(fs::path(path) / "manifest.json", "missed_cleavages", "1.5");
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::load(path, destination))
  TEST_EQUAL(destination.getRun("sentinel").getUuid(), sentinel_uuid)
  ID::InferenceResult detached;
  detached.identifier = "detached alias";
  detached.qualified_accessions["unused"] = {"db", "A"};
  data.addInferenceResult(detached);
  std::string output;
  NEW_TMP_FILE(output)
  RemoveDirectory cleanup_output {output};
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::store(output, data))
  TEST_FALSE(fs::exists(output))
}
END_SECTION

START_SECTION((empty inference tables and optional parent catalogue have distinct representations))
{
  ID data;
  data.addRun("no catalogue");
  data.addRun("empty catalogue").setParents(std::vector<ID::ParentRecord> {});
  ID::InferenceResult inference_result;
  inference_result.identifier = "empty result";
  data.addInferenceResult(inference_result);
  std::string path;
  NEW_TMP_FILE(path)
  RemoveDirectory cleanup {path};
  IdentificationDataFile::store(path, data);
  ID loaded;
  IdentificationDataFile::load(path, loaded);
  TEST_FALSE(loaded.getRun("no catalogue").getParents().has_value())
  TEST_TRUE(loaded.getRun("empty catalogue").getParents().has_value())
  TEST_EQUAL(loaded.getRun("empty catalogue").getParents()->size(), 0)
  TEST_EQUAL(loaded.getInferenceResults().size(), 1)
  TEST_EQUAL(loaded.getInferenceResults()[0].inputs.size(), 0)
  TEST_TRUE(loaded.getInferenceResults()[0].proteins == inference_result.proteins)
  for (const std::string table : {"inputs", "proteins", "groups"})
  {
    TEST_TRUE(fs::exists(fs::path(path) / (table + ".parquet")))
  }
  TEST_FALSE(fs::exists(fs::path(path) / "group_members.parquet"))
  TEST_FALSE(fs::exists(fs::path(path) / "input_members.parquet"))
  TEST_FALSE(fs::exists(fs::path(path) / "assignments.parquet"))
  fs::remove(fs::path(path) / "groups.parquet");
  const auto saved_uuid = loaded.getRun("no catalogue").getUuid();
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::load(path, loaded))
  TEST_EQUAL(loaded.getRun("no catalogue").getUuid(), saved_uuid)
  std::string filtered;
  NEW_TMP_FILE(filtered)
  RemoveDirectory cleanup_filtered {filtered};
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::filter(
                                            path, filtered, [](const std::string&, const IdentificationDataFile::MatchRecord&) { return true; },
                                            ID::InferencePolicy::PRESERVE))
  TEST_FALSE(fs::exists(filtered))
}
END_SECTION

START_SECTION((streaming preservation validates live allocation counters and copied parent payloads))
{
  ID data;
  auto& run = data.addRun("run");
  ID::ScoreDefinition score;
  score.name = "score";
  run.setPrimaryScore(run.addScore(score));
  const auto source = run.addSource({});
  const auto query = run.addIdentification(source, {});
  ID::MatchData peptide;
  peptide.representation = "PEPTIDE";
  peptide.charge = 2;
  const auto live = run.addMatch(query, peptide, {1.0});
  TEST_EQUAL(live.value, 1)
  ID::ParentRecord parent;
  parent.identity = {"db", "parent"};
  parent.sequence.assign(5000, 'A');
  run.setParents(std::vector<ID::ParentRecord> {parent});
  ID::InferenceResult result;
  result.identifier = "retained inference";
  ID::InferenceInput input;
  input.run_uuid = run.getUuid();
  input.run_identifier = run.getIdentifier();
  result.inputs.push_back(input);
  data.addInferenceResult(result);
  std::string path;
  NEW_TMP_FILE(path)
  RemoveDirectory cleanup {path};
  IdentificationDataFile::store(path, data);
  replaceJsonNumber(fs::path(path) / "manifest.json", "next_match_id", "1");
  std::string filtered;
  NEW_TMP_FILE(filtered)
  RemoveDirectory cleanup_filtered {filtered};
  auto keep = [](const std::string&, const IdentificationDataFile::MatchRecord&) { return true; };
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::filter(path, filtered, keep, ID::InferencePolicy::PRESERVE))
  TEST_FALSE(fs::exists(filtered))
  replaceJsonNumber(fs::path(path) / "manifest.json", "next_match_id", "2");
  IdentificationDataFile::Options options;
  options.max_record_bytes = 1024;
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataFile::filter(path, filtered, keep, ID::InferencePolicy::PRESERVE, false, options))
  TEST_FALSE(fs::exists(filtered))
}
END_SECTION

END_TEST
