// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <algorithm>
#include <filesystem>

using namespace OpenMS;
using ID = IdentificationData;
using Inference = IdentificationDataInference;
using Adapter = IdentificationDataAdapter;

namespace
{
void addRun(ID& data, const std::string& name, const std::string& database, bool pep, const std::string& file = "")
{
  auto& run = data.addRun(name);
  ID::RunSettings settings;
  settings.software = "engine";
  run.setSettings(settings);
  ID::Database search_database;
  search_database.path = database;
  const auto database_id = run.addDatabase(search_database);
  ID::ScoreDefinition definition;
  definition.name = pep ? "PEP" : "Posterior Probability";
  definition.higher_better = ! pep;
  definition.software = "engine";
  const auto score = run.addScore(definition);
  run.setPrimaryScore(score);
  ID::DatabaseSequence parent;
  parent.database = database_id;
  parent.accession = "P1";
  parent.sequence = "PEPTIDEOTHER";
  run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {parent});
  ID::SourceFile file_source;
  file_source.path = file;
  const auto source = run.addSource(file_source);
  auto query = run.addIdentification(source, {});
  ID::MatchData match;
  match.representation = database == "dbB" ? "OTHER" : "PEPTIDE";
  match.charge = 2;
  ID::SequenceEvidence evidence;
  evidence.database = database_id;
  evidence.accession = "P1";
  evidence.start = 0;
  evidence.end = 6;
  match.sequence_evidence = {evidence};
  run.addMatch(query, match, {pep ? 0.1 : 0.9});
  match.representation = database == "dbB" ? "PEPTIDE" : "OTHER";
  run.addMatch(query, match, {pep ? 0.5 : 0.5});
}
std::vector<Inference::Input> inputs(const ID& data)
{
  std::vector<Inference::Input> result;
  for (const auto& run : data.getRuns())
    result.push_back({run.getUuid(), *run.getPrimaryScore(),
                      run.getScoreDefinition(*run.getPrimaryScore()).higher_better ? Inference::ProbabilityType::POSTERIOR_PROBABILITY
                                                                                   : Inference::ProbabilityType::POSTERIOR_ERROR_PROBABILITY});
  return result;
}
} // namespace

START_TEST(IdentificationDataInference, "$Id$")

START_SECTION((static IdentificationData::InferenceResult infer(const IdentificationData&, const std::vector<Input>&, const std::string&)))
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbB", false);
  TEST_EXCEPTION(Exception::InvalidValue, Inference::infer(data, inputs(data), "mixed"))
  // Normalize both runs explicitly before combining their score columns.
  ID normalized;
  addRun(normalized, "A", "dbA", true);
  addRun(normalized, "B", "dbB", true);
  data = std::move(normalized);
  ID::ScoreDefinition pp_definition;
  pp_definition.name = "Posterior Probability";
  for (auto& run : data.getRuns())
  {
    // getRuns() is read-only; configure through the owning run accessor.
    auto& mutable_run = data.getRun(run.getIdentifier());
    const auto pp_score = mutable_run.addScore(pp_definition);
    for (const auto& source : mutable_run.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
          mutable_run.setScore(match.getId(), pp_score, 1.0 - *mutable_run.getScore(match.getId(), *mutable_run.getPrimaryScore()));
  }
  const auto pp = data.getRun("B").findScore(pp_definition);
  auto selected = inputs(data);
  // Inference may still explicitly consume a supplementary PP column.
  selected[1].score = pp;
  selected[1].probability = Inference::ProbabilityType::POSTERIOR_PROBABILITY;
  auto result = Inference::infer(data, selected, "pooled");
  TEST_EQUAL(result.inputs.size(), 2)
  TEST_EQUAL(result.inputs[0].run_uuid, data.getRuns()[0].getUuid())
  TEST_EQUAL(result.inputs[1].run_uuid, data.getRuns()[1].getUuid())
  TEST_TRUE(*result.inputs[0].score == data.getRuns()[0].getScoreDefinition(selected[0].score))
  TEST_EQUAL(result.proteins.getHits().size(), 2)
  TEST_NOT_EQUAL(result.proteins.getHits()[0].getAccession(), result.proteins.getHits()[1].getAccession())
  TEST_REAL_SIMILAR(result.proteins.getHits()[0].getScore(), 0.9)
  const auto& run = data.getRuns()[0];
  const auto& matches = run.getSources()[0].identifications[0].getMatches();
  TEST_EQUAL(matches.size(), 2)
  TEST_REAL_SIMILAR(*run.getScore(matches[0].getId(), selected[0].score), 0.1)
  TEST_REAL_SIMILAR(*run.getScore(matches[1].getId(), selected[0].score), 0.5)
  TEST_EQUAL(matches[0].sequence_evidence[0].accession, "P1")
  TEST_EQUAL(matches[1].sequence_evidence[0].accession, "P1")
  TEST_EQUAL(run.getDatabaseSequences()->size(), 1)
  data.addInferenceResult(result);
  data.getRun("A").eraseMatches([](const auto&) { return true; });
  TEST_EQUAL(data.getInferenceResults()[0].inputs[0].run_uuid, data.getRun("A").getUuid())
  TEST_EQUAL(data.getInferenceResults()[0].proteins.getHits().size(), 2)

  auto malformed = selected;
  malformed[0].probability = Inference::ProbabilityType::POSTERIOR_PROBABILITY;
  TEST_EXCEPTION(Exception::InvalidParameter, Inference::infer(data, malformed, "invalid"))
  TEST_EXCEPTION(Exception::InvalidParameter, Inference::infer(data, {selected[1], selected[1]}, "duplicate"))
}
END_SECTION

START_SECTION((static void retainProteins(IdentificationData::InferenceResult&, const std::set<IdentificationData::QualifiedAccession>&)))
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbB", true);
  auto result = Inference::infer(data, inputs(data), "pooled");
  ProteinIdentification::ProteinGroup group;
  group.probability = 0.75;
  for (const auto& hit : result.proteins.getHits())
    group.accessions.push_back(hit.getAccession());
  result.proteins.getProteinGroups().push_back(group);
  const auto original_input = result.inputs[0];
  Inference::retainProteins(result, {{"dbA", "P1"}});
  TEST_EQUAL(result.proteins.getHits().size(), 1)
  TEST_EQUAL(result.proteins.getProteinGroups().size(), 0)
  TEST_EQUAL(result.qualified_accessions.size(), 1)
  TEST_EQUAL(result.inputs[0].run_uuid, original_input.run_uuid)
  TEST_TRUE(result.inputs[0].score == original_input.score)
  TEST_EQUAL(result.inputs[0].selection, original_input.selection)
  TEST_EQUAL(data.getRuns()[1].getSources()[0].identifications[0].getMatches()[0].sequence_evidence.size(), 1)
  TEST_EQUAL(data.getRuns()[1].getDatabaseSequences()->size(), 1)
  data.addInferenceResult(result);
  std::string path;
  NEW_TMP_FILE(path)
  IdentificationDataFile::store(path, data);
  ID loaded;
  IdentificationDataFile::load(path, loaded);
  TEST_EQUAL(loaded.getInferenceResults()[0].proteins.getHits().size(), 1)
  std::filesystem::remove_all(path);
}
END_SECTION

START_SECTION([EXTRA] invalid probabilities are rejected without editing input values)
{
  ID data;
  addRun(data, "A", "dbA", true);
  auto selected = inputs(data);
  auto& run = data.getRun("A");
  const auto match = run.getSources()[0].identifications[0].getMatches()[0].getId();
  run.setScore(match, selected[0].score, 1.1);
  TEST_EXCEPTION(Exception::InvalidParameter, Inference::infer(data, selected, "invalid"))
  TEST_REAL_SIMILAR(*run.getScore(match, selected[0].score), 1.1)
  TEST_EQUAL(data.getInferenceResults().size(), 0)
}
END_SECTION


START_SECTION([EXTRA] pooled inference rejects inconsistent mappings of the same peptidoform)
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbB", true);
  auto& run = data.getRun("B");
  const auto match_id = run.getSources()[0].identifications[0].getMatches()[0].getId();
  auto match = run.getMatch(match_id).getData();
  match.representation = "PEPTIDE";
  run.replaceMatch(match_id, match, {0.1});
  TEST_EXCEPTION(Exception::InvalidParameter, Inference::infer(data, inputs(data), "conflicting"))
  TEST_EQUAL(data.getInferenceResults().size(), 0)
}
END_SECTION


START_SECTION([EXTRA] legacy export writes pooled inference as one merged protein run)
{
  // Two runs of the same search, as IDMerger or ProteinInference would merge them in the legacy model.
  ID data;
  addRun(data, "A", "db.fasta", true, "a.mzML");
  addRun(data, "B", "db.fasta", true, "b.mzML");
  const auto pooled = Inference::infer(data, inputs(data), "pooled");
  TEST_TRUE(pooled.protein_score.has_value())
  data.addInferenceResult(pooled);

  // Strict export represents the bridge output, including its score definitions.
  const auto exported = Adapter::toLegacy(data);
  TEST_EQUAL(exported.losses.size(), 0)
  TEST_EQUAL(exported.proteins.size(), 1)
  ABORT_IF(exported.proteins.size() != 1)
  const auto& merged = exported.proteins[0];
  TEST_EQUAL(merged.getIdentifier(), "pooled")
  StringList paths;
  merged.getPrimaryMSRunPath(paths);
  TEST_EQUAL(ListUtils::concatenate(paths, ","), "a.mzML,b.mzML")
  TEST_EQUAL(merged.getHits().size(), pooled.proteins.getHits().size())
  TEST_EQUAL(exported.peptides.size(), 2)
  for (Size i = 0; i < exported.peptides.size(); ++i)
  {
    TEST_EQUAL(exported.peptides[i].getIdentifier(), "pooled")
    TEST_EQUAL(static_cast<Int>(exported.peptides[i].getMetaValue("id_merge_index")), static_cast<Int>(i))
  }

  // Importing the legacy run (here from idXML) restores the files and the inference score definitions.
  std::string path;
  NEW_TMP_FILE(path)
  IdXMLFile().store(path, exported.proteins, exported.peptides);
  std::vector<ProteinIdentification> stored_proteins;
  PeptideIdentificationList stored_peptides;
  IdXMLFile().load(path, stored_proteins, stored_peptides);
  const auto imported = Adapter::importLegacy(stored_proteins, stored_peptides).data;
  TEST_EQUAL(imported.getRuns().size(), 1)
  TEST_EQUAL(imported.getInferenceResults().size(), 1)
  ABORT_IF(imported.getRuns().size() != 1 || imported.getRuns()[0].getSources().size() != 2 || imported.getInferenceResults().size() != 1)
  const auto& run = imported.getRuns()[0];
  TEST_EQUAL(run.getSources()[0].file.path, "a.mzML")
  TEST_EQUAL(run.getSources()[1].file.path, "b.mzML")
  TEST_EQUAL(run.getSettings().metaValueExists("identification:inference:protein_score:name"), false)
  const auto& restored = imported.getInferenceResults()[0];
  TEST_TRUE(restored.protein_score == pooled.protein_score)
  TEST_TRUE(restored.group_score == pooled.group_score)
  TEST_EQUAL(restored.inputs.size(), 1)
  ABORT_IF(restored.inputs.size() != 1)
  TEST_TRUE(restored.inputs[0].score == pooled.inputs[0].score)
  TEST_EQUAL(restored.proteins.metaValueExists("identification:inference:protein_score:name"), false)
  TEST_EQUAL(Adapter::toLegacy(imported).losses.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] legacy export never attaches a run to the search settings of another run)
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbB", true);
  data.addInferenceResult(Inference::infer(data, inputs(data), "pooled"));
  // Different databases cannot share one legacy protein run.
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(data))
  Adapter::ExportOptions options;
  options.loss_policy = Adapter::LossPolicy::ALLOW;
  const auto exported = Adapter::toLegacy(data, options);
  TEST_EQUAL(exported.proteins.size(), 2)
  ABORT_IF(exported.proteins.size() != 2 || exported.peptides.size() != 2)
  TEST_EQUAL(exported.proteins[0].getIdentifier(), "pooled")
  TEST_EQUAL(exported.proteins[0].getSearchParameters().db, "dbA")
  TEST_EQUAL(exported.proteins[1].getIdentifier(), "B")
  TEST_EQUAL(exported.proteins[1].getSearchParameters().db, "dbB")
  TEST_EQUAL(exported.peptides.size(), 2)
  TEST_EQUAL(exported.peptides[0].getIdentifier(), "pooled")
  TEST_EQUAL(exported.peptides[1].getIdentifier(), "B")
  TEST_TRUE(std::any_of(exported.losses.begin(), exported.losses.end(),
                        [](const auto& loss) { return loss.find("cannot share") != std::string::npos; }))

  // Settings that legacy merging accepts but that differ are a reported loss.
  ID close;
  addRun(close, "A", "db.fasta", true, "a.mzML");
  addRun(close, "B", "db.fasta", true, "b.mzML");
  auto settings = close.getRun("B").getSettings();
  settings.search.missed_cleavages = 2;
  close.getRun("B").setSettings(settings);
  close.addInferenceResult(Inference::infer(close, inputs(close), "pooled"));
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(close))
  const auto merged = Adapter::toLegacy(close, options);
  TEST_EQUAL(merged.proteins.size(), 1)
  ABORT_IF(merged.peptides.size() != 2)
  TEST_EQUAL(merged.peptides[1].getIdentifier(), "pooled")
  TEST_EQUAL(merged.losses.size(), 1)

  // Without primary MS files, the PSMs of the merged runs could not be told apart.
  ID unnamed;
  addRun(unnamed, "A", "db.fasta", true);
  addRun(unnamed, "B", "db.fasta", true);
  unnamed.addInferenceResult(Inference::infer(unnamed, inputs(unnamed), "pooled"));
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(unnamed))
}
END_SECTION

START_SECTION([EXTRA] an evidence - only run can share later database sequences)
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbA", true);
  data.getRun("A").setDatabaseSequences(std::nullopt);
  const auto result = Inference::infer(data, inputs(data), "catalogue");
  TEST_EQUAL(result.proteins.getHits().size(), 1)
  TEST_EQUAL(result.proteins.getHits()[0].getSequence(), "PEPTIDEOTHER")
  TEST_EQUAL(result.inputs.size(), 2)
}
END_SECTION

END_TEST
