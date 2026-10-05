// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <algorithm>
#include <filesystem>

using namespace OpenMS;
using ID = IdentificationData;
using Inference = IdentificationDataInference;

namespace
{
void addRun(ID& data, const std::string& name, const std::string& database, bool pep)
{
  auto& run = data.addRun(name);
  ProteinIdentification processing;
  processing.setIdentifier(name);
  processing.getSearchParameters().db = database;
  run.setProcessingMetadata(processing);
  ID::ScoreDefinition definition;
  definition.name = pep ? "PEP" : "Posterior Probability";
  definition.higher_better = ! pep;
  const auto score = run.addScore(definition);
  run.setPrimaryScore(score);
  ID::ParentRecord parent;
  parent.identity = {database, "P1"};
  parent.sequence = "PEPTIDEOTHER";
  run.setParents(std::vector<ID::ParentRecord> {parent});
  const auto source = run.addSource({});
  auto query = run.addIdentification(source, {});
  ID::MatchData match;
  match.representation = database == "dbB" ? "OTHER" : "PEPTIDE";
  match.charge = 2;
  ID::ParentEvidence evidence;
  evidence.parent = parent.identity;
  evidence.start = 0;
  evidence.end = 6;
  match.parent_evidence = {evidence};
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
  const auto selected = inputs(data);
  auto result = Inference::infer(data, selected, "pooled");
  TEST_EQUAL(result.inputs.size(), 2)
  TEST_EQUAL(result.inputs[0].matches.size(), 2)
  TEST_EQUAL(result.inputs[1].matches.size(), 2)
  TEST_TRUE(result.inputs[0].membership_known)
  TEST_TRUE(*result.inputs[0].score == data.getRuns()[0].getScoreDefinition(selected[0].score))
  TEST_EQUAL(result.proteins.getHits().size(), 2)
  TEST_NOT_EQUAL(result.proteins.getHits()[0].getAccession(), result.proteins.getHits()[1].getAccession())
  TEST_REAL_SIMILAR(result.proteins.getHits()[0].getScore(), 0.9)
  TEST_EQUAL(result.assignments.size(), 4)
  TEST_EQUAL(result.assignments[0].parents.size(), 1)
  TEST_EQUAL(result.assignments[1].parents.size(), 0)
  TEST_EQUAL(result.assignments[2].parents.size(), 1)
  TEST_EQUAL(result.assignments[3].parents.size(), 0)
  const auto& run = data.getRuns()[0];
  TEST_REAL_SIMILAR(*run.getScore(result.inputs[0].matches[0], selected[0].score), 0.1)
  TEST_EQUAL(run.getMatch(result.inputs[0].matches[0]).parent_evidence[0].parent.accession, "P1")
  TEST_EQUAL(run.getParents()->size(), 1)
  data.addInferenceResult(result);
  data.getRun("A").eraseMatches([](const auto&) { return true; });
  TEST_EQUAL(data.getInferenceResults()[0].inputs[0].matches.size(), 2)
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
  const auto original_members = result.inputs[0].matches;
  Inference::retainProteins(result, {{"dbA", "P1"}});
  TEST_EQUAL(result.proteins.getHits().size(), 1)
  TEST_EQUAL(result.proteins.getProteinGroups().size(), 0)
  TEST_EQUAL(result.qualified_accessions.size(), 1)
  TEST_EQUAL(result.assignments.size(), 4)
  TEST_EQUAL(result.assignments[2].parents.size(), 0)
  TEST_TRUE(result.inputs[0].matches == original_members)
  TEST_EQUAL(data.getRuns()[1].getSourceBlocks()[0].identifications[0].getMatches()[0].parent_evidence.size(), 1)
  TEST_EQUAL(data.getRuns()[1].getParents()->size(), 1)
  data.addInferenceResult(result);
  std::string path;
  NEW_TMP_FILE(path)
  IdentificationDataFile::store(path, data);
  ID loaded;
  IdentificationDataFile::load(path, loaded);
  TEST_EQUAL(loaded.getInferenceResults()[0].proteins.getHits().size(), 1)
  TEST_EQUAL(loaded.getInferenceResults()[0].assignments.size(), 4)
  TEST_TRUE(loaded.getInferenceResults()[0].assignments[2].parents.empty())
  std::filesystem::remove_all(path);
}
END_SECTION

START_SECTION([EXTRA] invalid probabilities are rejected without editing input values)
{
  ID data;
  addRun(data, "A", "dbA", true);
  auto selected = inputs(data);
  auto& run = data.getRun("A");
  const auto match = run.getSourceBlocks()[0].identifications[0].getMatches()[0].getId();
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
  const auto match_id = run.getSourceBlocks()[0].identifications[0].getMatches()[0].getId();
  auto match = run.getMatch(match_id).getData();
  match.representation = "PEPTIDE";
  run.replaceMatch(match_id, match, {0.1});
  TEST_EXCEPTION(Exception::InvalidParameter, Inference::infer(data, inputs(data), "conflicting"))
  TEST_EQUAL(data.getInferenceResults().size(), 0)
}
END_SECTION


START_SECTION([EXTRA] an evidence - only run can share a later parent catalogue)
{
  ID data;
  addRun(data, "A", "dbA", true);
  addRun(data, "B", "dbA", true);
  data.getRun("A").setParents(std::nullopt);
  const auto result = Inference::infer(data, inputs(data), "catalogue");
  TEST_EQUAL(result.proteins.getHits().size(), 1)
  TEST_EQUAL(result.proteins.getHits()[0].getSequence(), "PEPTIDEOTHER")
  TEST_EQUAL(result.inputs.size(), 2)
}
END_SECTION

END_TEST
