// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/IdentificationDataEdits.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

using namespace OpenMS;
using namespace std;

namespace
{
  using ID = IdentificationData;

  /**
    A consensus map with identification data: run "search" (q-values) with database sequences P1 and P2, identification
    scan=1 (PEPTIDEA/P1, PEPTIDEB/P2) that feature 0 links, scan=2 (PEPTIDEC/P2) and an inference result over the run
  */
  ConsensusMap makeMap()
  {
    ConsensusMap map;
    auto& data = map.getIdentificationData();
    auto& run = data.addRun("search");
    ID::ScoreDefinition score;
    score.name = "q-value";
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    ID::Database database;
    database.path = "proteins.fasta";
    const auto db = run.addDatabase(database);
    run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P1"}, {db, "P2"}});
    const auto source = run.addSource({});
    ID::Observation observation;
    observation.data_id = "scan=1";
    const auto first = run.addIdentification(source, observation);
    observation.data_id = "scan=2";
    const auto second = run.addIdentification(source, observation);
    const auto add = [&](ID::QueryId query, const std::string& sequence, const std::string& accession, double value) {
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      match.sequence_evidence.push_back({db, accession, std::nullopt, std::nullopt, 0, 0});
      return run.addMatch(query, match, {value});
    };
    const auto a = add(first, "PEPTIDEA", "P1", 0.01);
    const auto b = add(first, "PEPTIDEB", "P2", 0.05);
    add(second, "PEPTIDEC", "P2", 0.02);

    ID::InferenceResult result;
    result.identifier = "inference";
    ProteinHit hit;
    hit.setAccession("P1");
    hit.setScore(0.9);
    result.proteins.setScoreType("probability");
    result.proteins.setHits({hit});
    ID::ScoreDefinition probability;
    probability.name = "probability";
    probability.scope = ID::ScoreScope::PROTEIN;
    result.protein_score = probability;
    ID::InferenceInput input;
    input.run_identifier = run.getIdentifier();
    input.run_uuid = run.getUuid();
    result.inputs.push_back(input);
    data.addInferenceResult(result);

    ConsensusFeature feature;
    feature.setUniqueId(1);
    feature.addIDMatch({run.getUuid(), a});
    feature.addIDMatch({run.getUuid(), b});
    map.push_back(feature);
    return map;
  }
}

START_TEST(IdentificationDataEdits, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((explicit IdentificationDataEdits(IdentificationData& data)))
{
  auto map = makeMap();
  IdentificationDataEdits edits(map.getIdentificationData());
  TEST_EQUAL(&edits.getData() == &map.getIdentificationData(), true)
  TEST_EQUAL(edits.inferenceResults().size(), 1)
}
END_SECTION

START_SECTION((bool commit()))
{
  auto map = makeMap();
  auto& data = map.getIdentificationData();
  const auto& run = data.getRuns()[0];
  const auto uuid = run.getUuid();
  const auto& query = run.getSources()[0].identifications[0];
  const auto& match = query.getMatches()[0];
  const auto query_id = query.getId();
  const auto match_id = match.getId();

  // without changes, nothing changes
  {
    IdentificationDataEdits edits(data);
    edits.settings(run);
    edits.observation(run, query);
    edits.match(run, match);
    edits.databaseSequences(run);
    const auto before = data;
    TEST_EQUAL(edits.commit(), false)
    TEST_EQUAL(data == before, true)
  }

  // the edited records replace their originals; IDs and UUIDs stay
  {
    IdentificationDataEdits edits(data);
    edits.settings(run).software = "Comet";
    edits.settings(run).setMetaValue("note", "settings");
    edits.observation(run, query).setMetaValue("note", "identification");
    edits.match(run, match).setMetaValue("note", "match");
    edits.databaseSequences(run)[1].setMetaValue("note", "protein");
    edits.inferenceResults()[0].proteins.setScoreType("posterior probability");
    TEST_EQUAL(edits.commit(), true)
  }
  const auto& edited = data.getRuns()[0];
  TEST_EQUAL(edited.getUuid(), uuid)
  TEST_EQUAL(edited.getSettings().software, "Comet")
  TEST_EQUAL(edited.getSettings().getMetaValue("note"), "settings")
  TEST_EQUAL(edited.getIdentification(query_id).getMetaValue("note"), "identification")
  TEST_EQUAL(edited.getMatch(match_id).getMetaValue("note"), "match")
  TEST_EQUAL(edited.getScore(match_id, *edited.getPrimaryScore()).value(), 0.01)
  TEST_EQUAL((*edited.getDatabaseSequences())[1].getMetaValue("note"), "protein")
  // the protein score definition does not apply to proteins of another score type
  TEST_EQUAL(data.getInferenceResults()[0].proteins.getScoreType(), "posterior probability")
  TEST_EQUAL(data.getInferenceResults()[0].protein_score.has_value(), false)
  // the feature still links the matches
  const auto linked = map[0].getLinkedIdentifications(data);
  ABORT_IF(linked.size() != 1)
  TEST_EQUAL(linked[0].matches.size(), 2)
  TEST_EQUAL(linked[0].matches[0]->getMetaValue("note"), "match")
}
END_SECTION

START_SECTION(([EXTRA] commit() changes nothing if an edit is not valid))
{
  auto map = makeMap();
  auto& data = map.getIdentificationData();
  const auto& run = data.getRuns()[0];
  const auto before = data;
  IdentificationDataEdits edits(data);
  edits.match(run, run.getSources()[0].identifications[0].getMatches()[0]).setMetaValue("note", "match");
  // the files of a run are its sources, not settings metadata
  edits.settings(run).setMetaValue("spectra_data", ListUtils::create<std::string>("sample.mzML"));
  TEST_EXCEPTION(Exception::InvalidValue, edits.commit())
  TEST_EQUAL(data == before, true)
}
END_SECTION

START_SECTION((std::vector<IdentificationData::DatabaseSequence>& databaseSequences(const IdentificationData::Run& run)))
{
  ID data;
  data.addRun("search");
  IdentificationDataEdits edits(data);
  TEST_EXCEPTION(Exception::InvalidValue, edits.databaseSequences(data.getRuns()[0]))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
