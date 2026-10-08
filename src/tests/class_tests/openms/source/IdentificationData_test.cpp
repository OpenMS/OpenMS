// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

using namespace OpenMS;
using ID = IdentificationData;

namespace
{
ID::MatchData peptide(const std::string& sequence = "PEPTIDE")
{
  ID::MatchData data;
  data.representation = sequence;
  data.charge = 2;
  return data;
}
ID::ScoreDefinition score(const std::string& name = "raw", bool higher = true)
{
  ID::ScoreDefinition result;
  result.name = name;
  result.higher_better = higher;
  return result;
}
} // namespace

START_TEST(IdentificationData, "$Id$")
START_SECTION((consuming imports preserve ownership and dense validation rejects malformed scores))
{
  ID::Run run("import");
  auto primary = run.addScore(score());
  run.addScore(score("optional"));
  run.setPrimaryScore(primary);
  auto source = run.addSource({});
  ID::Database database;
  database.path = "database";
  const auto database_id = run.addDatabase(database);
  ID::Observation observation;
  observation.data_id = std::string(100, 'q');
  observation.setMetaValue("observation", "kept");
  auto q = run.importIdentification(source, ID::QueryId {7}, std::move(observation));
  auto payload = peptide(std::string(100, 'A'));
  payload.sequence_evidence.push_back({database_id, "protein", 1, 100, 'K', 'R'});
  payload.setMetaValue("list", StringList {"alpha", "beta"});
  auto copied = run.importMatch(q, ID::MatchId {8}, payload, {2.0});
  auto moved = run.importMatch(q, ID::MatchId {9}, std::move(payload), {3.0});
  TEST_EQUAL(run.getIdentification(q).data_id, std::string(100, 'q'))
  TEST_EQUAL(run.getIdentification(q).getMetaValue("observation"), "kept")
  TEST_EQUAL(run.getMatch(copied).representation, std::string(100, 'A'))
  TEST_EQUAL(run.getMatch(moved).sequence_evidence[0].accession, "protein")
  TEST_EQUAL(run.getMatch(moved).getMetaValue("list"), run.getMatch(copied).getMetaValue("list"))
  run.validate(); // missing supplementary values remain valid
  // Scores live in the run's columns; edits are checked, so the columns cannot become invalid.
  TEST_EXCEPTION(Exception::InvalidValue, run.setScore(moved, run.getScoreId(1), std::numeric_limits<double>::infinity()))
  TEST_EXCEPTION(Exception::InvalidValue, run.setScore(moved, run.getScoreId(0), std::nullopt))
  TEST_REAL_SIMILAR(*run.getScores(moved)[0], 3.0)
  TEST_FALSE(run.getScores(moved)[1].has_value())
  run.validate();
  TEST_EXCEPTION(Exception::InvalidValue, run.importMatch(q, moved, peptide(), {4.0}))
  TEST_EQUAL(run.getNumberOfMatches(), 2)
  run.importMatch(q, ID::MatchId {4}, peptide(), {1.0});
  run.importIdentification(source, ID::QueryId {1}, {});
  run.validate(); // imported IDs need not follow scientific record order
  auto& matches = const_cast<std::vector<ID::Match>&>(run.getIdentification(q).getMatches());
  matches[1] = matches[0];
  TEST_EXCEPTION(Exception::InvalidValue, run.validate())
}
END_SECTION

START_SECTION((dataset ordered score schema and atomic switching))
{
  ID data;
  ID::Run a("A"), b("B"), incompatible("bad");
  auto raw = score();
  auto pep = score("PEP", false);
  auto ar = a.addScore(raw);
  auto ap = a.addScore(pep);
  auto br = b.addScore(raw);
  auto bp = b.addScore(pep);
  a.setPrimaryScore(ar);
  b.setPrimaryScore(br);
  auto aq = a.addIdentification(a.addSource({}), {});
  auto bq = b.addIdentification(b.addSource({}), {});
  a.addMatch(aq, peptide(), {4.0, 0.1});
  auto bm = b.addMatch(bq, peptide(), {5.0, std::nullopt});
  data.addRun(a);
  data.addRun(b);
  data.addRun("placeholder");
  TEST_TRUE(data.getPrimaryScoreDefinition() == raw)
  TEST_TRUE(data.getScoreDefinitions() == a.getScoreDefinitions())
  incompatible.setPrimaryScore(incompatible.addScore(score("raw", false)));
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(incompatible))
  TEST_EQUAL(data.getRuns().size(), 3)
  ID mixed(data);
  mixed.getRun("B").setScore(bm, bp, 0.2);
  mixed.getRun("B").setPrimaryScore(bp);
  TEST_EXCEPTION(Exception::InvalidValue, mixed.validate())
  TEST_EXCEPTION(Exception::InvalidValue, data.setPrimaryScore(pep))
  TEST_TRUE(data.getRun("A").getPrimaryScore() == ar)
  TEST_TRUE(data.getRun("B").getPrimaryScore() == br)
  data.getRun("B").setScore(bm, bp, 0.2);
  data.setPrimaryScore(pep);
  TEST_TRUE(data.getPrimaryScoreDefinition() == pep)
  TEST_TRUE(data.getRun("A").getPrimaryScore() == ap)
  TEST_TRUE(data.getRun("B").getPrimaryScore() == bp)
  data.validate();
  // Mutable construction can be incomplete, but cannot cross validated boundaries.
  data.getRun("A").setPrimaryScore(ar);
  TEST_EXCEPTION(Exception::InvalidValue, data.validate())
  TEST_EXCEPTION(Exception::InvalidValue, data.getPrimaryScoreDefinition())
  data.setPrimaryScore(raw);
  data.validate();
  ID copy(data);
  copy.setPrimaryScore(pep);
  TEST_TRUE(data.getPrimaryScoreDefinition() == raw)
  TEST_TRUE(copy.getPrimaryScoreDefinition() == pep)
  auto changed = raw;
  changed.calibration = "different";
  ID::Run provenance("provenance");
  provenance.setPrimaryScore(provenance.addScore(changed));
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(provenance))
  ID::Run reordered("reordered");
  reordered.addScore(pep);
  reordered.setPrimaryScore(reordered.addScore(raw));
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(reordered))
  ID::Run missing("missing");
  missing.setPrimaryScore(missing.addScore(raw));
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(missing))
  ID extra(data);
  extra.getRun("B").addScore(score("extra"));
  TEST_EXCEPTION(Exception::InvalidValue, extra.validate())
  ID::Run supplementary("supplementary provenance");
  supplementary.setPrimaryScore(supplementary.addScore(raw));
  auto other_pep = pep;
  other_pep.calibration = "different";
  supplementary.addScore(other_pep);
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(supplementary))
  auto mismatched = data;
  mismatched.getRun("B").addScore(score("extra"));
  TEST_EXCEPTION(Exception::InvalidValue, mismatched.getScoreDefinitions())
  TEST_EXCEPTION(Exception::InvalidValue, mismatched.setPrimaryScore(pep))
  TEST_TRUE(mismatched.getRun("A").getPrimaryScore() == ar)
  TEST_TRUE(mismatched.getRun("B").getPrimaryScore() == br)
  ID::Run unscored("unscored");
  auto uq = unscored.addIdentification(unscored.addSource({}), {});
  unscored.addMatch(uq, peptide());
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(unscored))
  ID empty;
  TEST_EXCEPTION(Exception::InvalidValue, empty.addRun(unscored))
}
END_SECTION

START_SECTION((run - local scores, primary coverage and schema guards))
{
  ID::Run run("search");
  auto source = run.addSource({});
  auto raw = run.addScore(score());
  auto query = run.addIdentification(source, {});
  auto a = run.addMatch(query, peptide(), {3.0});
  auto b = run.addMatch(query, peptide("OTHER"), {2.0});
  run.setPrimaryScore(raw);
  TEST_REAL_SIMILAR(*run.getScore(a, raw), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, run.replaceMatch(a, peptide("DIFFERENT")))
  TEST_EXCEPTION(Exception::InvalidValue, run.replaceMatch(a, peptide("DIFFERENT"), {}))
  run.replaceMatch(a, peptide("DIFFERENT"), {4.0});
  TEST_REAL_SIMILAR(*run.getScore(a, raw), 4.0)
  run.replaceMatch(a, peptide(), {3.0});
  TEST_EXCEPTION(Exception::InvalidValue, run.addMatch(query, peptide()))
  TEST_EXCEPTION(Exception::InvalidValue, run.setScore(a, raw, std::nullopt))
  TEST_EXCEPTION(Exception::InvalidValue, run.setScore(a, raw, std::numeric_limits<double>::infinity()))
  auto bound = run.bindScore(raw);
  TEST_REAL_SIMILAR(*bound(run.getMatch(a)), 3.0)
  auto qvalue = run.addScore(score("q-value", false));
  TEST_EXCEPTION(Exception::InvalidValue, bound(run.getMatch(a)))
  TEST_REAL_SIMILAR(*run.getScore(a, raw), 3.0)
  TEST_EXCEPTION(Exception::InvalidValue, run.setPrimaryScore(qvalue))
  run.setScore(a, qvalue, 0.0);
  run.setScore(b, qvalue, 0.02);
  run.setPrimaryScore(qvalue);
  TEST_TRUE(run.getPrimaryScore() == qvalue)
  TEST_REAL_SIMILAR(*run.getScore(a, qvalue), 0.0)
  ID::Run copy(run);
  TEST_EQUAL(copy.getUuid(), run.getUuid())
  TEST_REAL_SIMILAR(*copy.getScore(a, raw), 3.0)
  auto copy_score = copy.addScore(score("copy-only"));
  auto original_score = run.addScore(score("original-only"));
  TEST_EQUAL(copy_score.value, original_score.value)
  TEST_EXCEPTION(Exception::InvalidValue, run.getScore(a, copy_score))
  TEST_EXCEPTION(Exception::InvalidValue, copy.getScore(a, original_score))
  auto foreign_source = copy.addSource({});
  run.addSource({});
  TEST_EXCEPTION(Exception::InvalidValue, run.addIdentification(foreign_source, {}))
  ID::Run independent("independent");
  auto foreign_score = independent.addScore(score());
  TEST_EXCEPTION(Exception::InvalidValue, run.getScore(a, foreign_score))
  run.validate();
  copy.validate();
}
END_SECTION

START_SECTION((filtering is atomic, keeps stable IDs and clears removed selections))
{
  ID::Run run("search");
  auto source = run.addSource({});
  auto raw = run.addScore(score());
  auto query = run.addIdentification(source, {});
  auto a = run.addMatch(query, peptide(), {1.0});
  auto b = run.addMatch(query, peptide(), {2.0});
  auto empty = run.addIdentification(source, {});
  run.setSelectedMatch(query, b);
  TEST_EXCEPTION(Exception::InvalidValue, run.setSelectedMatch(empty, b))
  Size calls = 0;
  TEST_EXCEPTION(std::runtime_error, run.filterMatches([&](const ID::Match&) {
    if (++calls == 2) throw std::runtime_error("stop");
    return false;
  }))
  TEST_EQUAL(run.getNumberOfMatches(), 2)
  TEST_TRUE(run.getIdentification(query).getSelectedMatch() == b)
  TEST_EXCEPTION(Exception::InvalidValue, run.filterMatches([&](const ID::Match&) {
    run.setScore(a, raw, 9.0);
    return true;
  }))
  TEST_EXCEPTION(Exception::InvalidValue, run.transformMatches([&](ID::MatchData&) { run.eraseMatches([](const ID::Match&) { return true; }); }))
  TEST_REAL_SIMILAR(*run.getScore(a, raw), 1.0)
  auto next = run.getNextMatchId();
  TEST_EQUAL(run.filterMatches([&](const ID::Match& match) { return match.getId() == a; }, true), 1)
  TEST_EQUAL(run.getNumberOfIdentifications(), 2)
  TEST_TRUE(run.findIdentification(empty) != nullptr)
  TEST_TRUE(run.findMatch(b) == nullptr)
  TEST_EXCEPTION(Exception::InvalidValue, run.importMatch(query, b, peptide(), {2.0}))
  TEST_TRUE(! run.getIdentification(query).getSelectedMatch())
  TEST_EQUAL(run.getNextMatchId(), next)
  auto c = run.addMatch(query, peptide(), {0.0});
  TEST_TRUE(c.value > b.value)
  TEST_TRUE(run.getMatch(a).getId() == a)
  TEST_EQUAL(run.retainBest(raw), 1)
  TEST_EQUAL(run.getNumberOfIdentifications(), 1)
  TEST_TRUE(run.findIdentification(empty) == nullptr)
  TEST_TRUE(run.findMatch(a) != nullptr)
  run.validate();
}
END_SECTION

START_SECTION((Size Run::eraseIdentifications(const std::function<bool(const Identification&)>& remove)))
{
  ID::Run run("search");
  auto first_file = run.addSource({});
  auto second_file = run.addSource({});
  auto raw = run.addScore(score());
  auto kept = run.addIdentification(first_file, {});
  auto a = run.addMatch(kept, peptide("PEPTIDEA"), {1.0});
  auto removed = run.addIdentification(first_file, {});
  auto b = run.addMatch(removed, peptide("PEPTIDEB"), {2.0});
  auto c = run.addMatch(removed, peptide("PEPTIDEC"), {3.0});
  auto empty_kept = run.addIdentification(second_file, {});
  auto empty_removed = run.addIdentification(second_file, {});
  auto last = run.addIdentification(second_file, {});
  auto d = run.addMatch(last, peptide("PEPTIDED"), {4.0});
  run.setSelectedMatch(removed, c);
  const auto view = run.bindScore(raw);
  // a throwing predicate changes nothing
  TEST_EXCEPTION(std::runtime_error, run.eraseIdentifications([](const ID::Identification& query) -> bool {
    if (query.getMatches().empty()) throw std::runtime_error("stop");
    return true;
  }))
  TEST_EQUAL(run.getNumberOfIdentifications(), 5)
  TEST_EQUAL(run.getNumberOfMatches(), 4)
  const auto revision = run.getRevision();
  TEST_EQUAL(run.eraseIdentifications([](const ID::Identification&) { return false; }), 0)
  TEST_EQUAL(run.getRevision(), revision)
  TEST_EQUAL(run.eraseIdentifications([&](const ID::Identification& query) { return query.getId() == removed || query.getId() == empty_removed; }), 2)
  TEST_TRUE(run.getRevision() > revision)
  TEST_EQUAL(run.getNumberOfIdentifications(), 3)
  TEST_EQUAL(run.getNumberOfMatches(), 2)
  TEST_TRUE(run.findIdentification(removed) == nullptr && run.findIdentification(empty_removed) == nullptr)
  TEST_TRUE(run.findMatch(b) == nullptr && run.findMatch(c) == nullptr)
  TEST_TRUE(run.findIdentification(empty_kept) != nullptr)
  TEST_EQUAL(run.getSources()[0].identifications.size(), 1)
  TEST_EQUAL(run.getSources()[1].identifications.size(), 2)
  // the remaining matches keep their scores; views of the earlier state are rejected
  TEST_REAL_SIMILAR(*run.getScore(a, raw), 1.0)
  TEST_REAL_SIMILAR(*run.getScore(d, raw), 4.0)
  TEST_EXCEPTION(Exception::InvalidValue, view(run.getMatch(d)))
  TEST_REAL_SIMILAR(*run.bindScore(raw)(run.getMatch(d)), 4.0)
  // IDs are not reused
  TEST_TRUE(run.addIdentification(first_file, {}).value > last.value)
  run.validate();
}
END_SECTION

START_SECTION((transform validation preserves original payloads and optional adduct ownership))
{
  ID::Run run("compounds", ID::MoleculeKind::COMPOUND);
  auto source = run.addSource({});
  auto query = run.addIdentification(source, {});
  ID::MatchData data;
  data.representation = "CCO";
  data.encoding = ID::Encoding::SMILES;
  data.charge = 1;
  data.details.emplace().adduct = AdductInfo::parseAdductString("M+H;1+");
  auto first = run.addMatch(query, data);
  auto second = run.addMatch(query, data);
  data.details.reset();
  TEST_TRUE(run.getMatch(first).details.value_or_default().adduct.has_value())
  Size calls = 0;
  TEST_EXCEPTION(Exception::InvalidValue, run.transformMatches([&](ID::MatchData& match) {
    match.details.emplace().name = "changed";
    if (++calls == 2) match.charge = 2;
  }))
  TEST_EQUAL(run.getMatch(first).details.value_or_default().name, "")
  TEST_EQUAL(run.getMatch(second).charge, 1)
  run.transformMatches([](ID::MatchData& match) { match.details.emplace().name = "ethanol"; });
  TEST_EQUAL(run.getMatch(first).details.value_or_default().name, "ethanol")
  ID::Run copy(run);
  auto changed = copy.getMatch(first).getData();
  changed.details->adduct.reset();
  TEST_EXCEPTION(Exception::InvalidValue, copy.replaceMatch(first, changed))
  copy.replaceMatch(first, changed, {});
  TEST_TRUE(run.getMatch(first).details.value_or_default().adduct.has_value())
  TEST_TRUE(! copy.getMatch(first).details.value_or_default().adduct.has_value())
  TEST_EXCEPTION(Exception::InvalidValue, run.addMatch(query, peptide()))
}
END_SECTION

START_SECTION((editing mixed optional adducts preserves identities scores and selections))
{
  ID::Run run("compounds", ID::MoleculeKind::COMPOUND);
  const auto primary = run.addScore(score());
  run.setPrimaryScore(primary);
  const auto source = run.addSource({});
  const auto query = run.addIdentification(source, {});
  const auto empty = run.addIdentification(source, {});
  const auto other_query = run.addIdentification(run.addSource({}), {});
  ID::MatchData data;
  data.representation = "CCO";
  data.encoding = ID::Encoding::SMILES;
  data.charge = 1;
  const auto removed = run.addMatch(query, data, {1.0});
  data.details.emplace().adduct = AdductInfo::parseAdductString("M+H;1+");
  const auto retained = run.addMatch(query, data, {2.0});
  const auto other = run.addMatch(other_query, data, {3.0});
  run.setSelectedMatch(query, removed);
  run.setSelectedMatch(other_query, other);
  run.prepareLookupIndexes();
  TEST_EQUAL(run.eraseMatches([&](const ID::Match& match) { return match.getId() == removed; }), 1)
  TEST_TRUE(run.findIdentification(empty) == nullptr)
  TEST_TRUE(! run.getIdentification(query).getSelectedMatch())
  TEST_TRUE(run.getIdentification(other_query).getSelectedMatch() == other)
  TEST_TRUE(run.getMatch(retained).details.value_or_default().adduct.has_value())
  TEST_REAL_SIMILAR(*run.getScore(retained, primary), 2.0)
  data.details.reset();
  run.replaceMatch(retained, data, {4.0});
  TEST_TRUE(! run.getMatch(retained).details.value_or_default().adduct.has_value())
  data.details.emplace().adduct = AdductInfo::parseAdductString("M+H;1+");
  run.replaceMatch(retained, data, {5.0});
  run.transformMatches([](ID::MatchData& match) { match.details.emplace().name = "ethanol"; });
  TEST_EQUAL(run.getMatch(retained).details.value_or_default().name, "ethanol")
  TEST_EQUAL(run.getMatch(other).details.value_or_default().name, "ethanol")
  TEST_TRUE(run.getMatch(retained).details.value_or_default().adduct.has_value())
  TEST_REAL_SIMILAR(*run.getScore(retained, primary), 5.0)
  TEST_REAL_SIMILAR(*run.getScore(other, primary), 3.0)
  run.validate();
}
END_SECTION

START_SECTION((run settings remain owned through copy move and replacement))
{
  ID::Run run("search");
  ID::RunSettings settings;
  settings.software = "original";
  run.setSettings(settings);
  ID::Run copy(run);
  TEST_TRUE(copy == run)
  settings.software = "replacement";
  copy.setSettings(settings);
  TEST_EQUAL(run.getSettings().software, "original")
  TEST_EQUAL(copy.getSettings().software, "replacement")
  ID::Run moved(std::move(copy));
  TEST_EQUAL(moved.getSettings().software, "replacement")
  copy.setSettings(run.getSettings());
  TEST_EQUAL(copy.getSettings().software, "original")
  ID::Run again(run);
  TEST_TRUE(again == run)
  again.setSettings(again.getSettings());
  TEST_TRUE(again == run)
  // The files of a run are its sources; only the raw files behind them may be listed here.
  settings.setMetaValue("spectra_data_raw", StringList {"raw.raw"});
  again.setSettings(settings);
  settings.setMetaValue("spectra_data", StringList {"sample.mzML"});
  TEST_EXCEPTION(Exception::InvalidValue, again.setSettings(settings))
  TEST_EQUAL(again.getSettings().metaValueExists("spectra_data_raw"), true)
  // Runs are edited in place; a copy can never be written back over a run.
  static_assert(! std::is_copy_assignable_v<ID::Run> && ! std::is_move_assignable_v<ID::Run>);
}
END_SECTION

START_SECTION((persistent import preserves order and reserves IDs without renumbering))
{
  ID::Run run("restored");
  auto source = run.addSource({});
  auto query = run.importIdentification(source, {50}, {});
  run.importIdentification(source, {2}, {});
  auto high = run.importMatch(query, {900}, peptide());
  auto low = run.importMatch(query, {3}, peptide());
  TEST_EQUAL(run.getIdentification(query).getMatches()[0].getId().value, 900)
  TEST_EQUAL(run.getIdentification(query).getMatches()[1].getId().value, 3)
  TEST_EXCEPTION(Exception::InvalidValue, run.importMatch(query, high, peptide()))
  TEST_EXCEPTION(Exception::InvalidValue, run.importIdentification(source, {0}, {}))
  const std::string uuid = "12345678-1234-4234-8234-123456789abc";
  run.restoreIdentity(uuid, 100, 1000);
  TEST_EQUAL(run.getUuid(), uuid)
  TEST_EQUAL(run.getNextQueryId(), 100)
  TEST_EQUAL(run.getNextMatchId(), 1000)
  TEST_EXCEPTION(Exception::InvalidValue, run.restoreIdentity(uuid, 1, 1))
  run.filterMatches([&](const ID::Match& match) { return match.getId() == low; });
  TEST_EQUAL(run.getMatch(low).getId().value, 3)
  TEST_EQUAL(run.addMatch(query, peptide()).value, 1000)
  run.reserveMatchId({5000});
  TEST_EQUAL(run.addMatch(query, peptide()).value, 5001)
  TEST_EXCEPTION(Exception::InvalidValue, run.reserveMatchId({std::numeric_limits<UInt64>::max()}))
  run.validate();
}
END_SECTION

START_SECTION((pooled inference is owned provenance and survives filtering))
{
  ID data;
  auto& first = data.addRun("A");
  auto source = first.addSource({});
  first.setPrimaryScore(first.addScore(score()));
  auto query = first.addIdentification(source, {});
  auto removed = first.addMatch(query, peptide(), {1.0});
  auto kept = first.addMatch(query, peptide(), {2.0});
  const auto uuid = first.getUuid();
  data.addRun("B");
  ID::InferenceResult inference;
  inference.identifier = "pooled";
  inference.inputs.push_back({"A", uuid, score(), "all candidates"});
  inference.inputs.push_back({"B", data.getRun("B").getUuid(), std::nullopt, {}});
  data.addInferenceResult(inference);
  TEST_EQUAL(data.getRun("A").getNextMatchId(), kept.value + 1)
  const auto* address = &data.getRun("A");
  TEST_EQUAL(data.filterMatches([&](const ID::Match& match) { return match.getId() == kept; }, ID::InferencePolicy::PRESERVE), 1)
  TEST_EQUAL(data.getInferenceResults().size(), 1)
  TEST_TRUE(address == &data.getRun("A"))
  TEST_EQUAL(data.getInferenceResults()[0].inputs[0].run_uuid, uuid)
  TEST_EQUAL(data.getInferenceResults()[0].inputs[0].selection, "all candidates")
  TEST_TRUE(data.getRun("A").findMatch(removed) == nullptr)
  // Filtering never lowers the ID counter, so removed IDs are not reused.
  TEST_EQUAL(data.getRun("A").getNextMatchId(), kept.value + 1)
  TEST_EXCEPTION(Exception::InvalidValue, data.filterMatches(
                                            [&](const ID::Match&) {
                                              data = ID();
                                              return true;
                                            },
                                            ID::InferencePolicy::PRESERVE))
  TEST_EQUAL(data.getRun("A").getNumberOfMatches(), 1)
  TEST_EQUAL(data.getInferenceResults().size(), 1)
  data.validate();
  data.filterMatches([](const ID::Match&) { return true; }, ID::InferencePolicy::DISCARD);
  TEST_EQUAL(data.getInferenceResults().size(), 0)
  TEST_EQUAL(data.getRun("A").getNumberOfMatches(), 1)
}
END_SECTION

START_SECTION((absent provenance runs cannot collide with later independent imports))
{
  ID data;
  ID::Run absent("absent");
  ID::InferenceResult inference;
  inference.identifier = "external";
  inference.inputs.push_back({"absent", absent.getUuid(), {}, {}});
  data.addInferenceResult(inference);
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(absent))
  data.validate();
}
END_SECTION

START_SECTION((void Run::shrinkToFit()))
{
  ID::Run run("compact");
  run.setPrimaryScore(run.addScore(score()));
  const auto source = run.addSource({});
  std::vector<ID::MatchId> ids;
  for (Size q = 0; q < 50; ++q)
  {
    const auto query = run.addIdentification(source, {});
    for (Size m = 0; m < 5; ++m)
      ids.push_back(run.addMatch(query, peptide(), {static_cast<double>(q * 10 + m)}));
  }
  run.prepareLookupIndexes();
  const ID::Run before(run);
  run.shrinkToFit();
  TEST_TRUE(run == before)
  const auto& queries = run.getSources().front().identifications;
  TEST_EQUAL(queries.capacity(), queries.size())
  TEST_EQUAL(queries.front().getMatches().capacity(), 5)
  // IDs and the lookup indexes built before stay valid.
  for (Size i = 0; i < ids.size(); ++i)
    TEST_REAL_SIMILAR(*run.getScores(ids[i])[0], static_cast<double>((i / 5) * 10 + i % 5))
  const auto query = run.addIdentification(source, {});
  TEST_EQUAL(run.getIdentification(query).getMatches().size(), 0)
}
END_SECTION

START_SECTION((void Run::removeScore(ScoreId score)))
{
  ID::Run run("columns");
  const auto raw = run.addScore(score("raw"));
  const auto pep = run.addScore(score("pep", false));
  const auto q = run.addScore(score("q", false));
  run.setPrimaryScore(q);
  const auto match = run.addMatch(run.addIdentification(run.addSource({}), {}), peptide("PEPTIDE"), {1.0, 0.1, 0.01});
  const auto view = run.bindScore(q);
  TEST_EXCEPTION(Exception::InvalidValue, run.removeScore(q)) // the primary score
  run.removeScore(raw);
  TEST_EQUAL(run.getScoreDefinitions().size(), 2)
  TEST_EQUAL(run.getScoreDefinitions()[0].name, "pep")
  // The primary score keeps its value; handles of moved and removed scores and views bound before are rejected.
  TEST_EQUAL(run.getScoreDefinition(*run.getPrimaryScore()).name, "q")
  TEST_REAL_SIMILAR(*run.getScore(match, *run.getPrimaryScore()), 0.01)
  TEST_EXCEPTION(Exception::InvalidValue, run.getScore(match, raw))
  TEST_EXCEPTION(Exception::InvalidValue, run.getScore(match, pep))
  TEST_EXCEPTION(Exception::InvalidValue, view(run.getMatch(match)))
  TEST_EQUAL(run.getScores(match).size(), 2)
  TEST_REAL_SIMILAR(*run.getScores(match)[0], 0.1)
  run.validate();

  // In a dataset, a score is removed from every run, or (if primary somewhere) from none.
  ID data;
  for (const auto* name : {"A", "B"})
  {
    auto& other = data.addRun(name);
    other.setPrimaryScore(other.addScore(score("raw")));
    other.addScore(score("pep", false));
    other.addMatch(other.addIdentification(other.addSource({}), {}), peptide("PEPTIDE"), {1.0, 0.1});
  }
  TEST_EXCEPTION(Exception::InvalidValue, data.removeScore(score("raw")))
  TEST_EQUAL(data.getScoreDefinitions().size(), 2)
  data.removeScore(score("pep", false));
  TEST_EQUAL(data.getScoreDefinitions().size(), 1)
  data.validate();
}
END_SECTION

START_SECTION((std::map<MatchReference, QueryReference> eraseMatches(...) and eraseIdentifications(...)))
{
  ID data;
  auto& run = data.addRun("A");
  run.setPrimaryScore(run.addScore(score("raw")));
  const auto source = run.addSource({});
  const auto first = run.addIdentification(source, {});
  const auto second = run.addIdentification(source, {});
  const auto kept = run.addMatch(first, peptide("PEPTIDE"), {1.0});
  const auto erased = run.addMatch(first, peptide("PEPTIDER"), {2.0});
  const auto alone = run.addMatch(second, peptide("PEPTIDEK"), {3.0});
  const auto removed = data.eraseMatches([&](const ID::Run&, const ID::Identification&, const ID::Match& match) {
    return match.getId() == erased || match.getId() == alone;
  });
  TEST_EQUAL(removed.size(), 2)
  TEST_TRUE(removed.at({run.getUuid(), alone}) == (ID::QueryReference {run.getUuid(), second}))
  // Identifications stay, also without matches.
  TEST_EQUAL(run.getNumberOfIdentifications(), 2)
  TEST_EQUAL(run.getNumberOfMatches(), 1)
  TEST_TRUE(run.findMatch(kept) != nullptr)
  const auto gone = data.eraseIdentifications([&](const ID::Run&, const ID::Identification& query) { return query.getId() == first; });
  TEST_TRUE(gone.first == (std::set<ID::QueryReference> {{run.getUuid(), first}}))
  TEST_TRUE(gone.second == (std::set<ID::MatchReference> {{run.getUuid(), kept}}))
  TEST_EQUAL(run.getNumberOfIdentifications(), 1)
  TEST_EQUAL(run.getNumberOfMatches(), 0)
  data.validate();
}
END_SECTION

START_SECTION((std::set<std::string> MatchData::extractProteinAccessionsSet() const))
{
  ID::MatchData match;
  TEST_TRUE(match.extractProteinAccessionsSet().empty())
  for (const auto* accession : {"P2", "P1", "", "P2"})
  {
    ID::SequenceEvidence evidence;
    evidence.accession = accession;
    match.sequence_evidence.push_back(evidence);
  }
  // Distinct and sorted; an empty accession is left out (as PeptideHit::extractProteinAccessionsSet())
  TEST_TRUE(match.extractProteinAccessionsSet() == (std::set<std::string> {"P1", "P2"}))
}
END_SECTION

START_SECTION((DatabaseSequence(DatabaseId database, std::string accession, TargetDecoy target_decoy)))
{
  ID::Run run("catalog");
  ID::Database fasta;
  fasta.path = "bsa_td.fasta";
  const auto db = run.addDatabase(fasta);
  run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P02769|ALBU_BOVIN", ID::TargetDecoy::TARGET},
                                                              {db, "DECOY_P02769|ALBU_BOVIN", ID::TargetDecoy::DECOY},
                                                              {db, "P1"}});
  const auto& sequences = *run.getDatabaseSequences();
  ABORT_IF(sequences.size() != 3)
  ID::DatabaseSequence expected;
  expected.database = db;
  expected.accession = "DECOY_P02769|ALBU_BOVIN";
  expected.target_decoy = ID::TargetDecoy::DECOY;
  TEST_TRUE(sequences[1] == expected)
  TEST_EQUAL(sequences[0].accession, "P02769|ALBU_BOVIN")
  TEST_TRUE(sequences[2].target_decoy == ID::TargetDecoy::UNKNOWN)
  TEST_TRUE(sequences[2].sequence.empty() && sequences[2].description.empty() && sequences[2].isMetaEmpty())
}
END_SECTION

START_SECTION((concurrent const lookups build the lazy indexes once))
{
  ID::Run run("lookups");
  run.setPrimaryScore(run.addScore(score()));
  const auto source = run.addSource({});
  std::vector<ID::MatchId> ids;
  for (Size i = 0; i < 2000; ++i)
    ids.push_back(run.addMatch(run.addIdentification(source, {}), peptide(), {static_cast<double>(i)}));
  const ID::Run copy(run); // fresh copy: no index built yet
  std::vector<int> found(8, 0);
  std::vector<std::thread> threads;
  for (Size t = 0; t < found.size(); ++t)
    threads.emplace_back([&, t] {
      for (const auto& id : ids)
        if (copy.findMatch(id) && copy.getIdentificationForMatch(id).getMatches().size() == 1) ++found[t];
    });
  for (auto& thread : threads)
    thread.join();
  for (int count : found)
    TEST_EQUAL(count, 2000)
}
END_SECTION

START_SECTION((merging appends runs, keeps existing run references valid and is atomic))
{
  ID data;
  auto& kept = data.addRun("A");
  auto primary = kept.addScore(score());
  kept.setPrimaryScore(primary);
  auto query = kept.addIdentification(kept.addSource({}), {});
  kept.addMatch(query, peptide(), {1.0});
  const auto* address = &kept;

  ID other;
  auto& incoming = other.addRun("A"); // repeated display name, independent run
  incoming.setPrimaryScore(incoming.addScore(score()));
  incoming.addMatch(incoming.addIdentification(incoming.addSource({}), {}), peptide("SEQUENCE"), {2.0});
  ID::InferenceResult result;
  result.identifier = "pooled";
  result.inputs.push_back({"A", incoming.getUuid(), {}, {}});
  other.addInferenceResult(result);

  data.merge(other);
  TEST_EQUAL(data.getRuns().size(), 2)
  // The existing run was neither copied nor moved.
  TEST_EQUAL(&data.getRun("A") == address, true)
  TEST_EQUAL(kept.getNumberOfMatches(), 1)
  TEST_EQUAL(data.getRun("A#2").getUuid(), incoming.getUuid())
  // Provenance follows the renamed run.
  TEST_EQUAL(data.getInferenceResults().at(0).inputs.at(0).run_identifier, "A#2")
  // Merging the same data again only re-checks equality and adds nothing.
  data.merge(other);
  TEST_EQUAL(data.getRuns().size(), 2)
  TEST_EQUAL(data.getInferenceResults().size(), 1)
  data.merge(data);
  TEST_EQUAL(data.getRuns().size(), 2)

  // A run with a different score definition is rejected and leaves the dataset unchanged.
  ID conflicting;
  auto& foreign = conflicting.addRun("B");
  foreign.setPrimaryScore(foreign.addScore(score("other", false)));
  foreign.addMatch(foreign.addIdentification(foreign.addSource({}), {}), peptide(), {0.5});
  TEST_EXCEPTION(Exception::InvalidValue, data.merge(conflicting))
  TEST_EQUAL(data.getRuns().size(), 2)
  TEST_EQUAL(&data.getRun("A") == address, true)
}
END_SECTION
START_SECTION(([EXTRA] scores are per-run columns, read through the run or bound views))
{
  // Compact records: no per-match score vector, molecule details and evidence flanks out of line or as characters.
  // (Loose bounds: they also hold for debug builds of other standard libraries.)
  TEST_TRUE(sizeof(ID::Match) <= 192)
  TEST_TRUE(sizeof(ID::SequenceEvidence) <= 80)
  ID::Run run("columns");
  const auto primary = run.addScore(score("raw"));
  run.setPrimaryScore(primary);
  const auto query = run.addIdentification(run.addSource({}), {});
  const auto first = run.addMatch(query, peptide("PEPTIDE"), {1.0});
  const auto second = run.addMatch(query, peptide("PEPTIDER"), {2.0});
  const auto view = run.bindScore(primary);
  TEST_REAL_SIMILAR(*view(run.getMatch(second)), 2.0)
  TEST_REAL_SIMILAR(*run.getScores(run.getMatch(first))[0], 1.0)
  const auto revision = run.getRevision();
  // A copy of a match taken now keeps its row while the run keeps the state of its columns.
  const ID::Match snapshot = run.getMatch(second);
  run.setScore(second, primary, 3.0);
  TEST_TRUE(run.getRevision() > revision)
  TEST_REAL_SIMILAR(*view(snapshot), 3.0)
  // Adding a score starts a new state: the old view and the snapshot are rejected, never read in another row.
  const auto extra = run.addScore(score("extra", false));
  TEST_EXCEPTION(Exception::InvalidValue, view(run.getMatch(second)))
  TEST_EXCEPTION(Exception::InvalidValue, run.getScores(snapshot))
  const auto rebound = run.bindScore(primary);
  TEST_REAL_SIMILAR(*rebound(run.getMatch(second)), 3.0)
  TEST_FALSE(run.bindScore(extra)(run.getMatch(second)).has_value())
  // Filtering compacts the columns; values stay with their matches.
  run.eraseMatches([&](const ID::Match& match) { return match.getId() == first; });
  TEST_EXCEPTION(Exception::InvalidValue, rebound(run.getMatch(second)))
  TEST_REAL_SIMILAR(*run.bindScore(primary)(run.getMatch(second)), 3.0)
  TEST_REAL_SIMILAR(*run.getScore(second, primary), 3.0)
  run.validate();
  // A copied run has columns of its own; views of the original reject its matches, and it keeps the revision.
  const ID::Run copy(run);
  const auto original_view = run.bindScore(primary);
  TEST_EXCEPTION(Exception::InvalidValue, original_view(copy.getMatch(second)))
  TEST_REAL_SIMILAR(*copy.bindScore(primary)(copy.getMatch(second)), 3.0)
  TEST_EQUAL(copy.getRevision(), run.getRevision())
  TEST_TRUE(copy == run)
  // A view shares the columns: it stays readable (in its state) when its run is gone.
  ID::ScoreView survivor;
  ID::Match kept;
  {
    ID::Run temporary(run);
    survivor = temporary.bindScore(primary);
    kept = temporary.getMatch(second);
  }
  TEST_REAL_SIMILAR(*survivor(kept), 3.0)
}
END_SECTION

START_SECTION(([EXTRA] molecule details are stored apart; an absent box equals empty details))
{
  ID::MatchData plain = peptide();
  ID::MatchData empty_details = plain;
  empty_details.details.emplace();
  TEST_TRUE(plain == empty_details)
  TEST_TRUE(plain.details.value_or_default().empty())
  ID::MatchData named = plain;
  named.details.emplace().name = "a peptide";
  TEST_FALSE(plain == named)
  ID::MatchData copy = named;
  copy.details->name = "changed";
  TEST_EQUAL(named.details->name, "a peptide") // copies are deep
  ID::Run run("details");
  run.setPrimaryScore(run.addScore(score()));
  const auto query = run.addIdentification(run.addSource({}), {});
  const auto stored = run.addMatch(query, empty_details, {1.0});
  TEST_FALSE(run.getMatch(stored).details.has_value()) // empty details are not kept
  // Single-character flanks and 32-bit positions.
  ID::SequenceEvidence evidence {run.addDatabase({}), "P1", 3, 9, 'K', 'R'};
  TEST_EQUAL(evidence.before, 'K')
  TEST_EQUAL(*evidence.end, 9)
}
END_SECTION

END_TEST
