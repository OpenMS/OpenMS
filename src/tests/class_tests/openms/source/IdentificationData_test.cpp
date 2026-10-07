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
  payload.sequence_evidence.push_back({database_id, "protein", 1, 100, "K", "R"});
  payload.setMetaValue("list", StringList {"alpha", "beta"});
  auto copied = run.importMatch(q, ID::MatchId {8}, payload, {2.0});
  auto moved = run.importMatch(q, ID::MatchId {9}, std::move(payload), {3.0});
  TEST_EQUAL(run.getIdentification(q).data_id, std::string(100, 'q'))
  TEST_EQUAL(run.getIdentification(q).getMetaValue("observation"), "kept")
  TEST_EQUAL(run.getMatch(copied).representation, std::string(100, 'A'))
  TEST_EQUAL(run.getMatch(moved).sequence_evidence[0].accession, "protein")
  TEST_EQUAL(run.getMatch(moved).getMetaValue("list"), run.getMatch(copied).getMetaValue("list"))
  run.validate(); // missing supplementary values remain valid
  auto& values = const_cast<std::vector<double>&>(run.getMatch(moved).getScoreValues());
  values[1] = std::numeric_limits<double>::infinity();
  TEST_EXCEPTION(Exception::InvalidValue, run.validate())
  values[1] = std::numeric_limits<double>::quiet_NaN();
  values[0] = std::numeric_limits<double>::quiet_NaN();
  TEST_EXCEPTION(Exception::InvalidValue, run.validate())
  values[0] = 3.0;
  values.pop_back();
  TEST_EXCEPTION(Exception::InvalidValue, run.validate())
  values.push_back(std::numeric_limits<double>::quiet_NaN());
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

START_SECTION((transform validation preserves original payloads and optional adduct ownership))
{
  ID::Run run("compounds", ID::MoleculeKind::COMPOUND);
  auto source = run.addSource({});
  auto query = run.addIdentification(source, {});
  ID::MatchData data;
  data.representation = "CCO";
  data.encoding = ID::Encoding::SMILES;
  data.charge = 1;
  data.adduct = AdductInfo::parseAdductString("M+H;1+");
  auto first = run.addMatch(query, data);
  auto second = run.addMatch(query, data);
  data.adduct.reset();
  TEST_TRUE(run.getMatch(first).adduct.has_value())
  Size calls = 0;
  TEST_EXCEPTION(Exception::InvalidValue, run.transformMatches([&](ID::MatchData& match) {
    match.name = "changed";
    if (++calls == 2) match.charge = 2;
  }))
  TEST_EQUAL(run.getMatch(first).name, "")
  TEST_EQUAL(run.getMatch(second).charge, 1)
  run.transformMatches([](ID::MatchData& match) { match.name = "ethanol"; });
  TEST_EQUAL(run.getMatch(first).name, "ethanol")
  ID::Run copy(run);
  auto changed = copy.getMatch(first).getData();
  changed.adduct.reset();
  TEST_EXCEPTION(Exception::InvalidValue, copy.replaceMatch(first, changed))
  copy.replaceMatch(first, changed, {});
  TEST_TRUE(run.getMatch(first).adduct.has_value())
  TEST_TRUE(! copy.getMatch(first).adduct.has_value())
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
  data.adduct = AdductInfo::parseAdductString("M+H;1+");
  const auto retained = run.addMatch(query, data, {2.0});
  const auto other = run.addMatch(other_query, data, {3.0});
  run.setSelectedMatch(query, removed);
  run.setSelectedMatch(other_query, other);
  run.prepareLookupIndexes();
  TEST_EQUAL(run.eraseMatches([&](const ID::Match& match) { return match.getId() == removed; }), 1)
  TEST_TRUE(run.findIdentification(empty) == nullptr)
  TEST_TRUE(! run.getIdentification(query).getSelectedMatch())
  TEST_TRUE(run.getIdentification(other_query).getSelectedMatch() == other)
  TEST_TRUE(run.getMatch(retained).adduct.has_value())
  TEST_REAL_SIMILAR(*run.getScore(retained, primary), 2.0)
  data.adduct.reset();
  run.replaceMatch(retained, data, {4.0});
  TEST_TRUE(! run.getMatch(retained).adduct.has_value())
  data.adduct = AdductInfo::parseAdductString("M+H;1+");
  run.replaceMatch(retained, data, {5.0});
  run.transformMatches([](ID::MatchData& match) { match.name = "ethanol"; });
  TEST_EQUAL(run.getMatch(retained).name, "ethanol")
  TEST_EQUAL(run.getMatch(other).name, "ethanol")
  TEST_TRUE(run.getMatch(retained).adduct.has_value())
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
    TEST_REAL_SIMILAR(*run.getMatch(ids[i]).getScores()[0], static_cast<double>((i / 5) * 10 + i % 5))
  const auto query = run.addIdentification(source, {});
  TEST_EQUAL(run.getIdentification(query).getMatches().size(), 0)
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
END_TEST
