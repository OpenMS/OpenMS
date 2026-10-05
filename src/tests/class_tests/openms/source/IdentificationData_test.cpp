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
START_SECTION((dataset primary score contract and atomic switching))
{
  ID data;
  ID::Run a("A"), b("B"), incompatible("bad");
  auto raw = score();
  auto pep = score("PEP", false);
  auto ar = a.addScore(raw);
  auto ap = a.addScore(pep);
  auto bp = b.addScore(pep); // Different local column order is valid.
  auto br = b.addScore(raw);
  a.setPrimaryScore(ar);
  b.setPrimaryScore(br);
  auto aq = a.addIdentification(a.addSource({}), {});
  auto bq = b.addIdentification(b.addSource({}), {});
  a.addMatch(aq, peptide(), {4.0, 0.1});
  auto bm = b.addMatch(bq, peptide(), {std::nullopt, 5.0});
  data.addRun(a);
  data.addRun(b);
  data.addRun("placeholder");
  TEST_TRUE(data.getPrimaryScoreDefinition() == raw)
  incompatible.setPrimaryScore(incompatible.addScore(score("raw", false)));
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(incompatible))
  TEST_EQUAL(data.getRuns().size(), 3)
  auto replacement = data.getRun("B");
  replacement.setScore(bm, bp, 0.2);
  replacement.setPrimaryScore(bp);
  TEST_EXCEPTION(Exception::InvalidValue, data.replaceRun(replacement))
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
  ID::Run unscored("unscored");
  auto uq = unscored.addIdentification(unscored.addSource({}), {});
  unscored.addMatch(uq, peptide());
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(unscored))
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
  TEST_EXCEPTION(Exception::InvalidValue, run.transformMatches([&](ID::MatchData&) { run = ID::Run("replacement"); }))
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

START_SECTION((pooled inference is owned provenance and survives filtering and replacement))
{
  ID data;
  auto& first = data.addRun("A");
  auto source = first.addSource({});
  first.addScore(score());
  auto query = first.addIdentification(source, {});
  auto removed = first.addMatch(query, peptide(), {1.0});
  auto kept = first.addMatch(query, peptide(), {2.0});
  const auto uuid = first.getUuid();
  const auto before_provenance = first;
  data.addRun("B");
  ID::InferenceResult inference;
  inference.identifier = "pooled";
  inference.inputs.push_back({"A", uuid, score(), {removed, kept, removed}, true, "all candidates"});
  inference.inputs.push_back({"B", data.getRun("B").getUuid(), std::nullopt, {}, false, {}});
  inference.assignments.push_back({"A", uuid, removed, 0, {{"db", "protein"}}});
  inference.assignments.push_back({"A", uuid, {8000}, 0, {}});
  data.addInferenceResult(inference);
  TEST_EQUAL(data.getRun("A").getNextMatchId(), 8001)
  data.getRun("A") = before_provenance;
  TEST_EQUAL(data.getRun("A").getNextMatchId(), 8001)
  TEST_EXCEPTION(Exception::InvalidValue, data.getRun("A").importMatch(query, {8000}, peptide(), {1.0}))
  auto old = data.getRun("A");
  const auto* address = &data.getRun("A");
  TEST_EQUAL(data.filterMatches([&](const ID::Match& match) { return match.getId() == kept; }, ID::InferencePolicy::PRESERVE), 1)
  TEST_EQUAL(data.getInferenceResults().size(), 1)
  TEST_TRUE(address == &data.getRun("A"))
  TEST_EQUAL(data.getInferenceResults()[0].inputs[0].matches.size(), 3)
  TEST_TRUE(data.getRun("A").findMatch(removed) == nullptr)
  TEST_EQUAL(data.getInferenceResults()[0].assignments[1].parents.size(), 0)
  old.filterMatches([&](const ID::Match& match) { return match.getId() == kept; });
  data.replaceRun(old);
  TEST_EQUAL(data.getRun("A").getNextMatchId(), 8001)
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
  inference.inputs.push_back({"absent", absent.getUuid(), {}, {{42}}, true, {}});
  data.addInferenceResult(inference);
  TEST_EXCEPTION(Exception::InvalidValue, data.addRun(absent))
  data.validate();
}
END_SECTION
END_TEST
