// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARecords.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <functional>
#include <limits>
#include <stdexcept>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
// Active in NDEBUG as well as Debug, exercising the production validators.
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif
using namespace OpenDIA;
namespace
{
  const RecordContext ctx{0, 1, 2, 3};
  const RecordLimits cap{100, 100};
  const RecordId large = std::numeric_limits<RecordId>::max();
  CandidateFamily specimen()
  {
    CandidateFamily f;
    ExtractionQuerySpec region;
    region.context = ctx;
    region.mz = {400, 500, EdgePolicy::CLOSED};
    f.region = ExtractionQuery(region);
    f.members = {{0, TargetDecoy::TARGET, -2, "[M-H]-"}, {large, TargetDecoy::DECOY, {}, {}}};
    f.alternatives = {{0, {}}, {1, {0}}, {large, {0, large}}};
    f.fragments = {{0, 7, FragmentRole::SHARED}, {large, 7, FragmentRole::SHARED}, {0, 8, FragmentRole::DISCRIMINATING}};
    f.provenance = {"library/rev", "model/rev"};
    return f;
  }
  EvidenceBlock support()
  {
    EvidencePoint missing;
    EvidencePoint zero;
    zero.precursor_id = large;
    zero.value = 0;
    return EvidenceBlock(ctx, cap, {missing, zero});
  }
  ScoreRecord result()
  {
    return {0, large, -3.5, 7.0, Quantity{0, "detector area"}, ScoreAmbiguity::UNRESOLVED, "raw/method-v0", {0, 1}};
  }
  bool rejects(const std::function<void()>& action)
  {
    try { action(); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }
  ScoreBatch produce()
  {
    auto f = specimen();
    FamilyBatch families(ctx, cap, {f});
    auto evidence = support();
    auto score = result();
    ScoreBatch batch(cap, families, evidence, {score});
    f.members.clear();
    score.evidence_rows.clear();
    families = FamilyBatch(ctx, {}, {});
    evidence = EvidenceBlock(ctx, {}, {});
    return batch;
  }
}
#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIAFamilies, "$Id$")
START_SECTION(private v0 families and raw scores)
#else
int main()
{
#endif
  auto batch = produce();
  CHECK(batch.context() == ctx)
  CHECK(batch.size() == 1)
  CHECK(batch.families().size() == 1)
  CHECK(batch.evidence().size() == 2)
  auto f = batch.families().family(0);
  CHECK(f.family_id == 0)
  CHECK(f.region.specification().mz.lower == 400)
  CHECK(f.members.size() == 2)
  CHECK(f.members[0].precursor_id == 0)
  CHECK(f.members[0].target_decoy == TargetDecoy::TARGET)
  CHECK(f.members[0].charge == -2)
  CHECK(f.members[0].adduct == "[M-H]-")
  CHECK(f.members[1].precursor_id == large)
  CHECK(f.members[1].target_decoy == TargetDecoy::DECOY)
  CHECK(!f.members[1].charge && !f.members[1].adduct)
  CHECK(f.alternatives.size() == 3)
  CHECK(f.alternatives[0].precursor_ids.empty())
  CHECK(f.alternatives[1].precursor_ids == std::vector<RecordId>{0})
  CHECK(f.alternatives[2].precursor_ids == (std::vector<RecordId>{0, large}))
  CHECK(f.fragments[0].fragment_id == f.fragments[1].fragment_id)
  CHECK(f.fragments[2].role == FragmentRole::DISCRIMINATING)
  CHECK(f.completeness == CandidateCompleteness::UNKNOWN)
  CHECK(f.provenance.library == "library/rev" && f.provenance.prediction == "model/rev")
  auto positive_charge = specimen();
  positive_charge.members[0].charge = 2;
  CHECK(FamilyBatch(ctx, cap, {positive_charge}).family(0).members[0].charge == 2)
  const auto retained_inputs = []
  {
    const auto parent = produce();
    return std::pair{parent.families(), parent.evidence()};
  }();
  CHECK(retained_inputs.first.family(0).members[1].precursor_id == large)
  CHECK(retained_inputs.second.point(1).value == 0)
  CHECK(!batch.evidence().point(0).value)
  CHECK(batch.evidence().point(1).value == 0)
  CHECK(batch.score(0).identification == -3.5)
  CHECK(batch.score(0).localization == 7)
  CHECK(batch.score(0).quantity->value == 0)
  CHECK(batch.score(0).quantity->unit == "detector area")
  CHECK(batch.score(0).scoring_method == "raw/method-v0")
  CHECK(batch.score(0).ambiguity == ScoreAmbiguity::UNRESOLVED)
  f.members.clear();
  auto detached = batch.score(0);
  detached.quantity->unit.clear();
  detached.evidence_rows.clear();
  CHECK(batch.families().family(0).members.size() == 2)
  CHECK(batch.score(0).evidence_rows.size() == 2)
  CHECK(batch.score(0).quantity->unit == "detector area")
  auto copy = batch;
  auto moved = std::move(batch);
  CHECK(rejects([&] { batch.context(); }))
  CHECK(rejects([&] { batch.size(); }))
  CHECK(rejects([&] { batch.score(0); }))
  CHECK(rejects([&] { batch.families(); }))
  CHECK(rejects([&] { batch.evidence(); }))
  batch = copy;
  copy = std::move(moved);
  CHECK(batch.score(0).identification == copy.score(0).identification)
  CHECK(rejects([&] { copy.score(1); }))
  auto families = copy.families();
  auto family_copy = families;
  auto family_move = std::move(families);
  CHECK(rejects([&] { families.context(); }))
  CHECK(rejects([&] { families.size(); }))
  CHECK(rejects([&] { families.family(0); }))
  families = family_copy;
  family_copy = std::move(family_move);
  CHECK(families.family(0).members.size() == 2)
  CHECK(rejects([&] { families.family(1); }))
  CHECK(rejects([&] { families.family(std::numeric_limits<OpenMS::Size>::max()); }))
  CHECK(rejects([&] { copy.score(std::numeric_limits<OpenMS::Size>::max()); }))
  FamilyBatch empty(ctx, {}, {});
  ScoreBatch empty_scores({}, empty, EvidenceBlock(ctx, {}, {}), {});
  CHECK(empty_scores.size() == 0 && empty_scores.context() == ctx)
  CHECK(empty_scores.families().size() == 0 && empty_scores.evidence().size() == 0)
  CHECK(rejects([&] { empty.family(0); }))
  CHECK(rejects([&] { empty_scores.score(0); }))
  CHECK(FamilyBatch(ctx, {1, 11}, {specimen()}).size() == 1)
  CHECK(rejects([&] { FamilyBatch x(ctx, {0, 100}, {specimen()}); }))
  CHECK(rejects([&] { FamilyBatch x(ctx, {1, 10}, {specimen()}); }))
  CHECK(rejects([&] { FamilyBatch x({1, 1, 2, 3}, {}, {}); }))
  CHECK(rejects([&] { FamilyBatch x(ctx, cap, {specimen(), specimen()}); }))
  const std::vector<std::function<void(CandidateFamily&)>> bad_families{
    [](auto& x) { auto s = x.region.specification(); ++s.context.run_id; x.region = ExtractionQuery(s); },
    [](auto& x) { auto s = x.region.specification(); ++s.context.lane_id; x.region = ExtractionQuery(s); },
    [](auto& x) { auto s = x.region.specification(); ++s.context.calibration_snapshot_id; x.region = ExtractionQuery(s); },
    [](auto& x) { x.completeness = static_cast<CandidateCompleteness>(99); },
    [](auto& x) { x.members.push_back(x.members[0]); },
    [](auto& x) { x.members[0].target_decoy = static_cast<TargetDecoy>(99); },
    [](auto& x) { x.members[0].charge = 0; },
    [](auto& x) { x.members[0].adduct = ""; },
    [](auto& x) { x.alternatives.push_back(x.alternatives[0]); },
    [](auto& x) { x.alternatives[1].precursor_ids = {42}; },
    [](auto& x) { x.alternatives[1].precursor_ids = {0, 0}; },
    [](auto& x) { x.fragments.push_back(x.fragments[0]); },
    [](auto& x) { x.fragments[0].precursor_id = 42; },
    [](auto& x) { x.fragments[0].role = static_cast<FragmentRole>(99); }
  };
  for (const auto& mutate : bad_families)
  {
    auto x = specimen(); mutate(x);
    CHECK(rejects([&] { FamilyBatch invalid(ctx, cap, {x}); }))
  }
  for (auto completeness : {CandidateCompleteness::UNKNOWN, CandidateCompleteness::PARTIAL, CandidateCompleteness::COMPLETE})
  {
    auto x = specimen(); x.completeness = completeness;
    CHECK(FamilyBatch(ctx, cap, {x}).family(0).completeness == completeness)
  }
  auto other = specimen(); other.family_id = 9;
  CHECK(FamilyBatch(ctx, cap, {specimen(), other}).size() == 2)
  // Family IDs are keys, not positions, and alternative IDs are family-local.
  auto first_family = specimen(); first_family.family_id = large;
  auto second_family = specimen(); second_family.family_id = 9;
  second_family.members[0].precursor_id = 42;
  second_family.members[1].precursor_id = 43;
  second_family.alternatives[1].precursor_ids = {42};
  second_family.alternatives[2].precursor_ids = {42, 43};
  second_family.fragments[0].precursor_id = 42;
  second_family.fragments[1].precursor_id = 43;
  second_family.fragments[2].precursor_id = 42;
  const FamilyBatch multiple(ctx, {2, 22}, {second_family, first_family});
  CHECK(multiple.size() == 2)
  CHECK(multiple.family(0).family_id == 9 && multiple.family(1).family_id == large)
  CHECK(multiple.family(0).alternatives[2].alternative_id == multiple.family(1).alternatives[2].alternative_id)
  CHECK(rejects([&] { FamilyBatch invalid(ctx, {2, 21}, {second_family, first_family}); }))
  // An ID elsewhere in the batch must not satisfy a local member reference.
  auto leaked_family = second_family;
  leaked_family.alternatives[1].precursor_ids = {0};
  CHECK(rejects([&] { FamilyBatch invalid(ctx, cap, {leaked_family, first_family}); }))
  leaked_family = second_family;
  leaked_family.fragments[0].precursor_id = 0;
  CHECK(rejects([&] { FamilyBatch invalid(ctx, cap, {leaked_family, first_family}); }))
  EvidencePoint first_point; first_point.family_id = large; first_point.precursor_id = 0;
  EvidencePoint second_point; second_point.family_id = 9; second_point.precursor_id = 42;
  const EvidenceBlock multi_evidence(ctx, cap, {first_point, second_point});
  auto first_score = result(); first_score.family_id = large; first_score.evidence_rows = {0, 0};
  auto second_score = result(); second_score.family_id = 9; second_score.evidence_rows = {1};
  const ScoreBatch multi_scores({2, 3}, multiple, multi_evidence, {first_score, second_score});
  CHECK(multi_scores.size() == 2 && multi_scores.evidence().size() == 2)
  CHECK(multi_scores.score(0).family_id == large && multi_scores.score(1).family_id == 9)
  CHECK(multi_scores.score(0).alternative_id == multi_scores.score(1).alternative_id)
  CHECK(multi_scores.score(0).evidence_rows == (std::vector<OpenMS::Size>{0, 0}))
  CHECK(multi_scores.score(1).evidence_rows == std::vector<OpenMS::Size>{1})
  CHECK(rejects([&] { ScoreBatch invalid({2, 2}, multiple, multi_evidence, {first_score, second_score}); }))
  auto leaked_point = second_point; leaked_point.precursor_id = 0;
  CHECK(rejects([&] { ScoreBatch invalid(cap, multiple, EvidenceBlock(ctx, cap, {first_point, leaked_point}), {first_score, second_score}); }))
  auto wrong_family_score = second_score; wrong_family_score.evidence_rows = {0};
  CHECK(rejects([&] { ScoreBatch invalid(cap, multiple, multi_evidence, {first_score, wrong_family_score}); }))
  for (auto changed : {RecordContext{0, 9, 2, 3}, RecordContext{0, 1, 9, 3}, RecordContext{0, 1, 2, 9}})
  {
    CHECK(rejects([&] { ScoreBatch invalid(cap, multiple, EvidenceBlock(changed, cap, {first_point, second_point}), {first_score, second_score}); }))
  }
  auto no_members = specimen(); no_members.members.clear(); no_members.fragments.clear(); no_members.alternatives = {{0, {}}};
  CHECK(FamilyBatch(ctx, {1, 1}, {no_members}).family(0).alternatives[0].precursor_ids.empty())
  auto no_alternatives = specimen(); no_alternatives.alternatives.clear();
  CHECK(FamilyBatch(ctx, cap, {no_alternatives}).family(0).alternatives.empty())
  CHECK(rejects([&] { ScoreBatch x({0, 2}, families, support(), {result()}); }))
  CHECK(rejects([&] { ScoreBatch x({1, 1}, families, support(), {result()}); }))
  CHECK(ScoreBatch({1, 2}, families, support(), {result()}).size() == 1)
  CHECK(rejects([&] { ScoreBatch x(cap, families, support(), {result(), result()}); }))
  const std::vector<std::function<void(ScoreRecord&)>> bad_scores{
    [](auto& x) { x.family_id = 42; },
    [](auto& x) { x.alternative_id = 42; },
    [](auto& x) { x.ambiguity = static_cast<ScoreAmbiguity>(99); },
    [](auto& x) { x.scoring_method.clear(); },
    [](auto& x) { x.quantity->unit.clear(); },
    [](auto& x) { x.evidence_rows = {2}; },
    [](auto& x) { x.evidence_rows = {std::numeric_limits<OpenMS::Size>::max()}; }
  };
  for (const auto& mutate : bad_scores)
  {
    auto x = result(); mutate(x);
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, support(), {x}); }))
  }
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
  {
    auto x = result(); x.identification = bad;
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, support(), {x}); }))
    x = result(); x.localization = bad;
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, support(), {x}); }))
    x = result(); x.quantity->value = bad;
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, support(), {x}); }))
  }
  for (auto changed : {RecordContext{0, 9, 2, 3}, RecordContext{0, 1, 9, 3}, RecordContext{0, 1, 2, 9}})
  {
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, EvidenceBlock(changed, {}, {}), {}); }))
  }
  for (bool wrong_family : {false, true})
  {
    EvidencePoint p;
    if (wrong_family) { p.family_id = 9; } else { p.precursor_id = 42; }
    auto x = result(); x.evidence_rows = {0};
    CHECK(rejects([&] { ScoreBatch invalid(cap, families, EvidenceBlock(ctx, cap, {p}), {x}); }))
  }
  auto null_score = result();
  null_score.alternative_id = 0; null_score.evidence_rows.clear();
  null_score.identification.reset(); null_score.localization = 0; null_score.quantity.reset();
  null_score.ambiguity = ScoreAmbiguity::RESOLVED;
  ScoreBatch null_batch({1, 0}, families, EvidenceBlock(ctx, {}, {}), {null_score});
  CHECK(!null_batch.score(0).identification && null_batch.score(0).localization == 0)
  CHECK(!null_batch.score(0).quantity && null_batch.score(0).evidence_rows.empty())
  CHECK(null_batch.score(0).ambiguity == ScoreAmbiguity::RESOLVED)
  auto absent = null_score; absent.localization.reset();
  CHECK(!ScoreBatch(cap, families, support(), {absent}).score(0).localization)
  auto identification_zero = absent;
  identification_zero.identification = 0;
  const ScoreBatch zero_identification(cap, families, support(), {identification_zero});
  CHECK(zero_identification.score(0).identification == 0 && !zero_identification.score(0).localization)
  CHECK(!ScoreBatch(cap, families, support(), {absent}).score(0).identification)
  auto repeated = result(); repeated.evidence_rows = {0, 0};
  CHECK(ScoreBatch(cap, families, support(), {repeated}).score(0).evidence_rows.size() == 2)
  auto dead = families; auto alive = std::move(dead);
  CHECK(rejects([&] { ScoreBatch invalid(cap, dead, support(), {}); }))
  auto dead_evidence = support(); auto alive_evidence = std::move(dead_evidence);
  CHECK(rejects([&] { ScoreBatch invalid(cap, alive, dead_evidence, {}); }))
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
