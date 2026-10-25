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
#include <utility>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
// The installed SDK does not export ClassTest. Keep exactly the same checks active,
// including in NDEBUG builds, without making a copy of production validation.
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif

using namespace OpenDIA;
using OpenMS::Size;

namespace
{
  bool invalid(const std::function<void()>& operation)
  {
    try { operation(); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }
  Coordinates coordinate()
  {
    return {-1.25, 501.23456789, 1234567, 9876543, 2, 490.125, 510.875, std::nullopt, std::nullopt};
  }
  RecordContext context() { return {0, 101, 202, 303}; }
  RawTile produceRaw()
  {
    auto c = coordinate();
    std::vector<RawPoint> points{{c, -10.5}};
    c.mobility = Mobility{0.0, MobilityUnit::MILLISECOND};
    c.faims_voltage = 0.0;
    points.push_back({c, 0.0});
    c.mobility = Mobility{1.125, MobilityUnit::INVERSE_REDUCED};
    c.faims_voltage = -45.0;
    points.push_back({c, 22.5});
    RawTile producer(context(), {3, 3}, std::move(points));
    auto lease = producer;
    return lease;
  }
  HypothesisPage produceHypotheses()
  {
    Hypothesis h{404, 505, 601.125, TargetDecoy::DECOY, {"library:revision-2", "model:sha256-abc"},
                 -2.5, Mobility{0.0, MobilityUnit::INVERSE_REDUCED}};
    Hypothesis target = h;
    target.precursor_id = 506;
    target.target_decoy = TargetDecoy::TARGET;
    target.predicted_rt_seconds.reset();
    target.predicted_mobility.reset();
    target.provenance.prediction.reset();
    HypothesisPage producer(context(), {2, 2}, {h, target}, {0, 2, 2}, {{606, 200.25, 0.0}, {607, 300.5, std::nullopt}});
    auto lease = producer;
    return lease;
  }
  EvidenceBlock produceEvidence()
  {
    EvidencePoint measured{404, 505, 606, coordinate(), 0.0};
    EvidencePoint missing = measured;
    missing.coordinates.mobility = Mobility{0.0, MobilityUnit::MILLISECOND};
    missing.coordinates.faims_voltage = 0.0;
    missing.value.reset();
    EvidenceBlock producer(context(), {2, 2}, {measured, missing});
    return EvidenceBlock(producer);
  }
}

#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIARecords, "$Id$")
START_SECTION(v0 owning records)
#else
int main()
{
#endif
  auto raw = produceRaw();
  auto raw_copy = raw;
  auto raw_last = std::move(raw_copy);
  CHECK(invalid([&] { raw_copy.size(); }))
  CHECK(invalid([&] { raw_copy.context(); }))
  CHECK(invalid([&] { raw_copy.point(0); }))
  raw = RawTile({}, {0, 0}, {});
  CHECK(raw_last.size() == 3)
  const auto c = raw_last.point(0).coordinates;
  CHECK(c.rt_seconds == -1.25 && c.mz == 501.23456789)
  CHECK(c.frame_id == 1234567 && c.scan_id == 9876543 && c.ms_level == 2)
  CHECK(c.isolation_lower_mz == 490.125 && c.isolation_upper_mz == 510.875)
  CHECK(!c.mobility && !c.faims_voltage)
  CHECK(raw_last.point(0).intensity == -10.5)
  CHECK(raw_last.point(1).coordinates.mobility->value == 0.0)
  CHECK(raw_last.point(1).coordinates.mobility->unit == MobilityUnit::MILLISECOND)
  CHECK(raw_last.point(1).coordinates.faims_voltage == 0.0 && raw_last.point(1).intensity == 0.0)
  CHECK(raw_last.point(2).coordinates.mobility->value == 1.125)
  CHECK(raw_last.point(2).coordinates.mobility->unit == MobilityUnit::INVERSE_REDUCED)
  CHECK(raw_last.point(2).coordinates.faims_voltage == -45.0)
  CHECK(raw_last.context().run_id == 101 && raw_last.context().lane_id == 202 && raw_last.context().calibration_snapshot_id == 303)
  auto detached = raw_last.point(0);
  detached.intensity = 999;
  CHECK(detached.intensity == 999)
  CHECK(raw_last.point(0).intensity == -10.5)
  raw_copy = raw_last;
  raw = std::move(raw_copy);
  CHECK(raw.size() == 3 && invalid([&] { raw_copy.size(); }))
  raw = RawTile({}, {}, {});

  auto page = produceHypotheses();
  auto page_copy = page;
  auto page_last = std::move(page_copy);
  page = HypothesisPage({}, {}, {}, {0}, {});
  CHECK(page_last.size() == 2 && page_last.fragmentCount(0) == 2 && page_last.fragmentCount(1) == 0)
  CHECK(page_last.hypothesis(0).family_id == 404 && page_last.hypothesis(0).precursor_id == 505)
  CHECK(page_last.hypothesis(0).precursor_mz == 601.125 && page_last.hypothesis(0).target_decoy == TargetDecoy::DECOY)
  CHECK(page_last.hypothesis(1).target_decoy == TargetDecoy::TARGET && page_last.hypothesis(1).precursor_id == 506)
  CHECK(page_last.hypothesis(0).provenance.library == "library:revision-2")
  CHECK(page_last.hypothesis(0).provenance.prediction == "model:sha256-abc")
  CHECK(page_last.hypothesis(0).predicted_rt_seconds == -2.5 && page_last.hypothesis(0).predicted_mobility->value == 0.0)
  CHECK(!page_last.hypothesis(1).predicted_rt_seconds && !page_last.hypothesis(1).predicted_mobility && !page_last.hypothesis(1).provenance.prediction)
  CHECK(page_last.fragment(0, 0).fragment_id == 606 && page_last.fragment(0, 0).mz == 200.25)
  CHECK(page_last.fragment(0, 0).predicted_intensity == 0.0 && !page_last.fragment(0, 1).predicted_intensity)
  CHECK(page_last.context().calibration_snapshot_id == 303)
  auto detached_hypothesis = page_last.hypothesis(0);
  detached_hypothesis.provenance.library = "changed";
  CHECK(detached_hypothesis.provenance.library == "changed" && page_last.hypothesis(0).provenance.library == "library:revision-2")
  CHECK(invalid([&] { page_copy.size(); }))
  CHECK(invalid([&] { page_copy.context(); }))
  CHECK(invalid([&] { page_copy.hypothesis(0); }))
  CHECK(invalid([&] { page_copy.fragmentCount(0); }))
  CHECK(invalid([&] { page_copy.fragment(0, 0); }))
  page_copy = page_last;
  page = std::move(page_copy);
  CHECK(page.size() == 2 && invalid([&] { page_copy.size(); }))
  page = HypothesisPage({}, {}, {}, {0}, {});

  auto evidence = produceEvidence();
  auto evidence_copy = evidence;
  auto evidence_last = std::move(evidence_copy);
  evidence = EvidenceBlock({}, {}, {});
  CHECK(evidence_last.size() == 2 && evidence_last.point(0).value == 0.0 && !evidence_last.point(1).value)
  CHECK(!evidence_last.point(0).coordinates.mobility && evidence_last.point(1).coordinates.mobility->value == 0.0)
  CHECK(evidence_last.point(1).coordinates.faims_voltage == 0.0)
  CHECK(evidence_last.point(0).family_id == 404 && evidence_last.point(0).precursor_id == 505 && evidence_last.point(0).fragment_id == 606)
  CHECK(evidence_last.context().run_id == 101)
  CHECK(invalid([&] { evidence_copy.size(); }))
  CHECK(invalid([&] { evidence_copy.context(); }))
  CHECK(invalid([&] { evidence_copy.point(0); }))
  evidence_copy = evidence_last;
  evidence = std::move(evidence_copy);
  CHECK(evidence.size() == 2 && invalid([&] { evidence_copy.size(); }))
  evidence = EvidenceBlock({}, {}, {});

  const auto maximum = std::numeric_limits<Size>::max();
  CHECK(invalid([&] { raw_last.point(maximum); }))
  CHECK(invalid([&] { evidence_last.point(2); }))
  CHECK(invalid([&] { page_last.hypothesis(2); }))
  CHECK(invalid([&] { page_last.fragmentCount(maximum); }))
  CHECK(invalid([&] { page_last.fragment(0, maximum); }))
  CHECK(invalid([&] { page_last.fragment(1, 0); }))
  for (const auto limits : {RecordLimits{0, 1}, RecordLimits{1, 0}})
  {
    CHECK(invalid([&] { RawTile({}, limits, {{coordinate(), 0}}); }))
    CHECK(invalid([&] { EvidenceBlock({}, limits, {{1, 2, 3, coordinate(), 0}}); }))
    CHECK(invalid([&] { HypothesisPage({}, limits, {Hypothesis{}}, {0, 1}, {Fragment{}}); }))
  }
  for (unsigned version : {1u, std::numeric_limits<unsigned>::max()})
  {
    RecordContext bad{version, 1, 2, 3};
    CHECK(invalid([&] { RawTile(bad, {}, {}); }))
    CHECK(invalid([&] { EvidenceBlock(bad, {}, {}); }))
    CHECK(invalid([&] { HypothesisPage(bad, {}, {}, {0}, {}); }))
  }
  for (const auto& offsets : std::vector<std::vector<Size>>{{}, {0}, {1, 1, 1}, {0, 2, 1}, {0, maximum, 1}, {0, 0, 0}})
  {
    CHECK(invalid([&] { HypothesisPage({}, {2, 1}, {Hypothesis{}, Hypothesis{}}, offsets, {Fragment{}}); }))
  }
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
  {
    for (auto field : {&Coordinates::rt_seconds, &Coordinates::mz, &Coordinates::isolation_lower_mz, &Coordinates::isolation_upper_mz})
    {
      auto bad_c = coordinate();
      bad_c.*field = bad;
      CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
      CHECK(invalid([&] { EvidenceBlock({}, {1, 1}, {{1, 2, 3, bad_c, 0}}); }))
    }
    auto bad_c = coordinate();
    bad_c.mobility = Mobility{bad, MobilityUnit::MILLISECOND};
    CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
    bad_c.mobility.reset();
    bad_c.faims_voltage = bad;
    CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
    CHECK(invalid([&] { RawTile({}, {1, 1}, {{coordinate(), bad}}); }))
    CHECK(invalid([&] { EvidenceBlock({}, {1, 1}, {{1, 2, 3, coordinate(), bad}}); }))
    auto h = Hypothesis{};
    h.precursor_mz = bad;
    CHECK(invalid([&] { HypothesisPage({}, {1, 0}, {h}, {0, 0}, {}); }))
    h.precursor_mz = 0;
    h.predicted_rt_seconds = bad;
    CHECK(invalid([&] { HypothesisPage({}, {1, 0}, {h}, {0, 0}, {}); }))
    h.predicted_rt_seconds.reset();
    h.predicted_mobility = Mobility{bad, MobilityUnit::INVERSE_REDUCED};
    CHECK(invalid([&] { HypothesisPage({}, {1, 0}, {h}, {0, 0}, {}); }))
    CHECK(invalid([&] { HypothesisPage({}, {1, 1}, {Hypothesis{}}, {0, 1}, {{1, bad, 0}}); }))
    CHECK(invalid([&] { HypothesisPage({}, {1, 1}, {Hypothesis{}}, {0, 1}, {{1, 0, bad}}); }))
  }
  auto bad_c = coordinate();
  bad_c.mobility = Mobility{0, static_cast<MobilityUnit>(99)};
  CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
  CHECK(invalid([&] { EvidenceBlock({}, {1, 1}, {{1, 2, 3, bad_c, 0}}); }))
  auto bad_h = Hypothesis{};
  bad_h.predicted_mobility = bad_c.mobility;
  CHECK(invalid([&] { HypothesisPage({}, {1, 0}, {bad_h}, {0, 0}, {}); }))
  bad_h.predicted_mobility.reset();
  bad_h.target_decoy = static_cast<TargetDecoy>(99);
  CHECK(invalid([&] { HypothesisPage({}, {1, 0}, {bad_h}, {0, 0}, {}); }))
  bad_c = coordinate();
  bad_c.ms_level = 0;
  CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
  bad_c = coordinate();
  bad_c.isolation_upper_mz = bad_c.isolation_lower_mz - 1;
  CHECK(invalid([&] { RawTile({}, {1, 1}, {{bad_c, 0}}); }))
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
