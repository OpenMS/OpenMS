// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARecords.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
// Runtime checks deliberately survive NDEBUG, using the same production implementation.
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif

using namespace OpenDIA;
namespace
{
  struct ConstantHash
  {
    std::size_t operator()(const ExtractionQuery&) const { return 0; }
  };
  ExtractionQuerySpec baseline()
  {
    ExtractionQuerySpec s;
    s.context = {0, 101, 202, 303};
    s.acquisitions = {{{1, 4}, {9, 7}}};
    s.ms_level = 2;
    s.isolation_mz = QueryInterval{490, 510, EdgePolicy::LOWER_CLOSED};
    s.rt_seconds = QueryInterval{-1, 20, EdgePolicy::CLOSED};
    s.mz = {499, 501, EdgePolicy::OPEN};
    s.mobility = MobilityInterval{{0, 2, EdgePolicy::OPEN}, MobilityUnit::MILLISECOND};
    s.faims_voltage = 0;
    s.tolerance = {2, ToleranceUnit::TH};
    s.isotope = {IsotopeMode::SELECTED, 1, 2};
    return s;
  }
  // Every mutation below changes exactly one nested scientific field from baseline.
  using Mutation = std::function<void(ExtractionQuerySpec&)>;
  std::vector<Mutation> mutations()
  {
    return {
      [](auto& s) { ++s.context.run_id; },
      [](auto& s) { ++s.context.lane_id; },
      [](auto& s) { ++s.context.calibration_snapshot_id; },
      [](auto& s) { s.acquisitions.reset(); },
      [](auto& s) { s.acquisitions->clear(); },
      [](auto& s) { ++s.acquisitions->front().first; },
      [](auto& s) { ++s.acquisitions->front().second; },
      [](auto& s) { ++s.ms_level; },
      [](auto& s) { s.isolation_mz.reset(); },
      [](auto& s) { s.isolation_mz->lower += 1; },
      [](auto& s) { s.isolation_mz->upper += 1; },
      [](auto& s) { s.isolation_mz->edges = EdgePolicy::UPPER_CLOSED; },
      [](auto& s) { s.rt_seconds.reset(); },
      [](auto& s) { s.rt_seconds->lower += 1; },
      [](auto& s) { s.rt_seconds->upper += 1; },
      [](auto& s) { s.rt_seconds->edges = EdgePolicy::OPEN; },
      [](auto& s) { s.mz.lower = std::nextafter(s.mz.lower, s.mz.upper); },
      [](auto& s) { s.mz.upper = std::nextafter(s.mz.upper, s.mz.lower); },
      [](auto& s) { s.mz.edges = EdgePolicy::CLOSED; },
      [](auto& s) { s.mobility.reset(); },
      [](auto& s) { s.mobility->interval.lower = 0.1; },
      [](auto& s) { s.mobility->interval.upper = 3; },
      [](auto& s) { s.mobility->interval.edges = EdgePolicy::LOWER_CLOSED; },
      [](auto& s) { s.mobility->unit = MobilityUnit::INVERSE_REDUCED; },
      [](auto& s) { s.faims_voltage.reset(); },
      [](auto& s) { s.faims_voltage = -45; },
      [](auto& s) { s.tolerance.full_width = std::nextafter(2.0, 3.0); },
      [](auto& s) { s.tolerance.unit = ToleranceUnit::PPM; },
      [](auto& s) { s.aggregation = QueryAggregation::SUM; },
      [](auto& s) { ++s.isotope.index; },
      [](auto& s) { ++s.isotope.charge; },
      [](auto& s) { s.isotope = {}; }
    };
  }
  bool rejected(ExtractionQuerySpec s)
  {
    try { ExtractionQuery query(std::move(s)); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }
}

#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIAQuery, "$Id$")
START_SECTION(v0 exact query identity)
#else
int main()
{
#endif
  const ExtractionQuery key(baseline());
  std::unordered_map<ExtractionQuery, unsigned, ConstantHash> values;
  CHECK(values.emplace(key, 0).second)
  std::vector<ExtractionQuery> keys{key};
  unsigned payload = 0;
  for (const auto& mutate : mutations())
  {
    auto s = baseline();
    mutate(s);
    ExtractionQuery changed(s);
    CHECK(changed != key)
    CHECK(!values.contains(changed))
    CHECK(values.emplace(changed, ++payload).second)
    keys.push_back(changed);
  }
  for (unsigned i = 0; i < keys.size(); ++i)
  {
    CHECK(values.at(keys[i]) == i)
  }
  CHECK(values.size() == mutations().size() + 1)
  OPENMS_LOG_INFO << "Retained " << values.size() << " distinct query payloads under one constant hash after "
                  << mutations().size() << " independent field mutations\n";
  auto equivalent = baseline();
  equivalent.acquisitions = {{{9, 7}, {1, 4}, {9, 7}}};
  equivalent.faims_voltage = -0.0;
  equivalent.mobility->interval.lower = -0.0;
  CHECK(ExtractionQuery(equivalent) == key)
  CHECK(!values.emplace(ExtractionQuery(equivalent), 999).second)
  CHECK(values.at(key) == 0)
  CHECK(key.specification().acquisitions == baseline().acquisitions)
  auto detached = key.specification();
  detached.context.calibration_snapshot_id++;
  detached.acquisitions->clear();
  detached.mz.lower = 0;
  CHECK(key == ExtractionQuery(baseline()))
  auto input = baseline();
  ExtractionQuery owned(input);
  input.acquisitions->clear();
  CHECK(owned == key)
  auto none = baseline();
  none.isotope = {};
  const ExtractionQuery no_isotope(none);
  none.isotope = {IsotopeMode::SELECTED, 0, 1};
  CHECK(ExtractionQuery(none) != no_isotope)
  // A mode-only change is invalid when inactive fields would become ambiguous.
  none.isotope.mode = IsotopeMode::NONE;
  CHECK(rejected(none))
  none = baseline();
  none.tolerance.full_width = 0;
  none.mz.upper = none.mz.lower;
  const ExtractionQuery zero_width(none);
  none.tolerance.full_width = -0.0;
  CHECK(ExtractionQuery(none) == zero_width)

  const std::vector<Mutation> malformed{
    [](auto& s) { s.context.schema_version = 1; },
    [](auto& s) { s.ms_level = 0; },
    [](auto& s) { s.mobility->unit = static_cast<MobilityUnit>(99); },
    [](auto& s) { s.tolerance.unit = static_cast<ToleranceUnit>(99); },
    [](auto& s) { s.aggregation = static_cast<QueryAggregation>(99); },
    [](auto& s) { s.tolerance.full_width = -1; },
    [](auto& s) { s.isotope.mode = static_cast<IsotopeMode>(99); },
    [](auto& s) { s.isotope.charge = 0; },
    [](auto& s) { s.isotope = {IsotopeMode::NONE, 1, 0}; },
    [](auto& s) { s.isotope = {IsotopeMode::NONE, 0, 1}; }
  };
  for (const auto& mutate : malformed)
  {
    auto s = baseline();
    mutate(s);
    CHECK(rejected(s))
  }
  // Exercise every interval's bounds, ordering and enum validation independently.
  const std::vector<std::function<QueryInterval&(ExtractionQuerySpec&)>> intervals{
    [](auto& s) -> QueryInterval& { return s.mz; },
    [](auto& s) -> QueryInterval& { return *s.isolation_mz; },
    [](auto& s) -> QueryInterval& { return *s.rt_seconds; },
    [](auto& s) -> QueryInterval& { return s.mobility->interval; }
  };
  for (const auto& axis : intervals)
  {
    auto s = baseline();
    axis(s).lower = axis(s).upper + 1;
    CHECK(rejected(s))
    s = baseline();
    axis(s).edges = static_cast<EdgePolicy>(99);
    CHECK(rejected(s))
    for (const auto edge : {EdgePolicy::OPEN, EdgePolicy::CLOSED, EdgePolicy::LOWER_CLOSED, EdgePolicy::UPPER_CLOSED})
    {
      s = baseline();
      axis(s).edges = edge;
      CHECK(ExtractionQuery(s).specification() == s)
    }
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
      for (auto field : {&QueryInterval::lower, &QueryInterval::upper})
      {
        s = baseline();
        axis(s).*field = bad;
        CHECK(rejected(s))
      }
    }
  }
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
  {
    auto s = baseline();
    s.faims_voltage = bad;
    CHECK(rejected(s))
    s = baseline();
    s.tolerance.full_width = bad;
    CHECK(rejected(s))
  }
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
