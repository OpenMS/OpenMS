// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIASparseEvidence.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <functional>
#include <stdexcept>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif

using namespace OpenDIA;
namespace
{
  RecordContext context() { return {0, 4, 5, 6}; }

  ExtractionQuery makeQuery(QueryInterval mz, std::optional<QueryInterval> rt = std::nullopt,
                            std::optional<QueryInterval> isolation = std::nullopt,
                            std::optional<MobilityInterval> mobility = std::nullopt,
                            std::optional<double> faims = std::nullopt,
                            std::optional<std::vector<std::pair<RecordId, RecordId>>> acquisitions = std::nullopt,
                            unsigned ms_level = 2)
  {
    ExtractionQuerySpec spec;
    spec.context = context();
    spec.acquisitions = std::move(acquisitions);
    spec.ms_level = ms_level;
    spec.isolation_mz = isolation;
    spec.rt_seconds = rt;
    spec.mz = mz;
    spec.mobility = mobility;
    spec.faims_voltage = faims;
    spec.tolerance = {mz.upper - mz.lower, ToleranceUnit::TH};
    return ExtractionQuery(std::move(spec));
  }

  RawPoint point(double mz, double rt = 15, double isolation_lower = 400, double isolation_upper = 410,
                 std::optional<Mobility> mobility = Mobility{1.0, MobilityUnit::INVERSE_REDUCED},
                 std::optional<double> faims = -45, RecordId frame = 1, RecordId scan = 2, unsigned ms_level = 2)
  {
    Coordinates c{rt, mz, frame, scan, ms_level, isolation_lower, isolation_upper, mobility, faims, "scan-id"};
    return {c, 12.5};
  }

  bool rejected(const std::function<void()>& action)
  {
    try { action(); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }
}

#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIASparseEvidence, "$Id$")
START_SECTION(interval-indexed sparse raw evidence)
#else
int main()
{
#endif
  std::vector<ExtractionQuery> selectors;
  selectors.push_back(makeQuery({100, 101, EdgePolicy::OPEN}, QueryInterval{10, 20, EdgePolicy::CLOSED},
                                QueryInterval{400, 410, EdgePolicy::CLOSED},
                                MobilityInterval{{0.5, 1.5, EdgePolicy::CLOSED}, MobilityUnit::INVERSE_REDUCED},
                                -45.0, std::vector<std::pair<RecordId, RecordId>>{{1, 2}}));
  selectors.push_back(makeQuery({100.5, 101.5, EdgePolicy::CLOSED}, QueryInterval{10, 20, EdgePolicy::CLOSED},
                                QueryInterval{405, 415, EdgePolicy::CLOSED},
                                MobilityInterval{{0.5, 1.5, EdgePolicy::CLOSED}, MobilityUnit::INVERSE_REDUCED},
                                -45.0, std::vector<std::pair<RecordId, RecordId>>{{1, 2}}));
  std::vector<FragmentQuery> fragment_queries{
    {10, 20, 30, 0, {100, 101, EdgePolicy::OPEN}},
    {11, 21, 31, 1, {100.5, 101.5, EdgePolicy::CLOSED}}
  };
  SparseEvidenceExtractor extractor(context(), {2, 20}, std::move(selectors), std::move(fragment_queries));
  CHECK(extractor.queryCount() == 2)

  std::vector<RawPoint> points{
    point(100.5), point(101.0, 20), point(100.7, 21), point(100.7, 15, 400, 410, std::nullopt),
    point(100.7, 15, 400, 410, Mobility{1.0, MobilityUnit::MILLISECOND}), point(100.7, 15, 400, 410, Mobility{1.0, MobilityUnit::INVERSE_REDUCED}, 0.0),
    point(100.7, 15, 400, 410, Mobility{1.0, MobilityUnit::INVERSE_REDUCED}, -45.0, 1, 3),
    point(100.7, 15, 420, 430), point(100.0), point(101.5)
  };
  const auto point_count = points.size();
  RawTile tile(context(), {point_count, point_count}, std::move(points));
  SparseEvidenceStats stats;
  const auto evidence = extractor.extract(tile, &stats);
  CHECK(evidence.size() == 4)
  CHECK(evidence[0] == SparseEvidenceRef{10, 20, 30, 0})
  CHECK(evidence[1] == SparseEvidenceRef{11, 21, 31, 0})
  CHECK(evidence[2] == SparseEvidenceRef{11, 21, 31, 1})
  CHECK(evidence[3] == SparseEvidenceRef{11, 21, 31, 9})
  CHECK(tile.point(evidence[0].raw_point_index).intensity == 12.5)
  CHECK(tile.point(evidence[0].raw_point_index).coordinates.native_id == "scan-id")
  CHECK(stats.raw_points_examined == 10 && stats.queries_matched == evidence.size())
  CHECK(stats.index_bytes > 0 && stats.output_bytes >= evidence.size() * sizeof(SparseEvidenceRef))
  CHECK(stats.retained_bytes == stats.index_bytes + stats.output_bytes)

  OpenMS::MSSpectrum spectrum;
  spectrum.setNativeID("scan=2");
  spectrum.setRT(15);
  spectrum.setMSLevel(2);
  for (const auto mz : {100.5, 101.0})
  {
    OpenMS::Peak1D peak;
    peak.setMZ(mz);
    peak.setIntensity(12.5);
    spectrum.push_back(peak);
  }
  const auto raw_region_bytes = RawRegion::requiredBytes({1, 2}, spectrum.getNativeID().size());
  ByteBudget raw_region_budget(raw_region_bytes, 0);
  auto raw_region = RawRegion::fromSpectrum(spectrum, 1, context(), raw_region_budget);
  SparseEvidenceExtractor zero_copy(context(), {2, 2}, {makeQuery({100, 101, EdgePolicy::OPEN},
                                                                  QueryInterval{10, 20, EdgePolicy::CLOSED})},
                                    {{12, 22, 32, 0, {100, 101, EdgePolicy::OPEN}}});
  SparseEvidenceStats zero_copy_stats;
  const auto native_evidence = zero_copy.extract(raw_region.scan(0), &zero_copy_stats);
  CHECK(native_evidence.size() == 1 && native_evidence.front().raw_point_index == 0)
  CHECK(zero_copy_stats.raw_points_examined == 2 && zero_copy_stats.queries_matched == 1)

  // A thousand distant intervals should not require a thousand checks for one observed point.
  std::vector<FragmentQuery> distant_queries;
  for (OpenMS::Size i = 0; i < 1024; ++i)
  {
    const double lower = static_cast<double>(i) * 2.0;
    distant_queries.push_back({1, i + 1, i + 2, 0, {lower, lower + 1, EdgePolicy::CLOSED}});
  }
  SparseEvidenceExtractor indexed(context(), {1024, 10}, {makeQuery({0, 1})}, std::move(distant_queries));
  RawTile one_point(context(), {1, 1}, {point(1000.5)});
  SparseEvidenceStats index_stats;
  const auto indexed_evidence = indexed.extract(one_point, &index_stats);
  CHECK(indexed_evidence.size() == 1 && indexed_evidence.front().raw_point_index == 0)
  CHECK(index_stats.interval_nodes_visited < 100)

  CHECK(rejected([&]
  {
    RawTile wrong_context({0, 4, 5, 7}, {1, 1}, {point(1000.5)});
    indexed.extract(wrong_context);
  }))
  CHECK(rejected([&]
  {
    SparseEvidenceExtractor too_small(context(), {1024, 0}, {makeQuery({0, 1})}, {{1, 2, 3, 0, {1000, 1001}}});
    too_small.extract(one_point);
  }))
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
