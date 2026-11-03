// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIASparseEvidence.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenDIA
{
  namespace
  {
    [[noreturn]] void invalid(const char* message)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "OpenDIA sparse evidence");
    }

    bool contains(const QueryInterval& interval, double value)
    {
      const bool lower_closed = interval.edges == EdgePolicy::CLOSED || interval.edges == EdgePolicy::LOWER_CLOSED;
      const bool upper_closed = interval.edges == EdgePolicy::CLOSED || interval.edges == EdgePolicy::UPPER_CLOSED;
      return (value > interval.lower || (lower_closed && value == interval.lower)) &&
             (value < interval.upper || (upper_closed && value == interval.upper));
    }

    bool overlapsIsolation(const QueryInterval& query, double lower, double upper)
    {
      // Isolation geometry is reported as physical closed bounds by raw readers.
      const bool query_lower_closed = query.edges == EdgePolicy::CLOSED || query.edges == EdgePolicy::LOWER_CLOSED;
      const bool query_upper_closed = query.edges == EdgePolicy::CLOSED || query.edges == EdgePolicy::UPPER_CLOSED;
      const bool left_overlap = upper > query.lower || (query_lower_closed && upper == query.lower);
      const bool right_overlap = lower < query.upper || (query_upper_closed && lower == query.upper);
      return left_overlap && right_overlap;
    }
  }

  SparseEvidenceExtractor::SparseEvidenceExtractor(RecordContext context, RecordLimits limits,
                                                   std::vector<ExtractionQuery> selectors,
                                                   std::vector<FragmentQuery> queries)
    : context_(context), limits_(limits)
  {
    if (context_.schema_version != 0) { invalid("Unsupported sparse evidence context"); }
    if (selectors.size() > limits_.records || queries.size() > limits_.records) { invalid("Fragment index cap exceeded"); }
    selectors_.reserve(selectors.size());
    for (const auto& selector : selectors)
    {
      auto specification = selector.specification();
      if (specification.context != context_) { invalid("Fragment selector context differs from index context"); }
      selectors_.push_back(std::move(specification));
    }
    queries_.reserve(queries.size());
    for (auto& query : queries)
    {
      if (query.selector_index >= selectors_.size()) { invalid("Fragment query references an unknown selector"); }
      if (!std::isfinite(query.mz.lower) || !std::isfinite(query.mz.upper) || query.mz.lower > query.mz.upper ||
          (query.mz.edges != EdgePolicy::OPEN && query.mz.edges != EdgePolicy::CLOSED &&
           query.mz.edges != EdgePolicy::LOWER_CLOSED && query.mz.edges != EdgePolicy::UPPER_CLOSED))
      {
        invalid("Invalid fragment m/z interval");
      }
      queries_.push_back({query.family_id, query.precursor_id, query.fragment_id, query.selector_index, query.mz});
    }

    order_.resize(queries_.size());
    for (OpenMS::Size i = 0; i < order_.size(); ++i) { order_[i] = i; }
    std::sort(order_.begin(), order_.end(), [&](OpenMS::Size lhs, OpenMS::Size rhs)
    {
      const auto& left = queries_[lhs].mz;
      const auto& right = queries_[rhs].mz;
      if (left.lower != right.lower) { return left.lower < right.lower; }
      if (left.upper != right.upper) { return left.upper < right.upper; }
      return lhs < rhs;
    });
    if (!order_.empty())
    {
      if (order_.size() > std::numeric_limits<OpenMS::Size>::max() / 4) { invalid("Fragment interval index size overflow"); }
      max_upper_.resize(order_.size() * 4);
      buildMaxUpper_(1, 0, order_.size());
    }
    index_bytes_ = checkedByteAdd(checkedByteMultiply(queries_.capacity(), sizeof(IndexedQuery)),
                                  checkedByteMultiply(selectors_.capacity(), sizeof(ExtractionQuerySpec)));
    index_bytes_ = checkedByteAdd(index_bytes_, checkedByteMultiply(order_.capacity(), sizeof(OpenMS::Size)));
    index_bytes_ = checkedByteAdd(index_bytes_, checkedByteMultiply(max_upper_.capacity(), sizeof(double)));
    for (const auto& selector : selectors_)
    {
      if (selector.acquisitions)
      {
        index_bytes_ = checkedByteAdd(index_bytes_,
          checkedByteMultiply(selector.acquisitions->capacity(), sizeof(std::pair<RecordId, RecordId>)));
      }
    }
  }

  void SparseEvidenceExtractor::buildMaxUpper_(OpenMS::Size node, OpenMS::Size begin, OpenMS::Size end)
  {
    if (end - begin == 1)
    {
      max_upper_[node] = queries_[order_[begin]].mz.upper;
      return;
    }
    const auto middle = begin + (end - begin) / 2;
    buildMaxUpper_(node * 2, begin, middle);
    buildMaxUpper_(node * 2 + 1, middle, end);
    max_upper_[node] = std::max(max_upper_[node * 2], max_upper_[node * 2 + 1]);
  }

  void SparseEvidenceExtractor::findMatches_(OpenMS::Size node, OpenMS::Size begin, OpenMS::Size end, double mz,
                                             std::vector<OpenMS::Size>& matches, SparseEvidenceStats& stats) const
  {
    ++stats.interval_nodes_visited;
    if (queries_[order_[begin]].mz.lower > mz || max_upper_[node] < mz) { return; }
    if (end - begin == 1)
    {
      const auto query_index = order_[begin];
      if (contains(queries_[query_index].mz, mz)) { matches.push_back(query_index); }
      return;
    }
    const auto middle = begin + (end - begin) / 2;
    findMatches_(node * 2, begin, middle, mz, matches, stats);
    findMatches_(node * 2 + 1, middle, end, mz, matches, stats);
  }

  bool SparseEvidenceExtractor::selectorsMatch_(const ExtractionQuerySpec& query, const Coordinates& coordinates) const
  {
    if (query.ms_level != 0 && query.ms_level != coordinates.ms_level) { return false; }
    if (query.rt_seconds && !contains(*query.rt_seconds, coordinates.rt_seconds)) { return false; }
    if (query.isolation_mz && !overlapsIsolation(*query.isolation_mz, coordinates.isolation_lower_mz, coordinates.isolation_upper_mz))
    {
      return false;
    }
    if (query.mobility)
    {
      if (!coordinates.mobility || coordinates.mobility->unit != query.mobility->unit ||
          !contains(query.mobility->interval, coordinates.mobility->value))
      {
        return false;
      }
    }
    if (query.faims_voltage && (!coordinates.faims_voltage || *query.faims_voltage != *coordinates.faims_voltage)) { return false; }
    if (query.acquisitions && !std::binary_search(query.acquisitions->begin(), query.acquisitions->end(),
                                                  std::pair(coordinates.frame_id, coordinates.scan_id)))
    {
      return false;
    }
    return true;
  }

  std::vector<SparseEvidenceRef> SparseEvidenceExtractor::extract(const RawTile& tile, SparseEvidenceStats* counters) const
  {
    if (tile.context() != context_) { invalid("Raw tile context differs from fragment index context"); }
    SparseEvidenceStats stats;
    std::vector<SparseEvidenceRef> output;
    std::vector<OpenMS::Size> matches;
    for (OpenMS::Size point_index = 0; point_index < tile.size(); ++point_index)
    {
      ++stats.raw_points_examined;
      const auto point = tile.point(point_index);
      matches.clear();
      if (!order_.empty()) { findMatches_(1, 0, order_.size(), point.coordinates.mz, matches, stats); }
      for (const auto query_index : matches)
      {
        const auto& query = queries_[query_index];
        if (!selectorsMatch_(selectors_[query.selector_index], point.coordinates)) { continue; }
        if (output.size() >= limits_.values) { invalid("Sparse evidence output cap exceeded"); }
        output.push_back({query.family_id, query.precursor_id, query.fragment_id, point_index});
      }
    }
    stats.queries_matched = output.size();
    stats.index_bytes = index_bytes_;
    stats.output_bytes = checkedByteMultiply(output.capacity(), sizeof(SparseEvidenceRef));
    stats.retained_bytes = checkedByteAdd(stats.index_bytes, stats.output_bytes);
    if (counters) { *counters = stats; }
    return output;
  }

  std::vector<SparseEvidenceRef> SparseEvidenceExtractor::extract(const RawRegion::Scan& scan, SparseEvidenceStats* counters) const
  {
    std::vector<SparseEvidenceRef> output;
    if (!tryExtract(scan, output, counters)) { invalid("Sparse evidence output cap exceeded"); }
    return output;
  }

  bool SparseEvidenceExtractor::tryExtract(const RawRegion::Scan& scan, std::vector<SparseEvidenceRef>& output,
                                            SparseEvidenceStats* counters) const
  {
    const auto acquisition = scan.acquisition();
    if (scan.context() != context_) { invalid("Raw scan context differs from fragment index context"); }
    const auto mz = scan.mz();
    const auto mobility = scan.mobility();
    const bool selectors_need_mobility = std::any_of(selectors_.begin(), selectors_.end(), [](const auto& selector)
    {
      return selector.mobility.has_value();
    });
    if (!mz.empty() && selectors_need_mobility && (! acquisition.point_mobility || mobility.empty()))
    {
      invalid("Requested point mobility is absent");
    }
    for (const auto& selector : selectors_)
    {
      if (selector.isolation_mz || selector.acquisitions)
      {
        invalid("Native scan view does not provide the requested isolation/acquisition selector");
      }
    }

    SparseEvidenceStats stats;
    output.clear();
    std::vector<OpenMS::Size> matches;
    for (OpenMS::Size point_index = 0; point_index < mz.size(); ++point_index)
    {
      ++stats.raw_points_examined;
      matches.clear();
      if (!order_.empty()) { findMatches_(1, 0, order_.size(), mz[point_index], matches, stats); }
      for (const auto query_index : matches)
      {
        const auto& query = queries_[query_index];
        const auto& selector = selectors_[query.selector_index];
        if (selector.ms_level != 0 && selector.ms_level != acquisition.ms_level) { continue; }
        if (selector.rt_seconds && !contains(*selector.rt_seconds, acquisition.rt_seconds)) { continue; }
        if (selector.faims_voltage && (!acquisition.faims_voltage || *selector.faims_voltage != *acquisition.faims_voltage)) { continue; }
        if (selector.mobility &&
            !contains(selector.mobility->interval, mobility[point_index]))
        {
          continue;
        }
        if (output.size() >= limits_.values)
        {
          output.clear();
          return false;
        }
        output.push_back({query.family_id, query.precursor_id, query.fragment_id, point_index});
      }
    }
    stats.queries_matched = output.size();
    stats.index_bytes = index_bytes_;
    stats.output_bytes = checkedByteMultiply(output.capacity(), sizeof(SparseEvidenceRef));
    stats.retained_bytes = checkedByteAdd(stats.index_bytes, stats.output_bytes);
    if (counters) { *counters = stats; }
    return true;
  }

  OpenMS::Size SparseEvidenceExtractor::queryCount() const noexcept { return queries_.size(); }
}
