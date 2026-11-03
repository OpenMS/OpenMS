// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include "OpenDIARecords.h"
#include "OpenDIARawRegion.h"
#include <vector>

namespace OpenDIA
{
  /// One fragment interval references a shared extraction selector instead of copying it.
  struct FragmentQuery
  {
    RecordId family_id = 0;
    RecordId precursor_id = 0;
    RecordId fragment_id = 0;
    OpenMS::Size selector_index = 0;
    QueryInterval mz;
  };

  /// A compact match referring to the measured point retained by its RawTile.
  struct SparseEvidenceRef
  {
    RecordId family_id = 0;
    RecordId precursor_id = 0;
    RecordId fragment_id = 0;
    OpenMS::Size raw_point_index = 0;
    bool operator==(const SparseEvidenceRef&) const = default;
  };

  struct SparseEvidenceStats
  {
    OpenMS::Size raw_points_examined = 0;
    OpenMS::Size interval_nodes_visited = 0;
    OpenMS::Size queries_matched = 0;
    /// Vector capacities retained by the index, including selector acquisition lists.
    ByteCount index_bytes = 0;
    /// Sparse match-vector capacity. Neither byte counter includes allocator overhead or the input tile.
    ByteCount output_bytes = 0;
    ByteCount retained_bytes = 0;
  };

  /** @brief Private m/z interval index that joins a raw tile to sparse fragment evidence.

      Queries are sorted and indexed once. Extraction visits only intervals containing each
      measured m/z, then applies the exact acquisition, RT, isolation, mobility and FAIMS
      selectors. Matches reference the RawTile by row index, avoiding coordinate copies.
  */
  class SparseEvidenceExtractor
  {
  public:
    /**
      @brief Build a bounded interval index for one immutable caller context.
      @param[in] context Run/lane/calibration identity shared by every query and input tile.
      @param[in] limits Maximum indexed queries and emitted evidence rows.
      @param[in] selectors Shared extraction constraints, usually one per precursor/acquisition scope.
      @param[in] queries Fragment m/z intervals and stable candidate identities that reference selectors.
      @throws OpenMS::Exception::InvalidValue for mismatched contexts, bad selector references or cap violations.
    */
    SparseEvidenceExtractor(RecordContext context, RecordLimits limits, std::vector<ExtractionQuery> selectors,
                           std::vector<FragmentQuery> queries);

    /**
      @brief Join a raw tile to all matching fragment queries without scanning the full library per point.
      @param[in] tile Immutable measured points with exact acquisition coordinates.
      @param[out] stats Optional deterministic work counters; reset on entry.
      @return Sparse row references into tile. An absent row is not a measured zero.
      @throws OpenMS::Exception::InvalidValue for context mismatch or output-cap exhaustion.
    */
    std::vector<SparseEvidenceRef> extract(const RawTile& tile, SparseEvidenceStats* stats = nullptr) const;

    /**
      @brief Join a native region scan without materializing a RawTile copy.
      @param[in] scan Immutable native point columns and acquisition metadata.
      @param[out] stats Optional deterministic work counters; reset on entry.
      @return Sparse row references into the scan's point arrays.
      @throws OpenMS::Exception::InvalidValue for unsupported selector fields or output-cap exhaustion.
    */
    std::vector<SparseEvidenceRef> extract(const RawRegion::Scan& scan, SparseEvidenceStats* stats = nullptr) const;

    /// Bounded variant for streaming callers. Returns false and clears output when the cap is exceeded.
    bool tryExtract(const RawRegion::Scan& scan, std::vector<SparseEvidenceRef>& output,
                    SparseEvidenceStats* stats = nullptr) const;

    OpenMS::Size queryCount() const noexcept;

  private:
    struct IndexedQuery
    {
      RecordId family_id;
      RecordId precursor_id;
      RecordId fragment_id;
      OpenMS::Size selector_index;
      QueryInterval mz;
    };

    void buildMaxUpper_(OpenMS::Size node, OpenMS::Size begin, OpenMS::Size end);
    void findMatches_(OpenMS::Size node, OpenMS::Size begin, OpenMS::Size end, double mz,
                      std::vector<OpenMS::Size>& matches, SparseEvidenceStats& stats) const;
    bool selectorsMatch_(const ExtractionQuerySpec& query, const Coordinates& coordinates) const;

    RecordContext context_;
    RecordLimits limits_;
    std::vector<IndexedQuery> queries_;
    std::vector<ExtractionQuerySpec> selectors_;
    std::vector<OpenMS::Size> order_;
    std::vector<double> max_upper_;
    ByteCount index_bytes_ = 0;
  };
}
