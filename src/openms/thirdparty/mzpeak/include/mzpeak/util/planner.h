/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <cstdint>
#include <memory>

#include "mzpeak/util/query.h"

namespace parquet {
class PageIndexReader;
}

namespace parquet {
class FileMetaData;
class Statistics;
} // namespace parquet

namespace parquet::arrow {
class FileReader;
}

namespace MzPeak::Util {

class Parquet;

/**
 * File-lifetime cache of per-row-group column statistics.
 *
 * The planner consults the min/max of a column in EVERY row group on EVERY
 * query, and materialising those from the file metadata allocates a
 * `RowGroupMetaData`, a `ColumnChunkMetaData` and a `Statistics` every time.
 * That is a few microseconds per row group: nothing on a 21-group file, and
 * about 18 ms per spectrum on a 3,518-group one -- a Bruker diaPASEF run --
 * before a single peak has been decoded.
 *
 * None of it changes while the file is open, so it is read once.  Entries fill
 * on demand because a run has many columns and a query touches one or two.
 */
class StatsIndex final {
public:
  /// Constructor.
  explicit StatsIndex(std::shared_ptr<parquet::FileMetaData>);

  /// Destructor.
  ~StatsIndex();

  /// Number of row groups in the file.
  int32_t row_group_count() const;

  /// Rows in the given row group.
  int64_t row_count(int32_t row_group) const;

  /**
   * Statistics for a (row group, leaf column) pair, or null when that column
   * chunk carries none.  Absence is cached too: a column without statistics
   * must not be re-examined once per query for the life of the file.
   */
  std::shared_ptr<parquet::Statistics> get(int32_t row_group, int32_t column) const;

  /// Does the row group DECLARE @p column sorted ascending, nulls last?  Only a
  /// file's own sorting-columns declaration counts -- a binary search over a
  /// column that is not actually sorted silently misses rows.
  bool sorted_ascending(int32_t row_group, int32_t column) const;

  /// Is @p column declared sorted ascending in EVERY row group?  Only then are
  /// the row groups matching an equality contiguous.
  bool all_sorted_ascending(int32_t column) const;

  /// The row group that last matched through this reader, and a setter for it.
  /// A hint only: it can change which groups are examined first, never which
  /// ones end up in the plan.
  int32_t group_hint() const;
  void set_group_hint(int32_t) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/**
 * Create a plan for the most efficient way to execute a query against
 * a Parquet file.
 */
class Planner final {
public:
  /// How to limit table reads.
  struct Range {
    int32_t row_group; /// The row group to request.
    int64_t offset;    /// The first row to read.
    int64_t length;    /// The number of rows to read.
  };

  /// Information about how to carry out a query.
  struct Plan {
    Query query;
    std::vector<Range> ranges;
  };

  /// Destructor.
  ~Planner();

  /**
   * Construct a query plan.
   */
  const Plan& plan() const;

private:
  friend class Parquet;

  /// Constructor.
  ///
  /// `stats` is shared with every other planner over the same file; passing
  /// null falls back to reading statistics from the file metadata each time.
  /// @param page_index  the file's page-index reader, built ONCE by the
  ///   caller. Obtaining it per query re-reads and re-parses index buffers;
  ///   it depends only on the file, so it is built once and shared.
  Planner(parquet::arrow::FileReader&,
          const Query&,
          std::shared_ptr<StatsIndex> stats,
          std::shared_ptr<parquet::PageIndexReader> page_index = nullptr);

  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::Util
