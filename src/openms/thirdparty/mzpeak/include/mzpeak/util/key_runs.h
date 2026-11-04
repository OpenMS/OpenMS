/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "mzpeak/util/types.h"

namespace arrow {
class DataType;
class Field;
class RecordBatch;
} // namespace arrow

namespace parquet {
class RowGroupMetaData;
class SchemaDescriptor;
} // namespace parquet

namespace parquet::arrow {
class FileReader;
}

namespace MzPeak::Util {

struct RowGroupBatches;

/**
 * A row group's declared-sorted key column, held as RUNS instead of as a column.
 *
 * In the point layout every peak row carries its spectrum index: 8 bytes of a
 * ~21-byte row, ~40% of a decoded group, for a column that says nothing but
 * "rows a..b belong to spectrum s".  When the column really is sorted, that is
 * all it can say, so a decoded group keeps one (value, first row) pair per
 * spectrum -- ~2,100 of them for a 1,048,576-row group -- and the key's place in
 * every batch is taken by a zero-byte null placeholder, so every other column
 * keeps its position.
 *
 * Why the runs are built FROM THE COLUMN, not from the metadata's per-spectrum
 * point counts: counts and footer statistics cannot prove where one spectrum
 * ends and the next begins.  A file whose counts for two neighbouring spectra
 * are swapped still sums to the table's row count and still agrees with every
 * row group's min/max, yet puts the boundary in the wrong place -- and then
 * the decode's own count check, the one thing that catches a mislocated
 * spectrum, would compare the ranges against the very counts they came from.
 * So the column is read once per decode, in 64 K-value chunks through the
 * low-level reader (never as an 8 MB array), and checked as it is read.
 *
 * A group is held this way only when all of these hold, and is otherwise
 * decoded in full exactly as before:
 *   - the footer declares the column the first sorting column, ascending,
 *     nulls last (run_key_leaf());
 *   - it is a non-repeated INT64 leaf that decodes to Arrow int64 or uint64,
 *     one or two levels deep, and not the table's only leaf;
 *   - it has no null: for a nullable key, the footer's null count says so,
 *     or, when the footer gives no null count Arrow accepts (no statistics,
 *     statistics is_stats_set() rejects or that do not decode, or no null
 *     count in them), every row read has a definition level that gives it a
 *     value.  Reading the values alone would not do: a dictionary page may
 *     zero-pad its last indices to eight, and the padding reads back as
 *     values.  A footer that counts no null and holds one is believed, as
 *     run_key_leaf() believes it;
 *   - its values never DECREASE, in the order of that Arrow type -- which is
 *     what makes each value's rows one contiguous run, and exactly the run a
 *     binary search over the column finds.  A file that declares sorted and is
 *     not keeps its column, so its behaviour does not change either.
 */
struct KeyRuns {
  int32_t leaf = -1;         ///< the key's Parquet leaf index
  int top = -1;              ///< its top-level column in a batch
  int child = -1;            ///< its child in that struct column; -1 if top-level
  bool is_unsigned = false;  ///< uint64 (else int64): the order the runs follow
  std::shared_ptr<arrow::Field> top_field; ///< the ORIGINAL top-level field
  std::shared_ptr<arrow::DataType> type;   ///< the key's Arrow type

  /// One entry per run, in row order: the value (its 64-bit pattern)...
  std::vector<uint64_t> key;
  /// ...and the run's first row, relative to the group.  One entry longer than
  /// `key`: the last is the group's row count, so run i is [start[i], start[i+1]).
  std::vector<int64_t> start;

  /// Rows [first, last) of the group holding @p value -- empty when no row
  /// does -- or nothing when @p value is not of the column's C++ type.
  std::optional<std::pair<int64_t, int64_t>> rows(const any_value_type& value) const;
};

/// The leaf a row group will hold as runs, judged from the footer alone: the
/// decode's own shape test, plus the footer's word that a nullable key has no
/// nulls.  Used to size a decoded group, so the cache charges what it will
/// actually hold.  Wrong only for a file that lies about sorting (under-count
/// by the key) or has no statistics (over-count).
std::optional<int32_t> run_key_leaf(const parquet::RowGroupMetaData&,
                                    const parquet::SchemaDescriptor&);

/// Read row group @p group's key column into runs, or null when the group must
/// keep its column (see KeyRuns for when).  The caller holds whatever lock
/// guards @p reader's file position.
std::shared_ptr<KeyRuns> scan_key_runs(parquet::arrow::FileReader& reader,
                                       int32_t group);

/// @p batch, read WITHOUT the key, with the placeholder put back in its place;
/// null when the batch does not have the shape the runs describe.
std::shared_ptr<arrow::RecordBatch>
insert_key_placeholder(const std::shared_ptr<arrow::RecordBatch>& batch,
                       const KeyRuns& runs);

/// The group as a full decode would have produced it: the key column rebuilt
/// from @p stripped's runs, and its per-batch key spans.  For the rare query
/// the runs cannot answer directly -- anything but an equality on the key, or
/// a projection that asks for the key itself.
std::shared_ptr<const RowGroupBatches> restore_key_column(const RowGroupBatches& stripped);

} // namespace MzPeak::Util
