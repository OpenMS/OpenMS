/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <arrow/record_batch.h>
#include <parquet/metadata.h>
#include <parquet/statistics.h>

#include "mzpeak/io/file.h"
#include "mzpeak/schema/file.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/util/executor.h"
#include "mzpeak/util/query.h"
#include "mzpeak/util/row_group_cache.h"

namespace MzPeak::Util {

/// Diagnostic counters for the per-entity read path.  Free-running and
/// relaxed; they exist to answer "how much work does one spectrum cost"
/// without a profiler, and cost an untaken increment when nobody reads them.
struct ReadCounters {
  long plan_group_evals = 0; ///< row groups examined while planning
  long batches_visited = 0;  ///< record batches walked in the executor
  long slices_made = 0;      ///< RecordBatch::Slice calls
  long plan_pruned = 0;      ///< groups rejected by column statistics
  long plan_full_scan = 0;   ///< groups planned as a WHOLE-GROUP range
  long plan_page_index = 0;  ///< groups narrowed by the page index
  long plan_pi_null = 0;     ///< groups with no usable page-index reader
  long plan_ranges = 0;      ///< ranges emitted
  long plan_range_rows = 0;  ///< rows covered by those ranges
  long ns_plan = 0;          ///< nanoseconds building + running the planner
  long ns_plan_ctor = 0;     ///< of which: constructing the Planner
  long ns_exec = 0;          ///< nanoseconds in Executor::execute
  long ns_rowgroup = 0;      ///< nanoseconds obtaining a decoded row group
  long ns_project = 0;       ///< nanoseconds in Executor::Impl::project
};
ReadCounters read_counters();
void count_plan_group_eval();
void count_plan_pruned();
void count_plan_full_scan();
void count_plan_page_index();
void count_plan_pi_null();
void count_plan_range(long rows);
void count_ns_plan(long ns);
void count_ns_plan_ctor(long ns);
void count_ns_exec(long ns);
void count_ns_rowgroup(long ns);
void count_ns_project(long ns);

/// A row group's DECODED footprint: what an Arrow RecordBatch of it costs.
///
/// NOT Parquet's total_byte_size, which is the ENCODED size and can be several
/// times smaller -- a sorted, repetitive column is RLE'd to almost nothing on
/// disk and expands straight back out in memory.  Anyone sizing a memory
/// budget in row groups must use this, or the budget is in the wrong unit.
std::size_t decoded_row_group_bytes(const parquet::RowGroupMetaData& group,
                                    const parquet::SchemaDescriptor& schema);

/**
 * Low-level wrapper around Parquet files.
 */
class Parquet final {
public:
  using file_metadata_t = std::shared_ptr<parquet::FileMetaData>;

  /// Constructor.
  /// @param cache  shared with every other Parquet over the same archive, so a
  ///   row group is decoded once for all of them; without it this object keeps
  ///   a private two-group cache.
  /// @param footer  this file's footer, already parsed by another Parquet over
  ///   the same member (see Manager::parquet()), so it is not read and parsed
  ///   again.  Null reads it from the file.  Read-only once parsed, so any
  ///   number of readers may share one, on any threads.
  Parquet(std::unique_ptr<MzPeak::IO::File>,
          Schema::File,
          std::shared_ptr<RowGroupCache> cache = nullptr,
          file_metadata_t footer = nullptr);

  /// Destructor.
  ~Parquet();

  /**
   * Return the file information from the MzPeak index.
   */
  const Schema::File& index_file() const;

  /**
   * Return the schema encoded as a map of Group objects.
   */
  const std::shared_ptr<Schema::GroupMap>& groups() const;

  /**
   * Return a Group and Field matching the given names.
   */
  std::optional<Schema::Column> field(std::string_view, std::string_view) const;

  /**
   * Access the file metadata.
   */
  file_metadata_t file_metadata() const;

  /**
   * Fetch a string value from the metadata key-value store.
   */
  std::optional<std::string> kv_string(const file_metadata_t&,
                                       std::string_view) const;

  /**
   * Fetch a std::size_t value from the metadata key-value store.
   */
  std::optional<std::size_t> kv_size_t(const file_metadata_t&,
                                       std::string_view) const;

  /**
   * Directly access the FileReader.  This reference is only valid
   * while this Parquet object exists.
   */
  parquet::arrow::FileReader& reader() const;

  /**
   * The record batches of one row group, decoded once and retained.
   *
   * A query that selects a single entity decodes the row group holding it and
   * slices that down to a few thousand rows.  Reading a run entity by entity
   * therefore decoded the same group hundreds of times -- on a 21-row-group
   * file, 13,009 decodes where 21 would do.  This hands back the same decoded
   * batches to every caller that lands in the group.
   *
   * Only the most recent few groups are kept -- two without a shared
   * RowGroupCache, what its budget allows with one: a decoded group is tens
   * of megabytes, and a forward pass never looks back.
   *
   * Returned by shared_ptr, not by reference: the cache evicts, and a reference
   * into it would dangle the moment a caller fetched a second group while still
   * holding the first.
   *
   * @note An arrow::Array sliced out of these batches SHARES their buffers, so
   * anything retaining such a slice keeps the whole decoded row group alive.
   * Spectrum copies its values out and drops the slice; Chromatogram keeps its
   * decoder, and so keeps a group resident for as long as the caller holds it.
   * That is bounded by how many entities the caller holds, not by the run, and
   * comes on top of a shared cache's budget.
   *
   * @note Reading through a shared RowGroupCache, this object also HOLDS the
   * last group it returned -- a memo, so asking for that group again skips
   * the cache -- until it returns another group or is destroyed, idle or not.
   * The cache neither evicts a held group nor counts it as room, so every
   * live Parquet reading through it keeps about one group on top of the
   * cache's budget: readers on distinct groups hold that many groups,
   * whatever the budget.  See RowGroupCache.
   *
   * @note A group whose declared-sorted key proves sorted is decoded WITHOUT
   * that column: RowGroupBatches::key_runs holds it as runs and a null
   * placeholder takes its place in the batches (key_runs.h).  The Executor
   * handles this; a direct caller must check key_runs before reading the key.
   */
  using RowGroupBatches = Util::RowGroupBatches;
  std::shared_ptr<const RowGroupBatches> row_group(int32_t);

  /**
   * Return a planner for the given query.
   */
  /// See StatsIndex::sorted_ascending.
  bool sorted_ascending(int32_t row_group, int32_t leaf_column) const;

  Planner planner(const Query&);

  /**
   * Return an executor that will capture the requested fields.
   */
  Executor executor(const Projection&);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::Util
