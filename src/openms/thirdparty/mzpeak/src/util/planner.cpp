/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include <atomic>
#include <deque>
#include "mzpeak/util/planner.h"

#include <algorithm>

#include "mzpeak/util/parquet.h"

#include <arrow/array.h>
#include <arrow/record_batch.h>
#include <bit>
#include <memory>
#include <mutex>
#include <parquet/api/reader.h>
#include <parquet/arrow/reader.h>
#include <parquet/page_index.h>
#include <parquet/statistics.h>
// OpenMS: <print> dropped (unused; not in libstdc++ 13)
#include <ranges>

#include "mzpeak/util/algorithm.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Util {

using namespace std::placeholders;

/******************************************************************************/
/*
 * The statistics cache.
 *
 * Keyed on (row group, leaf column).  The owning `ColumnChunkMetaData` is
 * retained alongside each `Statistics`: for the byte-array types the min/max
 * accessors hand back views into the chunk's thrift buffers, so dropping the
 * chunk would leave them dangling.
 */
struct StatsIndex::Impl {
  using key_type = std::pair<int32_t, int32_t>;
  using value_type = std::pair<std::shared_ptr<parquet::ColumnChunkMetaData>,
                               std::shared_ptr<parquet::Statistics>>;

  /// One lazily-filled, then lock-free, statistics slot.  See get().
  struct Slot {
    std::atomic<bool> ready{false};
    std::shared_ptr<parquet::ColumnChunkMetaData> chunk;
    std::shared_ptr<parquet::Statistics> stats;
  };

  explicit Impl(std::shared_ptr<parquet::FileMetaData> metadata)
      : metadata_(std::move(metadata))
      , rows_(static_cast<std::size_t>(metadata_->num_row_groups()), 0)
      , columns_(metadata_->num_columns())
      , slots_(static_cast<std::size_t>(metadata_->num_row_groups()) *
               static_cast<std::size_t>(metadata_->num_columns()))
  {
    // Row counts are two integers per group, are needed whenever a query falls
    // back to a full scan, and reading them here means the common path never
    // constructs a RowGroupMetaData at all.
    for (int32_t g : std::views::iota(0, metadata_->num_row_groups()))
      rows_[static_cast<std::size_t>(g)] = metadata_->RowGroup(g)->num_rows();
  }

  std::shared_ptr<parquet::Statistics> get(int32_t row_group, int32_t column)
  {
    // A FIXED slot per (row group, column), not a map behind a mutex.
    //
    // A planner is built per entity and evaluates the statistics of EVERY row
    // group, so this is the hottest call in the read path: measured at
    // 276,635,040 calls for one 762,016-spectrum run over a 363-group archive.
    // A map lookup under a lock_guard there cost a mutex pair and a hash per
    // call for a value that never changes once read.
    //
    // Filled lazily and published with release/acquire, so a filled slot is
    // read with one atomic load and no lock.  The mutex still guards FILLING,
    // which keeps the class safe for the shared use its header promises even
    // though each Parquet object owns its own index today.
    if (row_group < 0 || row_group >= metadata_->num_row_groups() ||
        column < 0 || column >= columns_) {
      return nullptr;
    }
    Slot& slot = slots_[static_cast<std::size_t>(row_group) *
                            static_cast<std::size_t>(columns_) +
                        static_cast<std::size_t>(column)];

    if (slot.ready.load(std::memory_order_acquire)) return slot.stats;

    std::lock_guard<std::mutex> guard(mutex_);
    if (slot.ready.load(std::memory_order_relaxed)) return slot.stats;

    std::shared_ptr<parquet::ColumnChunkMetaData> chunk =
        metadata_->RowGroup(row_group)->ColumnChunk(column);
    if (chunk && chunk->is_stats_set()) {
      std::shared_ptr<parquet::Statistics> stats(chunk->statistics());
      if (stats) {
        slot.chunk = std::move(chunk);
        slot.stats = std::move(stats);
      }
    }
    // Absence is published too: a column that simply has no statistics must
    // not be re-examined once per query for the life of the file.
    slot.ready.store(true, std::memory_order_release);
    return slot.stats;
  }

  bool sorted_ascending(int32_t row_group, int32_t column)
  {
    std::lock_guard<std::mutex> guard(mutex_);
    const std::pair<int32_t, int32_t> key{row_group, column};
    auto it = sorted_.find(key);
    if (it != sorted_.end()) return it->second;

    bool sorted = false;
    if (row_group >= 0 && row_group < metadata_->num_row_groups()) {
      for (const auto& cs : metadata_->RowGroup(row_group)->sorting_columns()) {
        if (cs.column_idx != column) continue;
        sorted = !cs.descending && !cs.nulls_first;
        break;
      }
    }
    sorted_.emplace(key, sorted);
    return sorted;
  }

  /// Is @p column declared sorted ascending in EVERY row group?  Computed
  /// once; the answer is a property of the file, and the fast path in plan()
  /// is only correct when it holds for all of them, not just the first.
  bool all_sorted_ascending(int32_t column)
  {
    std::lock_guard<std::mutex> guard(mutex_);
    if (auto it = all_sorted_.find(column); it != all_sorted_.end()) return it->second;
    bool all = metadata_->num_row_groups() > 0;
    for (int32_t g = 0; all && g < metadata_->num_row_groups(); ++g) {
      bool one = false;
      for (const auto& cs : metadata_->RowGroup(g)->sorting_columns()) {
        if (cs.column_idx != column) continue;
        one = !cs.descending && !cs.nulls_first;
        break;
      }
      all = one;
    }
    all_sorted_.emplace(column, all);
    return all;
  }

  /// The row group that last satisfied a query through this reader.  A reader
  /// walks a contiguous run of entities, so the next query almost always wants
  /// the same group or the one after it.
  /// Relaxed rather than under mutex_: the hint is a guess that plan() then
  /// validates, so it carries no ordering obligation, and taking the lock here
  /// would serialise every plan() call across threads.
  int32_t group_hint() const { return group_hint_.load(std::memory_order_relaxed); }
  void set_group_hint(int32_t g) { group_hint_.store(g, std::memory_order_relaxed); }

  std::shared_ptr<parquet::FileMetaData> metadata_;
  std::vector<int64_t> rows_;
  std::map<int32_t, bool> all_sorted_;
  std::atomic<int32_t> group_hint_{0};
  int32_t columns_ = 0;
  std::deque<Slot> slots_; ///< (row group, column) -> statistics; see get()
  std::map<std::pair<int32_t, int32_t>, bool> sorted_;
  std::mutex mutex_;
};

/******************************************************************************/
StatsIndex::StatsIndex(std::shared_ptr<parquet::FileMetaData> metadata)
    : impl_(std::make_unique<Impl>(std::move(metadata)))
{
}

/******************************************************************************/
StatsIndex::~StatsIndex() = default;

/******************************************************************************/
int32_t StatsIndex::row_group_count() const
{
  return static_cast<int32_t>(impl_->rows_.size());
}

/******************************************************************************/
int64_t StatsIndex::row_count(int32_t row_group) const
{
  if (row_group < 0 || static_cast<std::size_t>(row_group) >= impl_->rows_.size()) {
    return 0;
  }
  return impl_->rows_[static_cast<std::size_t>(row_group)];
}

/******************************************************************************/
std::shared_ptr<parquet::Statistics> StatsIndex::get(int32_t row_group,
                                                     int32_t column) const
{
  return impl_->get(row_group, column);
}

/******************************************************************************/
bool StatsIndex::sorted_ascending(int32_t row_group, int32_t column) const
{
  return impl_->sorted_ascending(row_group, column);
}

/******************************************************************************/
bool StatsIndex::all_sorted_ascending(int32_t column) const
{
  return impl_->all_sorted_ascending(column);
}

/******************************************************************************/
int32_t StatsIndex::group_hint() const { return impl_->group_hint(); }
void StatsIndex::set_group_hint(int32_t g) const { impl_->set_group_hint(g); }

/******************************************************************************/
// Helper cache for the parquet page index.
class IndexCache {
public:
  // Cache items.
  struct Item {
    Item() {}

    Item(const std::shared_ptr<parquet::ColumnIndex>& c,
         const std::shared_ptr<parquet::OffsetIndex>& o)
        : page_count(o->page_locations().size())
        , col_index(c)
        , off_index(o)
    {
    }

    // Return true if this item is valid.
    bool valid() const
    {
      return page_count > 0 && col_index != nullptr && off_index != nullptr;
    }

    std::size_t page_count = 0;
    std::shared_ptr<parquet::ColumnIndex> col_index = nullptr;
    std::shared_ptr<parquet::OffsetIndex> off_index = nullptr;
  };

  // Store pointers so we can keep the items updated.
  using value_type = std::shared_ptr<Item>;

  // Keyed on the column path.
  using key_type = int32_t;

  /// Constructor.
  explicit IndexCache(const std::shared_ptr<parquet::RowGroupMetaData>& rg,
                      const std::shared_ptr<parquet::RowGroupPageIndexReader>& ri)
      : rg_(rg)
      , ri_(ri)
      , cache_()
  {
  }

  // Construct a key for the given group/field.
  //
  // FIXME: Should probably move this into the destination code.
  key_type key(const Schema::Column& dest) const
  {
    return dest.second->absolute_index();
  }

  // Return an item from the cache, or compute it.
  std::optional<value_type> get(const Schema::Column& dest)
  {
    key_type k(key(dest));
    auto it = cache_.find(k);

    if (it != cache_.end()) {
      if (it->second->valid()) {
        return it->second;
      } else {
        return {};
      }
    } else {
      value_type v = load(k).value_or(std::make_shared<Item>());
      cache_[k] = v;

      if (v->valid()) {
        return v;
      } else {
        return {};
      }
    }
  }

  // Return the page index with the maximum number of pages.
  std::optional<value_type> max_page_count() const
  {
    if (cache_.empty()) return {};
    value_type result = nullptr;

    for (const auto& kv : cache_) {
      if (!kv.second->valid()) return {};

      if (result == nullptr || kv.second->page_count > result->page_count) {
        result = kv.second;
      }
    }

    if (result != nullptr) {
      return result;
    } else {
      return {};
    }
  }

  // Generate ranges and insert them into the given vector.
  void generate_ranges(const std::vector<bool>& pages,
                       std::vector<Planner::Range>& ranges,
                       int32_t row_group,
                       value_type major_index)
  {
    auto on_span = [&](std::size_t first_page, std::size_t last_page) -> void {
      auto page_map = major_index->off_index->page_locations();

      int64_t first_row_index = page_map[first_page].first_row_index;

      // The index of the last row is computed by using the first row
      // of the next page.  But if we are already on the last page we
      // just use the total number of rows in the row group to get the
      // last row index.
      int64_t last_row_index = (last_page + 1 < major_index->page_count)
                                   ? page_map[last_page + 1].first_row_index - 1
                                   : rg_->num_rows() - 1;

      ranges.push_back(
          {row_group, first_row_index, last_row_index - first_row_index + 1});
    };

    Algorithm::spans(pages, on_span);
  }

private:
  std::optional<value_type> load(key_type column_index)
  {
    try {
      std::shared_ptr<parquet::ColumnIndex> cindex =
          ri_->GetColumnIndex(column_index);
      if (cindex == nullptr) return {};

      std::shared_ptr<parquet::OffsetIndex> oindex =
          ri_->GetOffsetIndex(column_index);
      if (oindex == nullptr) return {};

      return std::make_shared<Item>(cindex, oindex);
    } catch (const parquet::ParquetException&) {
      // The page index is flaky, especially if the page index was
      // generated by a different version of Arrow/Parquet.
      return {};
    }
  }

  // Internal vars.
  std::shared_ptr<parquet::RowGroupMetaData> rg_;
  std::shared_ptr<parquet::RowGroupPageIndexReader> ri_;

  // The actual cache.
  std::map<key_type, value_type> cache_;
};

/******************************************************************************/
/**
 * Fetch the minimum and maximum values for a column using either
 * column statistics or the page index.
 */
struct ColMinMax final {
  // Use column statistics.
  ColMinMax(std::shared_ptr<parquet::Statistics> stats)
      : stats_(std::move(stats))
      , index_(nullptr)
  {
  }

  // Use the page index.
  ColMinMax(std::shared_ptr<parquet::ColumnIndex>& index, std::size_t page)
      : stats_(nullptr)
      , index_(index)
      , page_index_(page)
  {
  }

  /// Any of these may be null.
  std::shared_ptr<parquet::Statistics> stats_;
  std::shared_ptr<parquet::ColumnIndex> index_;

  /// When using a column index this keeps track of which page we are
  /// currently querying.
  std::size_t page_index_ = 0;

  template <Type T> Query::Result<Query::range_t> operator()() const
  {
    using V = typename type_traits<T>::value_type;
    using P = typename type_traits<T>::parquet_type;
    using R = Query::Result<Query::range_t>;

    if (index_ != nullptr) {
      using Index = parquet::TypedColumnIndex<P>;
      std::shared_ptr<Index> index = std::static_pointer_cast<Index>(index_);
      V min = safe_cast_or_copy<V>(index->min_values()[page_index_]);
      V max = safe_cast_or_copy<V>(index->max_values()[page_index_]);
      return R(std::make_pair(min, max));
    } else if (stats_ != nullptr) {
      using Stats = parquet::TypedStatistics<P>;
      if (!stats_->HasMinMax()) return R::skip();
      std::shared_ptr<Stats> stats = std::static_pointer_cast<Stats>(stats_);
      V min = safe_cast_or_copy<V>(stats->min());
      V max = safe_cast_or_copy<V>(stats->max());
      return R(std::make_pair(min, max));
    } else {
      return R::skip();
    }
  }
};

// We don't support string searches right now.
template <>
Query::Result<Query::range_t> ColMinMax::operator()<Util::Type::ByteArray>() const
{
  return Query::Result<Query::range_t>::fail();
}

// RDR-4a: UInt64 columns are stored with a signed INT64 physical type, so
// the min/max raw values from statistics/page-index are bit-pattern signed.
// bit_cast them to uint64_t so the range type matches the predicate value type.
template <>
Query::Result<Query::range_t> ColMinMax::operator()<Util::Type::UInt64>() const
{
  using R = Query::Result<Query::range_t>;

  if (index_ != nullptr) {
    using Index = parquet::TypedColumnIndex<parquet::Int64Type>;
    std::shared_ptr<Index> index = std::static_pointer_cast<Index>(index_);
    return R(
        std::make_pair(std::bit_cast<uint64_t>(index->min_values()[page_index_]),
                       std::bit_cast<uint64_t>(index->max_values()[page_index_])));
  } else if (stats_ != nullptr) {
    using Stats = parquet::TypedStatistics<parquet::Int64Type>;
    if (!stats_->HasMinMax()) return R::skip();
    std::shared_ptr<Stats> stats = std::static_pointer_cast<Stats>(stats_);
    return R(std::make_pair(std::bit_cast<uint64_t>(stats->min()),
                            std::bit_cast<uint64_t>(stats->max())));
  } else {
    return R::skip();
  }
}
/******************************************************************************/
/**
 * Helper class for planning which row groups/pages need to be read by
 * the query executor.
 */
class Planner::Impl final {
public:
  /// Constructor.
  Impl(parquet::ParquetFileReader& reader,
       const Query& query,
       std::shared_ptr<StatsIndex> stats,
       std::shared_ptr<parquet::PageIndexReader> page_index)
      : reader_(&reader)
      , metadata_(reader.metadata())
      , page_index_reader_(std::move(page_index))
      , page_index_ready_(page_index_reader_ != nullptr)
      , stats_(std::move(stats))
      , plan_({query, {}})
  {
    // A planner is constructed per query; the index is shared across all of
    // them.  Building a private one here would reintroduce exactly the
    // per-query metadata cost the index exists to remove, so it is only a
    // fallback for callers that did not supply one.
    if (!stats_) stats_ = std::make_shared<StatsIndex>(metadata_);
  }

  /// Do the actual planning.
  const Planner::Plan& plan();

private:
  /**
   * Create a plan for the given row group.
   */
  void plan_row_group(int32_t row_group_index);

  /**
   * Test an entire row group using column statistics.
   *
   * If the there are no statistics returns a failure.  If all the
   * rows are null returns null.  Otherwise returns the query
   * evaluation result.
   */
  Query::Result<bool> with_column_stats(int32_t row_group_index);

  /**
   * Test the pages inside a row group.
   */
  Query::Result<bool>
  with_page_index(const std::shared_ptr<parquet::RowGroupMetaData>&,
                  const std::shared_ptr<parquet::RowGroupPageIndexReader>&,
                  int32_t);

private:
  /**
   * Mark the given row group as needing a full scan.
   */
  void full_scan(int32_t);

  /// The page-index reader, obtained ON DEMAND.
  ///
  /// Asking the file reader for one costs real time -- 79 of the 95
  /// thread-seconds of planning, once the group scan was fixed -- and a query
  /// that resolves to a whole group never needs it.  Caching one for the life
  /// of the file was tried and is worse; see Parquet::planner.
  std::shared_ptr<parquet::PageIndexReader> page_index()
  {
    if (!page_index_ready_) {
      page_index_reader_ = reader_->GetPageIndexReader();
      page_index_ready_ = true;
    }
    return page_index_reader_;
  }

  parquet::ParquetFileReader* reader_ = nullptr;
  std::shared_ptr<parquet::FileMetaData> metadata_;
  std::shared_ptr<parquet::PageIndexReader> page_index_reader_;
  bool page_index_ready_ = false;
  /// Set by plan() when the query resolves to whole row groups, so
  /// plan_row_group() does not consult the page index at all.
  bool whole_groups_ = false;
  std::shared_ptr<StatsIndex> stats_;
  Planner::Plan plan_;
};

/******************************************************************************/
Query::Result<bool> Planner::Impl::with_column_stats(int32_t row_group_index)
{
  // Evaluation callback that tries to use column statistics.
  auto via_stats = [this, row_group_index](
                       const Schema::Column& dest) -> Query::Result<Query::range_t> {
    auto stats = stats_->get(row_group_index, dest.second->absolute_index());
    if (stats != nullptr && dest.second->type().has_value()) {
      return lift_type(dest.second->type().value(), ColMinMax(std::move(stats)));
    } else {
      return Query::Result<Query::range_t>::skip();
    }
  };

  return plan_.query.eval(via_stats);
}

/******************************************************************************/
Query::Result<bool> Planner::Impl::with_page_index(
    const std::shared_ptr<parquet::RowGroupMetaData>& rg,
    const std::shared_ptr<parquet::RowGroupPageIndexReader>& ri,
    int32_t row_group_index)
{
  IndexCache index_cache(rg, ri);
  std::size_t current_page = 0;

  auto via_page_index =
      [&](const Schema::Column& dest) -> Query::Result<Query::range_t> {
    if (!dest.second->type().has_value())
      return Query::Result<Query::range_t>::fail();

    auto index = index_cache.get(dest);
    if (!index.has_value()) return Query::Result<Query::range_t>::skip();

    if (current_page >= index.value()->page_count) {
      // No more pages for this column.
      return Query::Result<Query::range_t>::skip();
    }

    if (index.value()->col_index->null_pages()[current_page]) {
      // This page only contains nulls.
      return Query::Result<Query::range_t>::skip();
    }

    return lift_type(dest.second->type().value(),
                     ColMinMax(index.value()->col_index, current_page));
  };

  // We need to loop once for each page in the page index.  But we
  // don't know how many pages there are yet.
  Query::Result<bool> qres = plan_.query.eval(via_page_index);
  if (qres.failed()) return qres;

  auto major_index = index_cache.max_page_count();
  if (!major_index.has_value() || major_index.value()->page_count < 1) {
    return Query::Result<bool>::fail();
  }

  // Track the query result for each page.  Start with all pages
  // marked as needed (true) so any errors or skips result in that
  // page being included.
  std::size_t page_count = major_index.value()->page_count;
  std::vector<bool> page_results(page_count, true);
  page_results[0] = !qres.is(false);

  for (current_page = 1; !qres.failed() && current_page < page_count;
       ++current_page) {
    qres = plan_.query.eval(via_page_index);
    page_results[current_page] = !qres.is(false);
  }

  // Record the page results.
  index_cache.generate_ranges(page_results, plan_.ranges, row_group_index,
                              major_index.value());
  return qres.to(true);
}

/******************************************************************************/
void Planner::Impl::plan_row_group(int32_t row_group_index)
{
  // The column statistics, if present, can tell us if we are able to
  // skip an entire row group.
  //
  // This runs FIRST, before the page index for this group is touched and
  // before any RowGroupMetaData is built.  A query that selects one spectrum
  // rejects every row group but one or two, and reading their page indexes to
  // then discard them was the dominant cost of planning on a file with
  // thousands of row groups.
  auto stats_result = with_column_stats(row_group_index);
  if (stats_result.has_value() && !stats_result.value()) {
    count_plan_pruned();
    return;
  }

  // A query the executor can resolve inside the group by itself needs only to
  // know WHICH group.  Narrowing further is wasted: the page index of a
  // heavily-RLE'd sorted key column narrows to ~56% of a row group anyway, and
  // the reader now finds its record batch from cached per-batch key spans.
  if (whole_groups_) {
    count_plan_full_scan();
    count_plan_range(stats_->row_count(row_group_index));
    plan_.ranges.push_back({row_group_index, 0, stats_->row_count(row_group_index)});
    return;
  }

  std::shared_ptr<parquet::RowGroupPageIndexReader> row_index_reader = nullptr;
  auto page_index_reader_local = page_index();

  if (page_index_reader_local != nullptr) {
    row_index_reader = page_index_reader_local->RowGroup(row_group_index);
  }

  if (page_index_reader_local == nullptr || row_index_reader == nullptr) {
    count_plan_pi_null();
    // Record this row group if the statistic planner selected it.  If
    // the plan failed we fall back to a full row group scan and
    // record it as well.
    if (!stats_result.is(false)) {
      full_scan(row_group_index);
    }
  } else {
    count_plan_page_index();
    std::shared_ptr<parquet::RowGroupMetaData> rg(
        metadata_->RowGroup(row_group_index));

    std::size_t before_count = plan_.ranges.size();
    auto page_res = with_page_index(rg, row_index_reader, row_group_index);

    // The page index code should have already created the necessary
    // range record even if it failed.  This is a "just in case"
    // check.
    if (!page_res.has_value() && plan_.ranges.size() == before_count) {
      full_scan(row_group_index);
    } else {
      // Whatever the page index narrowed this group to: this is the number
      // that decides how many record batches the executor must slice.
      for (std::size_t r = before_count; r < plan_.ranges.size(); ++r) {
        count_plan_range(plan_.ranges[r].length);
      }
    }
  }
}

/******************************************************************************/
void Planner::Impl::full_scan(int32_t row_group_index)
{
  // A whole-row-group scan.  This is the correct fallback when a group has no
  // page index to narrow within it -- the reference writer emits none, so real
  // files take this path -- and the executor binary-searches a declared-sorted
  // index within the range, so it is no longer O(group) per spectrum.  No
  // warning: a fresh Planner is built per fetch, so a per-call stderr line here
  // is one write per spectrum (hundreds of thousands on a real run) for a path
  // that is neither an error nor slow.
  count_plan_full_scan();
  count_plan_range(stats_->row_count(row_group_index));
  plan_.ranges.push_back({row_group_index, 0, stats_->row_count(row_group_index)});
}

/******************************************************************************/
const Planner::Plan& Planner::Impl::plan()
{
  const int32_t groups = metadata_->num_row_groups();

  // FAST PATH: an equality on a column this file declares sorted ascending in
  // every row group.
  //
  // Then the row groups that can match are CONTIGUOUS, so the scan may stop at
  // the first non-match after a match instead of evaluating every group.  That
  // matters because the general loop below evaluates the statistics of ALL of
  // them for EVERY entity: measured at 276,635,040 evaluations for one
  // 762,016-spectrum run over a 363-group archive, about two thirds of the
  // whole parallel read phase.
  //
  // Contiguous, NOT unique: 64 of those 762,016 spectra straddle a group
  // boundary and legitimately produce two ranges, so stopping at the FIRST
  // match would silently drop points.  The scan stops only once a match has
  // been seen and the next group does not match.
  //
  // Starting from the last group that matched makes the common case one
  // evaluation: a reader walks entities in order and ~2,100 of them share a
  // group.  If nothing matches from there to the end, the groups before it are
  // still searched, so the result cannot depend on the hint -- only the cost.
  bool ordered = false;
  if (auto equality = plan_.query.as_equality()) {
    if (equality->first.second->type().has_value()) {
      ordered = stats_->all_sorted_ascending(equality->first.second->absolute_index());
    }
  }

  if (ordered && groups > 0) {
    whole_groups_ = true; // the executor narrows within the group by key
    bool found = false;
    auto scan = [&](int32_t from, int32_t to) {
      for (int32_t i = from; i < to; ++i) {
        const std::size_t before = plan_.ranges.size();
        plan_row_group(i);
        if (plan_.ranges.size() > before) {
          found = true;
          stats_->set_group_hint(i);
        } else if (found) {
          return; // matches are contiguous: nothing further can match
        }
      }
    };
    const int32_t hint = std::clamp(stats_->group_hint(), 0, groups - 1);
    scan(hint, groups);

    // A match AT the hint says nothing about the group BEFORE it.
    //
    // An entity straddling a group boundary lives in two groups, and the
    // forward scan above starts at the hint -- so fetching such an entity
    // after the hint has moved to its SECOND group found only that half and
    // returned a spectrum missing its first rows.  The forward scan cannot
    // see this: it stops at the first non-match after a match, and there is
    // no match before its starting point to stop at.  Reachable from any
    // non-monotonic access order (a subsampled or random-access read, a
    // re-read of the same spectrum) and silent whenever the file declares no
    // per-spectrum count to check the short array against.
    //
    // Matches are contiguous, so walking DOWN from the hint until a group
    // does not match costs one extra statistics evaluation in the common case
    // (the group before the hint holds a different entity and stops it at
    // once), and only runs when the hint itself matched.
    if (found && hint > 0) {
      for (int32_t i = hint - 1; i >= 0; --i) {
        const std::size_t before = plan_.ranges.size();
        plan_row_group(i);
        if (plan_.ranges.size() == before) break;
      }
    }

    if (!found) scan(0, hint);
  } else {
    for (int32_t i : std::views::iota(0, groups)) {
      plan_row_group(i);
    }
  }

  std::ranges::sort(plan_.ranges, {}, &Range::row_group);
  return plan_;
}

/******************************************************************************/
Planner::Planner(parquet::arrow::FileReader& reader,
                 const Query& query,
                 std::shared_ptr<StatsIndex> stats,
                 std::shared_ptr<parquet::PageIndexReader> page_index)
    : impl_(std::make_unique<Impl>(*reader.parquet_reader(), query, std::move(stats),
                                   std::move(page_index)))
{
}

/******************************************************************************/
Planner::~Planner() = default;

/******************************************************************************/
const Planner::Plan& Planner::plan() const { return impl_->plan(); }

} // namespace MzPeak::Util
