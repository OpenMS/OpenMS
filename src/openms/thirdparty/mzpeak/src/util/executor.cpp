/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include <atomic>
#include <chrono>
#include <limits>
#include <optional>

// Diagnostic counters. OFF by default: the planner one fires hundreds of
// millions of times per run, and a global atomic there is itself a cost.
// Build with -DMZPEAK_READ_COUNTERS to measure.
#ifdef MZPEAK_READ_COUNTERS
#define MZPEAK_COUNT(c) (c).fetch_add(1, std::memory_order_relaxed)
#else
#define MZPEAK_COUNT(c) ((void)0)
#endif
#include "mzpeak/util/executor.h"

#include <algorithm>
#include <arrow/array.h>
#include <arrow/record_batch.h>
#include <arrow/result.h>
#include <memory>
#include <parquet/arrow/reader.h>
#include <ranges>
#include <string_view>
#include <type_traits>

#include "mzpeak/util/algorithm.h"
#include "mzpeak/util/key_runs.h"
#include "mzpeak/util/parquet.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Util {

namespace
{
  /// The query value as an int64, when it is an integer at all.  Used only to
  /// pick which record batches can hold a key; a non-integer value simply
  /// falls back to walking them.
  std::optional<int64_t> integer_value(const Query::value_t& v)
  {
    if (auto p = std::get_if<int64_t>(&v)) return *p;
    if (auto p = std::get_if<uint64_t>(&v)) return static_cast<int64_t>(*p);
    if (auto p = std::get_if<int32_t>(&v)) return *p;
    if (auto p = std::get_if<uint32_t>(&v)) return *p;
    if (auto p = std::get_if<int8_t>(&v)) return *p;
    if (auto p = std::get_if<uint8_t>(&v)) return *p;
    return std::nullopt;
  }
} // namespace


namespace
{
  std::atomic<long> g_plan_group_evals{0}, g_batches_visited{0}, g_slices_made{0};
  std::atomic<long> g_plan_pruned{0}, g_plan_full_scan{0}, g_plan_page_index{0};
  std::atomic<long> g_plan_pi_null{0}, g_plan_ranges{0}, g_plan_range_rows{0};
  std::atomic<long> g_ns_plan_ctor{0};
  std::atomic<long> g_ns_plan{0}, g_ns_exec{0}, g_ns_rowgroup{0}, g_ns_project{0};
}

ReadCounters read_counters()
{
  return {g_plan_group_evals.load(std::memory_order_relaxed),
          g_batches_visited.load(std::memory_order_relaxed),
          g_slices_made.load(std::memory_order_relaxed),
          g_plan_pruned.load(std::memory_order_relaxed),
          g_plan_full_scan.load(std::memory_order_relaxed),
          g_plan_page_index.load(std::memory_order_relaxed),
          g_plan_pi_null.load(std::memory_order_relaxed),
          g_plan_ranges.load(std::memory_order_relaxed),
          g_plan_range_rows.load(std::memory_order_relaxed),
          g_ns_plan.load(std::memory_order_relaxed),
          g_ns_plan_ctor.load(std::memory_order_relaxed),
          g_ns_exec.load(std::memory_order_relaxed),
          g_ns_rowgroup.load(std::memory_order_relaxed),
          g_ns_project.load(std::memory_order_relaxed)};
}

void count_plan_group_eval() { MZPEAK_COUNT(g_plan_group_evals); }
void count_plan_pruned() { MZPEAK_COUNT(g_plan_pruned); }
void count_plan_full_scan() { MZPEAK_COUNT(g_plan_full_scan); }
void count_plan_page_index() { MZPEAK_COUNT(g_plan_page_index); }
void count_plan_pi_null() { MZPEAK_COUNT(g_plan_pi_null); }
void count_ns_plan([[maybe_unused]] long ns)
{
#ifdef MZPEAK_READ_COUNTERS
  g_ns_plan.fetch_add(ns, std::memory_order_relaxed);
#endif
}
void count_ns_plan_ctor([[maybe_unused]] long ns)
{
#ifdef MZPEAK_READ_COUNTERS
  g_ns_plan_ctor.fetch_add(ns, std::memory_order_relaxed);
#endif
}
void count_ns_exec([[maybe_unused]] long ns)
{
#ifdef MZPEAK_READ_COUNTERS
  g_ns_exec.fetch_add(ns, std::memory_order_relaxed);
#endif
}
void count_ns_rowgroup([[maybe_unused]] long ns)
{
#ifdef MZPEAK_READ_COUNTERS
  g_ns_rowgroup.fetch_add(ns, std::memory_order_relaxed);
#endif
}
void count_ns_project([[maybe_unused]] long ns)
{
#ifdef MZPEAK_READ_COUNTERS
  g_ns_project.fetch_add(ns, std::memory_order_relaxed);
#endif
}

void count_plan_range([[maybe_unused]] long rows)
{
  MZPEAK_COUNT(g_plan_ranges);
#ifdef MZPEAK_READ_COUNTERS
  g_plan_range_rows.fetch_add(rows, std::memory_order_relaxed);
#endif
}


/******************************************************************************/
struct Executor::Impl {
  Impl(Parquet& source, const Projection& fields)
      : source_(source)
      , projection_(fields)
      , slice_(nullptr)
  {
  }

  /// Helper to check a result and throw an error if necessary.
  template <typename T> T check(std::string_view msg, arrow::Result<T> r)
  {
    if (!r.ok()) {
      std::string error("while executing a query: " + std::string(msg) + ": ");
      throw ParquetError(error + r.status().ToString());
    }

    return std::move(r.ValueOrDie());
  }

  /// Run a query on the given record batch and return a vector
  /// indicating which rows should be kept.
  ///
  std::vector<bool> filter(const Planner::Plan&,
                           std::shared_ptr<arrow::RecordBatch>&,
                           bool sorted_column,
                           std::pair<int64_t, int64_t>& run,
                           bool& have_run);

  /// Capture the requested columns.
  void project(std::shared_ptr<arrow::RecordBatch>&);

  /// The rows [first, last) of a group held as runs that @p plan selects --
  /// only when the plan is the one the sorted binary search below would answer
  /// (an equality on this key, declared sorted in this group, in the key's own
  /// type) and nothing projects the key.  Otherwise nothing.
  std::optional<std::pair<int64_t, int64_t>>
  key_run(const Planner::Plan&, const KeyRuns&, bool sorted_column) const;

  /// Return the array associated with the given field.
  std::shared_ptr<arrow::Array> array(std::shared_ptr<arrow::RecordBatch>& batch,
                                      const Schema::Column& field);

  /// Helper to extract a single value from an array.
  struct ArrayValueHelper {
    template <Type T> Query::Result<Query::value_t> operator()();
    Executor::Impl& impl_;
    const Schema::Column& field_;
    std::shared_ptr<arrow::RecordBatch>& batch_;
    int64_t row_index_;
  };

  Parquet& source_;
  Projection projection_;
  std::unique_ptr<Slice> slice_;
};

/******************************************************************************/
template <Type T>
Query::Result<Query::value_t> Executor::Impl::ArrayValueHelper::operator()()
{
  std::shared_ptr<arrow::Array> a(impl_.array(batch_, field_));

  if (row_index_ < 0 || row_index_ >= a->length() || a->IsNull(row_index_)) {
    return Query::Result<Query::value_t>::skip();
  }

  auto casted = std::static_pointer_cast<typename type_traits<T>::array_type>(a);
  auto value = Decoders::unsafe_array_value<T>(casted, row_index_);
  return Query::Result<Query::value_t>(value);
}

/******************************************************************************/
/*
 * Mark the rows of one typed column equal to `wanted`.
 *
 * The general path below resolves the column, dispatches on its type and walks
 * the predicate tree once PER ROW.  For the query this library actually runs --
 * "the entity index equals N" -- that is several hundred nanoseconds of
 * machinery around a single integer comparison, repeated across a whole page.
 */
struct EqualityScan {
  std::shared_ptr<arrow::Array> array;
  Query::value_t wanted;
  std::vector<bool>& want_rows;
  bool& handled;
  /// The row group DECLARES this column sorted ascending, nulls last.
  bool sorted;
  /// When the sorted path is taken, the matching run as [begin, end).
  std::pair<int64_t, int64_t>& run;
  bool& have_run;

  template <Type T> void operator()() const
  {
    using value_type = typename type_traits<T>::value_type;

    // A predicate whose value type does not match the column's is not something
    // this scan can answer; fall back rather than silently selecting nothing.
    if (!std::holds_alternative<value_type>(wanted)) return;

    const value_type target = std::get<value_type>(wanted);
    auto typed =
        std::static_pointer_cast<typename type_traits<T>::array_type>(array);
    const int64_t length = typed->length();

    // Binary-search path.  When the file DECLARES this column sorted ascending
    // with nulls last AND the batch has no nulls, the matching rows are one
    // contiguous run, found in O(log n).  This is what keeps a file with NO
    // page index -- so the planner hands over a whole 1M-row group -- from
    // costing O(group) per spectrum; the linear scan below is ~630x that on a
    // real Astral run.
    //
    // Guarded two ways, both required.  On the declaration, because searching a
    // column that is not actually sorted finds one run and silently misses the
    // rest.  On null_count, because the general path KEEPS null rows and a
    // sorted column's nulls form a SECOND run at the end -- two disjoint
    // intervals a single [begin, end) cannot express.  Neither guard holds for
    // an entity index, which is required and sorted, so the fast path is taken
    // for exactly the query this reader runs.
    if constexpr (std::is_arithmetic_v<value_type>) {
      if (sorted && typed->null_count() == 0) {
        const value_type* values = typed->raw_values();
        const value_type* lo = std::lower_bound(values, values + length, target);
        const value_type* hi = std::upper_bound(lo, values + length, target);
        run = {lo - values, hi - values};
        have_run = true;
        handled = true;
        return;
      }
    }

    want_rows.assign(static_cast<std::size_t>(length), true);
    for (int64_t i = 0; i < length; ++i) {
      // A null never compares equal, and the general path treats an unreadable
      // value as "no opinion" -- it yields no `false`, so the row stays
      // selected.  Matched here so the two paths agree row for row.
      want_rows[static_cast<std::size_t>(i)] =
          typed->IsNull(i) || typed->Value(i) == target;
    }
    handled = true;
  }
};

/******************************************************************************/
std::vector<bool> Executor::Impl::filter(const Planner::Plan& plan,
                                         std::shared_ptr<arrow::RecordBatch>& batch,
                                         bool sorted_column,
                                         std::pair<int64_t, int64_t>& run,
                                         bool& have_run)
{
  // Each path sizes the mask itself, so a fast path that declines costs no
  // allocation before the general path takes over.
  std::vector<bool> want_rows;

  if (auto equality = plan.query.as_equality()) {
    const auto& [field, wanted] = *equality;
    if (field.second->type().has_value()) {
      std::shared_ptr<arrow::Array> column(array(batch, field));
      bool handled = false;
      lift_type(field.second->type().value(),
                EqualityScan{column, wanted, want_rows, handled, sorted_column, run,
                             have_run});
      if (have_run) return {};
      if (handled) return want_rows;
    }
  }

  want_rows.assign(static_cast<std::size_t>(batch->num_rows()), true);

  auto get_value =
      [&](int64_t row_index,
          const Schema::Column& field) -> Query::Result<Query::value_t> {
    if (field.second->type().has_value()) {
      return lift_type(field.second->type().value(),
                       ArrayValueHelper{*this, field, batch, row_index});
    } else {
      return Query::Result<Query::value_t>::fail();
    }
  };

  for (int64_t i : std::views::iota(0, batch->num_rows())) {
    Query::Result<bool> res =
        plan.query.eval(std::bind(get_value, i, std::placeholders::_1));
    want_rows[i] = !res.is(false);
  }

  return want_rows;
}

/******************************************************************************/
std::optional<std::pair<int64_t, int64_t>>
Executor::Impl::key_run(const Planner::Plan& plan, const KeyRuns& runs,
                        bool sorted_column) const
{
  if (!sorted_column) return std::nullopt;
  auto equality = plan.query.as_equality();
  if (!equality) return std::nullopt;
  const auto& [field, wanted] = *equality;
  if (field.second->absolute_index() != runs.leaf) return std::nullopt;

  // The type the binary search would compare in (EqualityScan), which must be
  // the order the runs were verified in.
  const auto& type = field.second->type();
  if (!type.has_value() || *type != (runs.is_unsigned ? Type::UInt64 : Type::Int64)) {
    return std::nullopt;
  }
  for (const auto& column : projection_.get()) {
    if (column.second->absolute_index() == runs.leaf) return std::nullopt;
  }
  return runs.rows(wanted);
}

/******************************************************************************/
void Executor::Impl::project(std::shared_ptr<arrow::RecordBatch>& batch)
{
  for (const auto& field : slice_->fields()) {
    slice_->append(field, array(batch, field));
  }
}

/******************************************************************************/
std::shared_ptr<arrow::Array>
Executor::Impl::array(std::shared_ptr<arrow::RecordBatch>& batch,
                      const Schema::Column& field)
{
  // Boxed here from the batch's ArrayData, NOT through RecordBatch::column()
  // and StructArray::field().  Those two keep a boxed Array cached inside the
  // batch and read that cache with std::atomic_load on a shared_ptr, which
  // libstdc++ implements with a pool of sixteen PROCESS-WIDE mutexes: every
  // reader thread took the same sixteen locks several times per spectrum,
  // whatever group it read, and slept on them once enough threads did.
  // Measured on a 717,924-spectrum archive with decode waits out of the way:
  // a read cost ~12 us at 16 threads and 115-135 us at 256, with half a
  // million voluntary context switches.  Boxing a fresh Array costs an
  // allocation and takes no lock.
  if (field.first->is_root()) {
    return arrow::MakeArray(batch->column_data(field.second->absolute_index()));
  }

  const std::shared_ptr<arrow::ArrayData>& col = batch->column_data(field.first->index());
  if (!col || col->type->id() != arrow::Type::STRUCT) {
    throw ParquetError("column not in batch: " + field.first->path(*field.second));
  }

  // Exactly what StructArray::field() does: a struct's children are not
  // sliced along with it, so the struct's window is applied here.
  std::shared_ptr<arrow::ArrayData> child =
      col->child_data[static_cast<std::size_t>(field.second->relative_index())];
  if (col->offset != 0 || child->length != col->length) {
    child = child->Slice(col->offset, col->length);
  }
  return arrow::MakeArray(child);
}

/******************************************************************************/
Executor::Executor(Parquet& source, const Projection& fields)
    : impl_(std::make_unique<Impl>(source, fields))
{
}

/******************************************************************************/
Executor::~Executor() = default;

/******************************************************************************/
std::unique_ptr<Executor::Slice> Executor::execute(const Planner::Plan& plan)
{
#ifdef MZPEAK_READ_COUNTERS
  const auto t_exec0 = std::chrono::steady_clock::now();
  struct ExecTimer {
    std::chrono::steady_clock::time_point t0;
    ~ExecTimer()
    {
      count_ns_exec(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count());
    }
  } exec_timer{t_exec0};
#endif
  impl_->slice_ = std::unique_ptr<Slice>(new Slice(impl_->projection_.get()));
  std::map<Schema::Group::index_type, std::vector<Planner::Range>> ranges;

  for (const auto& range : plan.ranges) {
    auto it = ranges.find(range.row_group);

    if (it != ranges.end()) {
      it->second.push_back(range);
    } else {
      ranges[range.row_group] = {range};
    }
  }

  // Reduce a batch to just those rows that were selected by the query
  // and then capture the projected columns.
  auto project_wanted = [&](std::size_t i, std::size_t j,
                            std::shared_ptr<arrow::RecordBatch>& batch) -> void {
    int64_t offset = i;
    int64_t length = j - i + 1;
    if (length <= 0) return;

    auto sliced = batch->Slice(offset, length);
    impl_->project(sliced);
  };

  for (const auto& row_group : ranges) {
    int64_t row_group_start = 0;

    // The last row any range of this group wants.  Batches past it hold
    // nothing for this query.
    int64_t wanted_end = 0;
    // ...and the first row any range wants.  Batches that end before it hold
    // nothing either, and skipping them is an integer add.
    int64_t wanted_begin = std::numeric_limits<int64_t>::max();
    for (const auto& range : row_group.second) {
      wanted_end = std::max(wanted_end, range.offset + range.length);
      wanted_begin = std::min(wanted_begin, range.offset);
    }

    // Does THIS row group declare the predicate's column sorted?  Asked once
    // per group; the answer lets filter() binary-search instead of scanning the
    // whole group when the file carries no page index to narrow it.
    bool sorted_column = false;
    if (auto equality = plan.query.as_equality()) {
      sorted_column = impl_->source_.sorted_ascending(
          row_group.first, equality->first.second->absolute_index());
    }

    // Decoded once per row group and retained, so reading a run entity by
    // entity costs one decode per group rather than one per entity.
#ifdef MZPEAK_READ_COUNTERS
    const auto t_rg0 = std::chrono::steady_clock::now();
#endif
    auto batches = impl_->source_.row_group(row_group.first);
#ifdef MZPEAK_READ_COUNTERS
    count_ns_rowgroup(std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - t_rg0)
                          .count());
#endif

    // A group decoded WITHOUT its sorted key (KeyRuns): the runs name the
    // wanted rows directly.  The walk below is the sorted-equality path further
    // down, run for run -- same batches, same ranges, same order, same slices
    // -- with the run found in the group's runs instead of by a binary search
    // over the column it replaces.  Any other query, or one that projects the
    // key itself, gets the group back as a full decode would have produced it.
    if (batches->key_runs) {
      if (auto run = impl_->key_run(plan, *batches->key_runs, sorted_column)) {
        const auto [first, last] = *run;
        const auto& offsets = batches->row_offset;
        auto bi = static_cast<std::size_t>(
            std::upper_bound(offsets.begin(), offsets.end(), first) - offsets.begin());
        for (bi = bi > 0 ? bi - 1 : 0;
             first < last && bi < batches->batches.size() && offsets[bi] < last; ++bi) {
          const std::shared_ptr<arrow::RecordBatch>& batch = batches->batches[bi];
          const int64_t start = offsets[bi];
          const int64_t rows = batch->num_rows();
          MZPEAK_COUNT(g_batches_visited);
          for (const auto& range : row_group.second) {
            if (range.offset + range.length <= start || range.offset >= start + rows) {
              continue;
            }
            const auto [offset, length] =
                Algorithm::intersect_range(range.offset, range.length, start, rows);
            if (length == 0) continue;
            const int64_t lo = std::max(first - start, offset);
            const int64_t hi = std::min(last - start, offset + length);
            if (hi > lo) {
              MZPEAK_COUNT(g_slices_made);
              auto matched = batch->Slice(lo, hi - lo);
              impl_->project(matched);
            }
          }
        }
        continue;
      }
      batches = restore_key_column(*batches);
    }

    // WHICH batches can hold the wanted key, decided arithmetically.
    //
    // The decoded group carries the first and last key of every batch, taken
    // once by the thread that decoded it. When the query is an equality on the
    // sorted column, two binary searches over those give exactly the batches
    // that can match -- typically one of ~17 -- and the rest are never touched
    // at all: no shared_ptr copy, no Arrow accessor, no filter.
    //
    // That last part is the point. Arrow's RecordBatch::column() and
    // StructArray::field() are not free once boxed: libstdc++ implements the
    // atomic shared_ptr load they use with a pool of 16 process-wide mutexes,
    // so every visited batch was contending with every other thread's visits
    // regardless of which row group each was reading. 9.67 visits per spectrum
    // over 762,016 spectra is 7.4 million of them.
    std::size_t b_begin = 0;
    std::size_t b_end = batches->batches.size();
    if (sorted_column && batches->has_keys()) {
      if (auto equality = plan.query.as_equality()) {
        // The spans describe the column the FILE declares sorted. If the query
        // asks about a different column, they say nothing about it.
        const bool same_column =
            equality->first.second->absolute_index() == batches->key_leaf;
        if (auto wanted = same_column ? integer_value(equality->second)
                                      : std::nullopt) {
          // key_last is non-decreasing across batches: the first batch that
          // could contain the key is the first whose last key reaches it.
          const auto it = std::lower_bound(batches->key_last.begin(),
                                           batches->key_last.end(), *wanted);
          b_begin = static_cast<std::size_t>(it - batches->key_last.begin());
          b_end = b_begin;
          while (b_end < batches->key_first.size() &&
                 batches->key_first[b_end] <= *wanted) {
            ++b_end;
          }
        }
      }
    }

    for (std::size_t bi = b_begin; bi < b_end; ++bi) {
      const auto& cached = batches->batches[bi];
      row_group_start = batches->row_offset.empty()
                            ? row_group_start
                            : batches->row_offset[bi];
      if (row_group_start >= wanted_end) break;
      const int64_t rows = cached->num_rows();

      // Walk PAST the leading batches rather than through them.  A group holds
      // ~17 batches of 65,536 rows and one spectrum's few hundred points sit
      // inside one of them, but the traversal started at batch 0 every time:
      // measured at 13.0 batches visited and 9.7 slices made per spectrum, for
      // a run of 762,016 spectra.  Each visit copied the batch's shared_ptr
      // and each slice copied the ownership of all three columns -- and Arrow
      // shares one datatype object per primitive type, so those refcounts are
      // TRUE sharing between every thread, whatever group it is reading.
      if (row_group_start + rows <= wanted_begin) {
        row_group_start += rows;
        continue;
      }
      MZPEAK_COUNT(g_batches_visited);

      // By reference: the cache owns these batches for the whole call.
      const std::shared_ptr<arrow::RecordBatch>& batch = cached;

      for (const auto& range : row_group.second) {
        // Is this batch within this range?
        if (range.offset + range.length <= row_group_start ||
            range.offset >= row_group_start + rows) {
          continue;
        }

        const auto [offset, length] = Algorithm::intersect_range(
            range.offset, range.length, row_group_start, rows);
        if (length == 0) continue;

        // FAST PATH (under diagnosis): filter the UNSLICED batch and slice
        // only the run that matched. See doc/HANDOFF-page-index-resolution.md.
        if (sorted_column && plan.query.as_equality()) {
          std::pair<int64_t, int64_t> run{0, 0};
          bool have_run = false;
          std::shared_ptr<arrow::RecordBatch> whole = batch;
          (void)impl_->filter(plan, whole, sorted_column, run, have_run);
          if (have_run) {
            const int64_t lo = std::max(run.first, offset);
            const int64_t hi = std::min(run.second, offset + length);
            if (hi > lo) {
              MZPEAK_COUNT(g_slices_made);
              auto matched = batch->Slice(lo, hi - lo);
              impl_->project(matched);
            }
            continue;
          }
        }

        MZPEAK_COUNT(g_slices_made);
        auto sliced = batch->Slice(offset, length);

        std::pair<int64_t, int64_t> run{0, 0};
        bool have_run = false;
        auto want_rows = impl_->filter(plan, sliced, sorted_column, run, have_run);

        if (have_run) {
          // The sorted fast path handed back one contiguous run; project it
          // directly rather than materialising and re-scanning an n-bit mask.
          if (run.second > run.first) {
            project_wanted(static_cast<std::size_t>(run.first),
                           static_cast<std::size_t>(run.second - 1), sliced);
          }
        } else {
          Algorithm::spans(want_rows, project_wanted, sliced);
        }
      }

      row_group_start += rows;
    }
  }

  return std::move(impl_->slice_);
}

} // namespace MzPeak::Util
