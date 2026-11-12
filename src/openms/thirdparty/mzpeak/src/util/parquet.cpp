/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/parquet.h"

#include <arrow/array.h>
#include <arrow/record_batch.h>
#include <arrow/util/key_value_metadata.h>
#include <boost/json.hpp>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <parquet/api/reader.h>
#include <parquet/arrow/reader.h>
#include <ranges>

#include <arrow/memory_pool.h>

#include "mzpeak/exception.h"
#include "mzpeak/util/arrow.h"
#include "mzpeak/util/decode_pool.h"
#include "mzpeak/util/key_runs.h"

namespace MzPeak::Util {

/******************************************************************************/
std::optional<std::string> get_kv_string(const Parquet::file_metadata_t& fmd,
                                         std::string_view key)
{
  auto result(fmd->key_value_metadata()->Get(key));

  if (result.ok()) {
    return result.ValueOrDie();
  } else {
    return {};
  }
}

/******************************************************************************/
std::optional<std::size_t> get_kv_uint(const Parquet::file_metadata_t& fmd,
                                       std::string_view key)
{
  return get_kv_string(fmd, key).and_then(
      [](const std::string& s) -> std::optional<std::size_t> {
        std::size_t r{};
        auto [ptr, ec]{std::from_chars(s.data(), s.data() + s.size(), r)};

        // Require a clean, FULL parse: a partial match like "0junk" must not
        // be accepted as a present 0 (that would suppress the record_count
        // statistics fallback and report a wrong count -- RDR-11).
        if (ec == std::errc() && ptr == s.data() + s.size()) {
          return r;
        } else {
          return std::nullopt;
        }
      });
}

/******************************************************************************/
struct Parquet::Impl {
  Impl(std::unique_ptr<IO::File> data,
       Schema::File file,
       std::shared_ptr<RowGroupCache> shared_cache,
       Parquet::file_metadata_t footer)
      : file_(std::move(file))
      , arrow_(std::make_unique<Arrow>(std::move(data)))
      , reader_(nullptr)
      , groups_(std::make_shared<Schema::GroupMap>())
      , shared_cache_(std::move(shared_cache))
  {
    auto raf = arrow_->reader();

    auto reader_builder = parquet::arrow::FileReaderBuilder();
    // Decoded row groups come from the DecodePool, not Arrow's default.
    //
    // A decoded group outlives the thread that decoded it: the cache is
    // shared across readers (RowGroupCache), so a buffer allocated on a
    // worker is routinely freed later, often after that worker has exited --
    // a consumer using OpenMP tears its workers down at the end of each
    // parallel region.  Arrow 25's bundled mimalloc/jemalloc crash on that
    // free (EXC_BAD_ACCESS in _mi_arenas_page_abandon).  Plain malloc
    // survives it but strands the freed group in the decoding thread's arena
    // (3.56 GB at 128 threads); the DecodePool recycles it across threads.
    reader_builder.memory_pool(&decode_pool());
    // A footer parsed before is handed in rather than read and parsed again;
    // see Manager::parquet().  The page and decompression buffers come from
    // the same pool; every other reader property keeps its default.
    auto status = reader_builder.Open(
        std::move(raf), parquet::ReaderProperties(&decode_pool()), std::move(footer));

    if (!status.ok()) {
      std::string msg("while opening file: " + file_.file_name() + ": ");
      throw ParquetError(msg + status.ToString());
    }

    // NOT set_use_threads(true): measured flat (8.83 vs 8.78 ms/spectrum over
    // three runs each).  Arrow parallelises across COLUMNS within a read, and a
    // point signal table has three leaves, so there is almost nothing to
    // overlap.  Worth revisiting only for a file with many auxiliary arrays.
    std::unique_ptr<parquet::arrow::FileReader> reader;
    status = reader_builder.Build(&reader);

    if (!status.ok()) {
      std::string msg("while reading file: " + file_.file_name() + ": ");
      throw ParquetError(msg + status.ToString());
    }

    reader_ = std::move(reader);
    stats_ = std::make_shared<StatsIndex>(reader_->parquet_reader()->metadata());
    group_bytes_.assign(
        static_cast<std::size_t>(reader_->parquet_reader()->metadata()->num_row_groups()),
        0);
    parse_schema();
  }

  ~Impl() = default;

  void error(const std::string& error)
  {
    std::string msg("file accessing " + file_.file_name() + ": " + error);
    throw ParquetError(msg);
  }

  /// Load the schema.
  void parse_schema();

  /// Decode a row group in full and return its batches, caching the result:
  /// in the archive-wide shared cache when this object was given one, else
  /// in the private two-group cache below.
  std::shared_ptr<const Parquet::RowGroupBatches> row_group(int32_t index);

  /// The decode itself, under this object's lock; see decode_group_().
  std::shared_ptr<const Parquet::RowGroupBatches> decode_group_(int32_t index);

  /// decoded_row_group_bytes() of @p index, worked out once per group: it
  /// walks every column of the schema and cannot change.  Caller holds
  /// cache_mutex_.
  std::size_t group_bytes_locked_(int32_t index);

  Schema::File file_;
  std::unique_ptr<Arrow> arrow_;
  std::shared_ptr<parquet::arrow::FileReader> reader_;
  std::shared_ptr<Schema::GroupMap> groups_;

  /// Shared with every Parquet over the same archive (Manager::parquet()), so
  /// sixteen readers decode a group once, not sixteen times.  Null when this
  /// object was built outside a Manager; then cache_ below is used.
  std::shared_ptr<RowGroupCache> shared_cache_;

  /// The group this object last took from shared_cache_, HELD, so asking for
  /// it again does not go back to the cache.  See row_group().
  int32_t memo_group_ = -1;
  std::shared_ptr<const Parquet::RowGroupBatches> memo_batches_;

  /// group_bytes_locked_()'s answers, one per row group; 0 until asked.
  std::vector<std::size_t> group_bytes_;

  /// Private fallback: most-recently-decoded row groups, newest last.  Bounded
  /// because a decoded group is tens of megabytes; two is enough for a forward
  /// pass, where the only backward reach is an entity straddling a group
  /// boundary.
  static constexpr std::size_t kCachedGroups = 2;
  std::vector<std::pair<int32_t, std::shared_ptr<const Parquet::RowGroupBatches>>>
      cache_;

  /// Guards cache_ (private path), and the memo and group_bytes_ (shared
  /// path).  See row_group().
  std::mutex cache_mutex_;
  /// Guards the decode, which shares one file position; taken on both paths.
  std::recursive_mutex decode_mutex_;

  // Shared by every planner over this file; see StatsIndex.
  std::shared_ptr<StatsIndex> stats_;
};

/******************************************************************************/
void Parquet::Impl::parse_schema()
{
  auto fmd = reader_->parquet_reader()->metadata();
  auto root = fmd->schema()->group_node();
  int32_t offset = 0;

  std::shared_ptr<Schema::Group> root_group =
      std::make_shared<Schema::Group>(*root, file_);

  if (!root_group->fields().empty()) {
    (*groups_)[root_group->name()] = root_group;
  }

  for (int32_t i : std::views::iota(0, root->field_count())) {
    auto node = root->field(i);

    if (node->is_group() && !node->logical_type()->is_list()) {
      std::shared_ptr<parquet::schema::GroupNode> group =
          std::static_pointer_cast<parquet::schema::GroupNode>(node);

      std::shared_ptr<Schema::Group> s =
          std::make_shared<Schema::Group>(*group, file_, i, offset);
      (*groups_)[s->name()] = s;
      offset += group->field_count();
    }
  }
}

/******************************************************************************/
/// A row group's DECODED footprint: what an Arrow RecordBatch of it costs.
///
/// NOT Parquet's total_byte_size, which is the ENCODED size and can be several
/// times smaller.  Arrow materialises fixed-width leaves flat, so the decoded
/// size is rows x the leaf widths; on disk a sorted, repetitive column (a
/// spectrum index, say) is RLE'd to almost nothing and expands straight back
/// out in memory.  Measured on one archive: 4.5 MB encoded, 20.0 MB decoded,
/// a 4.4x under-count that let the shared cache hold 4.4x its budget.
///
/// ponytail: fixed-width leaves only.  BYTE_ARRAY has no width until it is
/// read, so a schema containing one falls back to the encoded figure and
/// under-counts that column; mzPeak's signal tables are all fixed-width.
std::size_t decoded_row_group_bytes(const parquet::RowGroupMetaData& group,
                                    const parquet::SchemaDescriptor& schema)
{
  const auto rows = static_cast<std::size_t>(group.num_rows());
  std::size_t width = 0;
  bool exact = true;
  // A key the decode will hold as runs costs nothing per row (KeyRuns).
  // ponytail: judged from the footer, so a file that declares a sorted key it
  // does not have is under-counted by that key; its runs (~16 B a spectrum)
  // are not counted at all.
  const std::optional<int32_t> runs = run_key_leaf(group, schema);

  for (int i : std::views::iota(0, schema.num_columns())) {
    if (runs && i == *runs) continue;
    const parquet::ColumnDescriptor* column = schema.Column(i);
    switch (column->physical_type()) {
      case parquet::Type::BOOLEAN: width += 1; break;
      case parquet::Type::INT32:
      case parquet::Type::FLOAT: width += 4; break;
      case parquet::Type::INT64:
      case parquet::Type::DOUBLE: width += 8; break;
      case parquet::Type::INT96: width += 12; break;
      case parquet::Type::FIXED_LEN_BYTE_ARRAY:
        width += static_cast<std::size_t>(column->type_length());
        break;
      default: exact = false; break; // BYTE_ARRAY: length unknown until read
    }
    ++width; // Arrow's per-column validity bitmap, rounded up to a byte a row
  }

  const std::size_t estimate = rows * width;
  const auto encoded = static_cast<std::size_t>(group.total_byte_size());
  return exact ? estimate : std::max(estimate, encoded);
}

/******************************************************************************/
std::shared_ptr<const Parquet::RowGroupBatches>
Parquet::Impl::row_group(int32_t index)
{
  if (shared_cache_) {
    // A memo in front of the shared cache.  A group holds thousands of
    // spectra and a reader walks a contiguous range, so nearly every call asks
    // for the group the previous one did -- and without the memo each such
    // call took the cache's one archive-wide mutex and copied a shared_future
    // every other thread in that group was copying too.  On a 347-group,
    // 717,924-spectrum archive that grew user CPU per spectrum 57% from 16 to
    // 128 threads, and voluntary context switches fivefold.
    //
    // Under THIS object's lock, not in plain members: e066604 memoised the
    // same thing unsynchronised and cbabd65 had to drop it, because one reader
    // shared across threads (one_shared_reader_across_threads_agrees) could be
    // handed another group's batches.  With one Spectra per thread, the
    // documented pattern, the lock is private to that thread and uncontended.
    //
    // HELD, not the weak_ptr e066604 used: a memo hit skips get(), so it does
    // not refresh the group's LRU recency, and a weak handle could lose the
    // group a reader is still walking and send it back to decode it again.
    //
    // ponytail: the memo keeps its group alive past the cache's budget -- one
    // group per Parquet that reads through the shared cache, i.e. per reader
    // thread (two with a separate peaks file) -- and the cache neither evicts
    // a held group nor counts it as room.  Bounded by the thread count, not
    // the budget: readers idle on distinct groups hold that many groups
    // whatever the budget says, as row_group_cache.h documents.
    // stats().held_bytes includes them.  Memo hits are not counted in
    // stats().hits; stats().decodes still counts every decode.
    std::size_t bytes = 0;
    {
      std::lock_guard<std::mutex> guard(cache_mutex_);
      if (index == memo_group_) return memo_batches_;
      bytes = group_bytes_locked_(index);
    }
    // Outside the lock: get() can wait on another reader's decode, and a
    // thread sharing this object must not be held up behind that.
    //
    // With what lets this reader decode AHEAD, on its own file handle, where
    // it would otherwise wait for another reader's decode of `index`.  The
    // callbacks take cache_mutex_ and decode_mutex_; neither is held here.
    RowGroupCache::Ahead ahead;
    ahead.groups = static_cast<int32_t>(group_bytes_.size());
    ahead.bytes = [this](int32_t g) {
      std::lock_guard<std::mutex> guard(cache_mutex_);
      return group_bytes_locked_(g);
    };
    ahead.decode = [this](int32_t g) { return decode_group_(g); };
    auto batches = shared_cache_->get(file_.file_name(), index, bytes,
                                      [this, index] { return decode_group_(index); },
                                      &ahead);
    std::lock_guard<std::mutex> guard(cache_mutex_);
    memo_group_ = index;
    memo_batches_ = batches;
    return batches;
  }

  // The lock is held across the DECODE as well as the bookkeeping.
  //
  // Guarding only the cache would leave the real hazard untouched: the decode
  // reads through one parquet::arrow::FileReader whose source is a single
  // seek-and-read handle (Util::ArrowFile_, and for an archive one zip_file_t
  // driven by zip_fseek).  Two threads decoding different row groups would race
  // on that file position, and Parquet pages carry no CRC by default, so the
  // result is silently wrong rather than an error.
  //
  // The cost is real but the right way round: concurrent readers serialise on
  // decode instead of corrupting each other, an uncontended lock costs nothing,
  // and a thread that waits finds the group already cached rather than decoding
  // it a second time.
  std::lock_guard<std::mutex> guard(cache_mutex_);

  for (const auto& [cached, batches] : cache_) {
    if (cached == index) return batches;
  }

  // Evict BEFORE decoding, not after.  Holding the outgoing group while the
  // incoming one is built makes the transient peak three groups where two
  // suffice -- about 21 MB for a 1,048,576-row group of three leaves -- for no
  // benefit: nothing reads the evicted entry between here and the insert below,
  // and anything still using it holds its own shared_ptr.
  if (cache_.size() >= kCachedGroups) cache_.erase(cache_.begin());

  auto batches = decode_group_(index);
  cache_.emplace_back(index, batches);
  return batches;
}

/******************************************************************************/
std::size_t Parquet::Impl::group_bytes_locked_(int32_t index)
{
  // Out of range: no slot, and RowGroup() below throws as it always did.
  const bool slot = index >= 0 && static_cast<std::size_t>(index) < group_bytes_.size();
  if (slot && group_bytes_[static_cast<std::size_t>(index)] != 0) {
    return group_bytes_[static_cast<std::size_t>(index)];
  }

  const auto metadata = reader_->parquet_reader()->metadata();
  const std::size_t bytes =
      decoded_row_group_bytes(*metadata->RowGroup(index), *metadata->schema());
  if (slot) group_bytes_[static_cast<std::size_t>(index)] = bytes;
  return bytes;
}

/******************************************************************************/
std::shared_ptr<const Parquet::RowGroupBatches>
Parquet::Impl::decode_group_(int32_t index)
{
  // Recursive: row_group() above holds cache_mutex_ on the private path and
  // does not on the shared path, and the file position needs the lock either
  // way.
  std::unique_lock<std::recursive_mutex> guard(decode_mutex_);

  auto batches = std::make_shared<Parquet::RowGroupBatches>();

  // Read the batches: without the declared-sorted key when that column proves
  // sorted (KeyRuns: ~40% of a point group, held as runs instead), in full
  // otherwise.  `strip` is null for a full read; false means the batches did
  // not come back in the shape the runs describe.
  auto read = [&](const KeyRuns* strip) -> bool {
    batches->batches.clear();
    std::vector<int> leaves;
    for (int c = 0; strip && c < reader_->parquet_reader()->metadata()->num_columns(); ++c) {
      if (c != strip->leaf) leaves.push_back(c);
    }
    auto reader_result = strip ? reader_->GetRecordBatchReader({index}, leaves)
                               : reader_->GetRecordBatchReader({index});
    if (!reader_result.ok()) {
      throw ParquetError("read row group " + std::to_string(index) + ": " +
                         reader_result.status().ToString());
    }
    int64_t rows = 0;
    for (auto maybe_batch : **reader_result) {
      if (!maybe_batch.ok()) {
        throw ParquetError("read row group " + std::to_string(index) + ": " +
                           maybe_batch.status().ToString());
      }
      std::shared_ptr<arrow::RecordBatch> batch = *maybe_batch;
      if (strip && !(batch = insert_key_placeholder(batch, *strip))) return false;
      rows += batch->num_rows();
      batches->batches.push_back(std::move(batch));
    }
    return strip == nullptr || rows == strip->start.back();
  };
  std::shared_ptr<KeyRuns> runs = scan_key_runs(*reader_, index);
  bool stripped = false;
  if (runs) {
    try {
      stripped = read(runs.get());
    } catch (const ParquetError&) {
      // The full read below meets the same fault, if it is one, and reports it.
    }
  }
  if (!stripped) {
    runs.reset();
    read(nullptr);
  }
  batches->key_runs = std::move(runs);

  // The leaf index of the column this file declares sorted, if any.
  int32_t sort_leaf = -1;
  {
    const auto rg = reader_->parquet_reader()->metadata()->RowGroup(index);
    for (const auto& cs : rg->sorting_columns()) {
      if (!cs.descending && !cs.nulls_first) sort_leaf = cs.column_idx;
      break;
    }
  }
  int64_t rows_seen = 0;
  std::vector<std::pair<int64_t, int64_t>> keyed;

  // Materialise Arrow's LAZY per-array boxing while this group is still
  // private to the decoding thread.
  //
  // RecordBatch::column(i) and StructArray::field(i) both box their child on
  // first access and cache it inside the array object.  A decoded group is
  // then shared, unsynchronised, by every reader thread through the cache, so
  // the first access from N threads at once is N concurrent lazy
  // initialisations of the same object.  Whether that is a data race depends
  // on the Arrow build; relying on it is not worth the crash it produces when
  // the answer is no.  Doing it here, once, under the decode lock, means every
  // reader afterwards only READS a pointer that is already set, and the cache
  // publishes the group through a future, which supplies the happens-before.
  //
  // The same walk collects the LEAVES in schema order, which is the order
  // Parquet numbers them in, so the file's declared sorting column can be
  // found by its leaf index.
  std::vector<std::shared_ptr<arrow::Array>> leaves;
  for (const auto& batch : batches->batches) {
    leaves.clear();
    for (int c = 0; c < batch->num_columns(); ++c) {
      std::shared_ptr<arrow::Array> col = batch->column(c);
      if (col && col->type_id() == arrow::Type::STRUCT) {
        const auto& sa = static_cast<const arrow::StructArray&>(*col);
        for (int f = 0; f < sa.num_fields(); ++f) leaves.push_back(sa.field(f));
      } else {
        leaves.push_back(std::move(col));
      }
    }
    batches->row_offset.push_back(rows_seen);
    rows_seen += batch->num_rows();

    // The sort-key span of this batch, when the file declares one and it is an
    // int64 leaf. A reader can then binary-search its way to the batch holding
    // a wanted key instead of walking and filtering every batch in the plan's
    // range -- which is what makes the skipped batches cost nothing at all.
    if (sort_leaf >= 0 && sort_leaf < static_cast<int32_t>(leaves.size()) &&
        batch->num_rows() > 0) {
      const auto& leaf = leaves[static_cast<std::size_t>(sort_leaf)];
      // NOT null_count(): on a freshly decoded array that is often
      // kUnknownNullCount (-1) until something forces the count, so comparing
      // it to zero silently disables the fast path. Ask whether a validity
      // bitmap exists at all instead -- that needs no counting.
      //
      // BOTH widths: an entity index is written unsigned, so this column
      // decodes as UInt64, and checking only for Int64 disabled the whole fast
      // path while looking like it worked.
      if (leaf && !leaf->data()->MayHaveNulls()) {
        if (leaf->type_id() == arrow::Type::INT64) {
          const auto& a = static_cast<const arrow::Int64Array&>(*leaf);
          keyed.push_back({a.Value(0), a.Value(a.length() - 1)});
          continue;
        }
        if (leaf->type_id() == arrow::Type::UINT64) {
          const auto& a = static_cast<const arrow::UInt64Array&>(*leaf);
          keyed.push_back({static_cast<int64_t>(a.Value(0)),
                           static_cast<int64_t>(a.Value(a.length() - 1))});
          continue;
        }
      }
    }
    keyed.clear();
    sort_leaf = -1; // one unusable batch disables the fast path for the group
  }
  if (sort_leaf >= 0 && keyed.size() == batches->batches.size()) {
    batches->key_leaf = sort_leaf;
    for (const auto& [f, l] : keyed) {
      batches->key_first.push_back(f);
      batches->key_last.push_back(l);
    }
  }
  return batches;
}

/******************************************************************************/
Parquet::Parquet(std::unique_ptr<IO::File> data,
                 Schema::File file,
                 std::shared_ptr<RowGroupCache> cache,
                 file_metadata_t footer)
    : impl_(std::make_unique<Impl>(
          std::move(data), std::move(file), std::move(cache), std::move(footer)))
{
}

/******************************************************************************/
Parquet::~Parquet() = default;

/******************************************************************************/
const Schema::File& Parquet::index_file() const { return impl_->file_; }

/******************************************************************************/
const std::shared_ptr<Schema::GroupMap>& Parquet::groups() const
{
  return impl_->groups_;
}

/******************************************************************************/
std::optional<Schema::Column> Parquet::field(std::string_view group_name,
                                             std::string_view field_name) const
{
  auto group_ptr = impl_->groups_->find(std::string{group_name});
  if (group_ptr == impl_->groups_->end()) return {};

  auto field_ptr = group_ptr->second->field(std::move(field_name));
  if (!field_ptr.has_value()) return {};

  return std::make_pair(group_ptr->second, field_ptr.value());
}

/******************************************************************************/
Parquet::file_metadata_t Parquet::file_metadata() const
{
  return impl_->reader_->parquet_reader()->metadata();
}

/******************************************************************************/
std::optional<std::string> Parquet::kv_string(const file_metadata_t& fmd,
                                              std::string_view key) const
{
  return get_kv_string(fmd, key);
}

/******************************************************************************/
std::optional<std::size_t> Parquet::kv_size_t(const file_metadata_t& fmd,
                                              std::string_view key) const
{
  return get_kv_uint(fmd, key);
}

/******************************************************************************/
parquet::arrow::FileReader& Parquet::reader() const { return *impl_->reader_; }

/******************************************************************************/
std::shared_ptr<const Parquet::RowGroupBatches> Parquet::row_group(int32_t index)
{
  return impl_->row_group(index);
}

/******************************************************************************/
bool Parquet::sorted_ascending(int32_t row_group, int32_t leaf_column) const
{
  return impl_->stats_->sorted_ascending(row_group, leaf_column);
}

/******************************************************************************/
Planner Parquet::planner(const Query& q)
{
  // NOT a cached page-index reader. Keeping one for the life of the file was
  // tried and is 48% SLOWER in the parallel phase (4.98 -> 7.36 s at 64
  // threads, 5.36 -> 7.67 at 192): a long-lived reader accumulates the parsed
  // index of every group it has touched, and that costs more than rebuilding
  // a fresh one per query. The constructor cost is real -- 79 of 95
  // thread-seconds of planning -- but this is not the way to remove it.
  return Planner(*impl_->reader_, q, impl_->stats_);
}

/******************************************************************************/
Executor Parquet::executor(const Projection& p) { return Executor(*this, p); }

} // namespace MzPeak::Util
