/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <condition_variable>
#include <functional>
#include <future>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace arrow {
class RecordBatch;
}

namespace MzPeak::Util {

struct KeyRuns;

/// A decoded row group: the record batches Parquet hands back for it, plus
/// what a reader needs to find a row in them WITHOUT touching Arrow.
///
/// The extra vectors are filled once, by the thread that decoded the group,
/// and then only read.  That matters twice over: Arrow's accessors
/// (RecordBatch::column, StructArray::field) are not free even when the child
/// is already boxed -- libstdc++ implements the atomic shared_ptr load they
/// use with a pool of 16 process-wide mutexes -- and a reader that can pick
/// its batch arithmetically never calls them for the batches it skips.
struct RowGroupBatches {
  std::vector<std::shared_ptr<arrow::RecordBatch>> batches;

  /// First row index of each batch, relative to the row group.
  std::vector<int64_t> row_offset;

  /// First and last value of the column this file declares sorted, per batch.
  /// Empty unless the group declares exactly one sorted column and it decodes
  /// to int64 -- the only case the reader's fast path uses.
  std::vector<int64_t> key_first;
  std::vector<int64_t> key_last;

  /// Which LEAF column key_first/key_last describe.  A reader must check this
  /// against the column it is querying: the spans belong to the column the
  /// FILE declares sorted, which is not necessarily the one being asked for,
  /// and searching one column's spans with another's value silently returns
  /// the wrong rows.
  int32_t key_leaf = -1;

  /// Set when the group's declared-sorted key column was verified sorted and
  /// NOT kept: its values are these runs, and its place in every batch holds a
  /// zero-byte null placeholder.  A reader must check this before touching
  /// that column.  See key_runs.h.
  std::shared_ptr<const KeyRuns> key_runs;

  bool has_keys() const { return !key_first.empty(); }
};

/******************************************************************************/
/// One decoded copy of each row group, shared by every reader over an archive.
///
/// A reader per thread is the library's concurrency model, and each reader
/// used to keep its own two decoded groups.  With sixteen readers over a
/// seven-group file that is up to thirty-two decoded copies of seven groups:
/// a gigabyte of duplicates.  This cache lives on the archive's Manager, so
/// the readers share it, and a group is decoded ONCE however many readers
/// want it.  The second reader that asks while the first is still decoding
/// waits for that decode rather than starting its own.
///
/// Decoding itself is not serialised here: every Parquet object has its own
/// file handle, so two readers decoding two different groups run in parallel.
/// The lock covers only the bookkeeping.
///
/// Bounded by BYTES in two places, because eviction alone cannot hold a
/// budget: a group still being decoded is not evictable, so N readers
/// decoding at once used to pin N groups however small the budget was.
///
/// So a decode is also ADMITTED: a reader that misses waits until the bytes
/// already in flight leave room for its group (half the budget is reserved
/// for them, the rest stays available to hold decoded groups).  A single
/// reader is always admitted, so a group larger than the budget still makes
/// progress.
///
/// What the budget does NOT bound is the groups readers hold.  A group a
/// reader still holds is never evicted (dropping it would free nothing), nor
/// counted as room to decode ahead into.  And every Parquet that reads
/// through this cache holds the last group it read -- its memo, see
/// parquet.cpp -- until it reads another group or is destroyed, idle or not.
/// So peak decoded memory is about the budget PLUS one group per live
/// Parquet reading through the cache: one per reader, two for a reader with
/// a separate peaks file.  Readers on distinct groups hold that many groups
/// whatever the budget, so size memory as budget + readers x the largest
/// group (x2 with a separate peaks file); stats().held_bytes includes them.
/// They also take room decoding ahead needs (see get()): once the groups
/// readers hold fill the budget, decoding ahead stops, so size the budget
/// above readers x the largest group to keep it.
/// A caller that keeps a slice of a group (a Chromatogram does, through its
/// decoder) keeps that group's buffers alive on top, evicted or not.
///
/// A budget smaller than the working set thrashes (a group is decoded again
/// each time a reader comes back to it) but stays correct.
///
/// A reader that finds its group IN FLIGHT does not simply wait for it when it
/// can decode AHEAD instead: it claims the next group of the same file that
/// nobody has claimed, decodes that, and only then comes back for its own.
/// See get().  A group decoded that way is evicted like any other, least
/// recently used first.  And a reader does not decode ahead while its group
/// is one the cache evicted lately and is decoding again: the budget is
/// then too small for what readers come back to, and decoding ahead would
/// evict more of it.
class RowGroupCache final {
public:
  using Batches = std::shared_ptr<const RowGroupBatches>;
  using Decode = std::function<Batches()>;

  /// What a reader needs to decode groups OTHER than the one it asked for:
  /// the file's group count, and each group's budget size and decode, on the
  /// reader's own file handle.  Parquet hands one to every get() it makes,
  /// whichever way its reader moves through the file; get() calls these
  /// under the cache's lock (bytes) and outside it (decode).
  struct Ahead {
    int32_t groups = 0;
    std::function<std::size_t(int32_t)> bytes;
    std::function<Batches(int32_t)> decode;
  };

  /// 1 GiB: room for ~30 groups of a typical archive, two of a chunked Astral
  /// one.  Consumers that know their thread count set their own.
  static constexpr std::size_t kDefaultBudget = std::size_t(1) << 30;

  explicit RowGroupCache(std::size_t budget_bytes = kDefaultBudget);

  /// The decoded group @p group of @p file: cached, in flight (wait for it),
  /// or decoded now through @p decode.  @p bytes is the group's size for the
  /// budget; Parquet's uncompressed row-group size is the right estimate.
  /// Rethrows whatever @p decode threw, for waiters too, and forgets the
  /// entry so a later call tries again.
  ///
  /// With @p ahead, a call that would wait for another reader's decode of
  /// @p group first decodes later groups of @p file that nobody has claimed,
  /// within the budget -- unless that decode brings back a group evicted
  /// lately.  A failed decode AHEAD is not this call's error: it reaches that
  /// group's own waiters and is forgotten like any other.
  Batches get(const std::string& file,
              int32_t group,
              std::size_t bytes,
              const Decode& decode,
              const Ahead* ahead = nullptr);

  /// The bytes of decoded groups eviction holds the cache to, and the room
  /// decodes are admitted into (half of it in flight at once).  Not a cap on
  /// decoded memory: groups readers hold are never evicted, so they take it
  /// past the budget, by about one group per live Parquet reading through
  /// the cache (see the class comment).  Lowering it evicts at once, down to
  /// what readers hold.
  void set_budget(std::size_t bytes);
  std::size_t budget() const;

  struct Stats {
    std::size_t decodes = 0;   ///< groups decoded (each miss, once)
    std::size_t hits = 0;      ///< served from a decoded group
    std::size_t waits = 0;     ///< served by waiting for another reader's decode
    std::size_t admission_waits = 0; ///< decodes held back to stay in budget
    std::size_t ahead_decodes = 0; ///< of decodes, those done AHEAD by a reader that would have waited
    std::size_t evictions = 0;
    std::size_t held_bytes = 0; ///< budget currently accounted for
    /// Thread-nanoseconds spent decoding, waiting for another reader's decode,
    /// and waiting for admission.  Timed per cache call, never per row.
    std::uint64_t decode_ns = 0;
    std::uint64_t wait_ns = 0;
    std::uint64_t admission_ns = 0;
  };
  Stats stats() const;

private:
  using Key = std::pair<std::string, int32_t>;
  struct Entry {
    std::shared_future<Batches> future;
    std::size_t bytes = 0;
    std::uint64_t used = 0;
    bool ready = false;
    bool ahead = false; ///< decoded ahead, and no reader has asked for it yet
    bool redecode = false; ///< claimed within lately_locked_() of its eviction; see get()
  };

  /// How many ticks of the clock after its eviction a group claimed again
  /// counts as a re-decode: a few turns of the clock over every entry.  A
  /// working set that fits the budget, read at random, asks for each of its
  /// groups again within about one turn; a walk asks again for a group it
  /// left behind only on its next pass, after the rest of the file.  Caller
  /// holds mutex_.
  std::uint64_t lately_locked_() const { return 4 * entries_.size() + 16; }

  /// Is a READY entry still held by a reader (its per-reader memo, or a read
  /// in progress)?  Evicting such an entry frees nothing -- the reader keeps
  /// the batches alive -- and sends the next reader that asks back to decode
  /// it again.
  static bool held_by_reader_(const Entry& e) { return e.future.get().use_count() > 1; }

  /// Claim @p key for a decode by the caller.  Caller holds mutex_.
  std::promise<Batches> claim_locked_(const Key& key, std::size_t bytes, bool ahead);

  /// Run @p decode for a claimed @p key outside the lock and publish it (or
  /// its exception) to every waiter.
  Batches decode_claimed_(const Key& key,
                          std::size_t bytes,
                          std::promise<Batches>& promise,
                          const Decode& decode);

  /// Decode groups after @p group of @p file that nobody has claimed, while
  /// the budget has room for them and @p group is still in flight.  Returns
  /// with mutex_ held, as it was called.
  void decode_ahead_locked_(std::unique_lock<std::mutex>& lock,
                            const std::string& file,
                            int32_t group,
                            const Ahead& ahead);

  /// Bytes a decode ahead must not count as room: groups in flight and groups
  /// a reader holds, which eviction cannot reclaim, and every group decoded
  /// ahead that nobody has asked for yet, from this burst or an earlier one,
  /// until it is asked for or evicted -- eviction could reclaim those, but a
  /// burst must not push them out to go on.  Caller holds mutex_.
  std::size_t resident_locked_() const;

  /// Evict least-recently-used READY entries, never @p keep and never one a
  /// reader still holds, until the budget holds or nothing evictable is left.
  /// Groups decoded ahead are ordered by the same clock as the rest.  Each
  /// eviction is recorded in evicted_ for claim_locked_.  Caller holds mutex_.
  void evict_locked_(const Key& keep);

  mutable std::mutex mutex_;
  std::condition_variable admit_;
  std::map<Key, Entry> entries_;
  std::map<Key, std::uint64_t> evicted_; ///< tick_ at each group's last eviction, lately
  std::size_t budget_;
  std::size_t held_ = 0;
  std::size_t in_flight_ = 0; ///< bytes of entries still decoding
  std::uint64_t tick_ = 0;
  Stats stats_;
};

} // namespace MzPeak::Util
