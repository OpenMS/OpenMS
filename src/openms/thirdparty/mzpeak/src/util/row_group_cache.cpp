/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/row_group_cache.h"

#include <algorithm>
#include <chrono>
#include <exception>

namespace MzPeak::Util {

namespace
{
  std::uint64_t ns_since(std::chrono::steady_clock::time_point t0)
  {
    return static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                          std::chrono::steady_clock::now() - t0)
                                          .count());
  }

  /// How far past its own group a waiting reader may look for one to decode.
  ///
  /// ponytail: a fixed window, not a tuned one.  The budget is what really
  /// bounds decoding ahead (see decode_ahead_locked_); this only bounds the
  /// scan past groups other readers have already claimed.
  constexpr int32_t kAheadWindow = 64;
} // namespace

/******************************************************************************/
RowGroupCache::RowGroupCache(std::size_t budget_bytes) : budget_(budget_bytes) {}

/******************************************************************************/
RowGroupCache::Batches RowGroupCache::get(const std::string& file,
                                          int32_t group,
                                          std::size_t bytes,
                                          const Decode& decode,
                                          const Ahead* ahead)
{
  const Key key(file, group);
  bool decoded_ahead = false;

  std::unique_lock<std::mutex> lock(mutex_);
  // Both on arrival and again after an admission wait or decoding ahead:
  // another reader may have claimed, finished or dropped this very group
  // meanwhile.
  for (;;) {
    if (auto it = entries_.find(key); it != entries_.end()) {
      it->second.used = ++tick_;
      it->second.ahead = false; // asked for now

      // In flight on another reader.  Waiting for it is how readers that walk
      // a file together lose their time: they all arrive at the group the
      // first of them is decoding, and every one of them then idles for that
      // decode.  Measured on a 347-group, 717,924-spectrum archive at 128
      // threads: 4,013 waits, 112 of 271 thread-seconds, and a read path that
      // took ~2 s at ANY thread count from 16 up, because only the readers
      // that happened to arrive first ever decoded.
      //
      // So decode the next groups nobody has claimed instead, which a reader
      // walking forward needs next anyway, and come back for this one.  The
      // waiting threads become the decoders the walk was short of.
      //
      // Only WHILE waiting.  Also keeping a lead -- the first reader of a group
      // decoded ahead decoding one more -- was tried and is worse (128 threads:
      // 581 decodes of 347 groups, waits 809 -> 1,763): a lead the cache
      // cannot see the walk's working set behind pushes out groups that
      // readers between two chunks do not hold at that moment but still need.
      //
      // And not while this decode brings back a group evicted lately
      // (Entry::redecode).  The budget is then already losing groups readers
      // come back to, and each group decoded ahead would push out another:
      // readers that do not walk forward -- random access, say -- then collide
      // on the groups they lost, decode ahead again and lose more.  Measured
      // with 32 threads reading 24 scattered groups at random under a
      // 32-group budget, 2 ms decodes: 215-5,546 decodes without this check,
      // 53-58 with it, where 24 would do.
      if (!it->second.ready && !it->second.redecode && ahead && !decoded_ahead) {
        decoded_ahead = true;
        decode_ahead_locked_(lock, file, group, *ahead);
        continue;
      }

      const bool ready = it->second.ready;
      if (ready) {
        ++stats_.hits;
      } else {
        ++stats_.waits;
      }
      std::shared_future<Batches> future = it->second.future;
      lock.unlock();
      // Outside the lock: a wait here must not stop other readers decoding
      // OTHER groups meanwhile.
      if (ready) return future.get();
      const auto t0 = std::chrono::steady_clock::now();
      Batches got = future.get();
      const std::uint64_t ns = ns_since(t0);
      std::lock_guard<std::mutex> guard(mutex_);
      stats_.wait_ns += ns;
      return got;
    }

    // Miss.  Admission control, and the reason the bytes in flight follow
    // the budget rather than the reader count: an entry being decoded is not
    // evictable, so without this N readers pin N groups whatever the budget
    // says.  (Groups readers HOLD are not bounded by it; see the class
    // comment.)  Half the budget is reserved for decodes in flight; the rest
    // stays available to HOLD decoded groups, which is what stops a fully
    // admitted cache from evicting everything it just decoded.
    //
    // The in_flight_ == 0 escape guarantees progress: a group bigger than
    // the budget is still decoded, alone.
    if (in_flight_ == 0 || in_flight_ + bytes <= budget_ / 2) break;
    ++stats_.admission_waits;
    const auto t0 = std::chrono::steady_clock::now();
    admit_.wait(lock);
    stats_.admission_ns += ns_since(t0);
  }

  // Claim the slot before decoding, so a second reader arriving during the
  // decode finds it in flight and waits instead of decoding it again.
  std::promise<Batches> promise = claim_locked_(key, bytes, false);
  lock.unlock();
  return decode_claimed_(key, bytes, promise, decode);
}

/******************************************************************************/
void RowGroupCache::decode_ahead_locked_(std::unique_lock<std::mutex>& lock,
                                         const std::string& file,
                                         int32_t group,
                                         const Ahead& ahead)
{
  const Key wanted(file, group);
  const int32_t last = std::min(ahead.groups - 1, group + kAheadWindow);

  for (int32_t g = group + 1; g <= last; ++g) {
    // Only while there is nothing better to do: once the group this reader
    // came for is ready (or its decode failed and it is gone), it has work.
    if (auto it = entries_.find(wanted); it == entries_.end() || it->second.ready) return;

    const Key key(file, g);
    if (entries_.contains(key)) continue; // claimed by someone else

    // Within the budget, both halves of it: the same in-flight room a miss is
    // admitted into (without the progress escape -- nothing waits on this
    // decode), and room to KEEP the group next to everything resident_locked_
    // counts, groups decoded ahead and not yet asked for included.  Eviction
    // takes those in plain LRU order, so without counting them a burst would
    // go on by evicting the groups it decoded first: work thrown away, and a
    // wait moved to the reader that wanted them.  It counts every one of them,
    // an earlier burst's too, until it is asked for or evicted; those leave
    // this burst less room, which errs toward decoding less ahead.
    const std::size_t bytes = ahead.bytes(g);
    if (in_flight_ + bytes > budget_ / 2) return;
    if (resident_locked_() + bytes > budget_) return;

    std::promise<Batches> promise = claim_locked_(key, bytes, true);
    lock.unlock();
    try {
      (void)decode_claimed_(key, bytes, promise, [&] { return ahead.decode(g); });
    } catch (...) {
      // Not this reader's error: the group's own waiters got it, and the
      // entry is gone, so whoever needs the group next tries again.
    }
    lock.lock();
  }
}

/******************************************************************************/
std::promise<RowGroupCache::Batches>
RowGroupCache::claim_locked_(const Key& key, std::size_t bytes, bool ahead)
{
  std::promise<Batches> promise;
  ++stats_.decodes;
  if (ahead) ++stats_.ahead_decodes;
  bool redecode = false;
  if (auto it = evicted_.find(key); it != evicted_.end()) {
    redecode = tick_ - it->second <= lately_locked_();
    evicted_.erase(it);
  }
  entries_.emplace(key,
                   Entry{promise.get_future().share(), bytes, ++tick_, false, ahead, redecode});
  held_ += bytes;
  in_flight_ += bytes;
  evict_locked_(key);
  return promise;
}

/******************************************************************************/
RowGroupCache::Batches RowGroupCache::decode_claimed_(const Key& key,
                                                      std::size_t bytes,
                                                      std::promise<Batches>& promise,
                                                      const Decode& decode)
{
  Batches batches;
  const auto t0 = std::chrono::steady_clock::now();
  try {
    batches = decode();
  } catch (...) {
    {
      std::lock_guard<std::mutex> guard(mutex_);
      in_flight_ -= bytes;
      if (auto it = entries_.find(key); it != entries_.end()) {
        held_ -= it->second.bytes;
        entries_.erase(it);
      }
    }
    admit_.notify_all();
    // Waiters holding the future get the same exception.
    promise.set_exception(std::current_exception());
    throw;
  }

  // The value is published BEFORE the entry is marked ready: held_by_reader_()
  // reads a ready entry's value under the lock, and must never block there.
  promise.set_value(batches);
  {
    std::lock_guard<std::mutex> guard(mutex_);
    stats_.decode_ns += ns_since(t0);
    in_flight_ -= bytes;
    if (auto it = entries_.find(key); it != entries_.end()) it->second.ready = true;
    // Now evictable: the entry that could not be reclaimed while it decoded.
    evict_locked_(key);
  }
  admit_.notify_all();
  return batches;
}

/******************************************************************************/
std::size_t RowGroupCache::resident_locked_() const
{
  std::size_t bytes = 0;
  for (const auto& [key, e] : entries_) {
    if (!e.ready || e.ahead || held_by_reader_(e)) bytes += e.bytes;
  }
  return bytes;
}

/******************************************************************************/
void RowGroupCache::evict_locked_(const Key& keep)
{
  while (held_ > budget_) {
    auto victim = entries_.end();
    for (auto it = entries_.begin(); it != entries_.end(); ++it) {
      const Entry& e = it->second;
      if (!e.ready || it->first == keep || held_by_reader_(e)) continue;
      // Oldest first, groups decoded ahead included.  One is claimed with a
      // fresh tick, so it outlives every group used before it was decoded and
      // then ages out like any other.  Keeping them until nothing else was
      // left let one collision off a reader that does not walk forward fill
      // the budget with groups nobody asks for, and every later miss then
      // evicted a group of the working set instead.
      if (victim == entries_.end() || e.used < victim->second.used) victim = it;
    }
    if (victim == entries_.end()) return; // everything left is in use; overshoot
    // Remembered while it can still mark a re-decode (claim_locked_).  Once
    // there are more records than lately_locked_() ticks, the older ones go,
    // so they stay a few times the entries in number.
    if (const std::uint64_t lately = lately_locked_(); evicted_.size() > lately) {
      std::erase_if(evicted_, [&](const auto& r) { return tick_ - r.second > lately; });
    }
    evicted_[victim->first] = tick_;
    held_ -= victim->second.bytes;
    entries_.erase(victim);
    ++stats_.evictions;
  }
}

/******************************************************************************/
void RowGroupCache::set_budget(std::size_t bytes)
{
  std::lock_guard<std::mutex> guard(mutex_);
  budget_ = bytes;
  evict_locked_(Key());
  admit_.notify_all(); // a raised budget may admit a waiting decode
}

std::size_t RowGroupCache::budget() const
{
  std::lock_guard<std::mutex> guard(mutex_);
  return budget_;
}

RowGroupCache::Stats RowGroupCache::stats() const
{
  std::lock_guard<std::mutex> guard(mutex_);
  Stats s = stats_;
  s.held_bytes = held_;
  return s;
}

} // namespace MzPeak::Util
