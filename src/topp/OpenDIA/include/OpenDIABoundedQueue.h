// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <optional>
#include <stop_token>
#include <type_traits>
#include <utility>
#include <vector>

namespace OpenDIA
{
  /// Coherent queue accounting. Capacity counts items, not payload bytes.
  struct QueueSnapshot
  {
    OpenMS::Size capacity{}, size{}, peak_size{}, pushed{}, popped{}, discarded{};
    bool closed{}, cancelled{};
  };

  /**
    @brief Private fixed-capacity queue with graceful close and immediate cancellation.

    Callers must separately use a byte budget. The fixed acquisition protocol is
    caller responsibility: reserve minimal same-budget drain space BEFORE starting
    producers, reserve the whole peak job before materialization, and make no more
    blocking admissions while holding that job. The consumer reuses the drain.
    Destruction requires that no threads are using this object.
  */
  template <typename T>
  class BoundedQueue
  {
    static_assert(std::is_nothrow_move_constructible_v<T>, "Queue payload must be nothrow move constructible");

  public:
    explicit BoundedQueue(OpenMS::Size capacity) : ring_(checkedCapacity_(capacity))
    {
      state_.capacity = capacity;
    }

    BoundedQueue(const BoundedQueue&) = delete;
    BoundedQueue& operator=(const BoundedQueue&) = delete;
    BoundedQueue(BoundedQueue&&) = delete;
    BoundedQueue& operator=(BoundedQueue&&) = delete;

    /// Moves value only after successful admission. Stop affects this call only.
    bool push(T&& value, std::stop_token stop = {})
    {
      std::unique_lock lock(mutex_);
      if (!changed_.wait(lock, stop, [&] { return state_.closed || state_.size < state_.capacity; }) ||
          stop.stop_requested() || state_.closed)
      {
        return false;
      }
      ring_[tail_].emplace(std::move(value));
      tail_ = next_(tail_);
      ++state_.size;
      ++state_.pushed;
      state_.peak_size = std::max(state_.peak_size, state_.size);
      lock.unlock();
      changed_.notify_all();
      return true;
    }

    /// Returns EOF after graceful drain, cancellation, or this call's stop.
    std::optional<T> pop(std::stop_token stop = {})
    {
      std::unique_lock lock(mutex_);
      if (!changed_.wait(lock, stop, [&] { return state_.closed || state_.size != 0; }) ||
          stop.stop_requested() || state_.cancelled || state_.size == 0)
      {
        return std::nullopt;
      }
      std::optional<T> result(std::in_place, std::move(*ring_[head_]));
      ring_[head_].reset();
      head_ = next_(head_);
      --state_.size;
      ++state_.popped;
      lock.unlock();
      changed_.notify_all();
      return result;
    }

    /// Reject further pushes but preserve queued items for draining.
    void close()
    {
      {
        std::lock_guard lock(mutex_);
        state_.closed = true;
      }
      changed_.notify_all();
    }

    /// Discard queued items outside the mutex, allowing reentrant destructors.
    void cancel()
    {
      std::vector<std::optional<T>> discarded;
      {
        std::lock_guard lock(mutex_);
        state_.closed = true;
        state_.cancelled = true;
        state_.discarded += state_.size;
        state_.size = 0;
        ring_.swap(discarded);
      }
      changed_.notify_all();
    }

    QueueSnapshot snapshot() const
    {
      std::lock_guard lock(mutex_);
      return state_;
    }

  private:
    static OpenMS::Size checkedCapacity_(OpenMS::Size capacity)
    {
      if (capacity == 0 || capacity > std::numeric_limits<OpenMS::Size>::max() / sizeof(std::optional<T>) ||
          capacity > std::vector<std::optional<T>>().max_size())
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "Invalid bounded queue capacity", "OpenDIA");
      }
      return capacity;
    }

    OpenMS::Size next_(OpenMS::Size index) const { return index == state_.capacity - 1 ? 0 : index + 1; }

    mutable std::mutex mutex_;
    std::condition_variable_any changed_;
    std::vector<std::optional<T>> ring_;
    QueueSnapshot state_;
    OpenMS::Size head_{}, tail_{};
  };
}
