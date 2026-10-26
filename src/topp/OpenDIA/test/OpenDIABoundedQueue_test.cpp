// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIABoundedQueue.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <future>
#include <latch>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
// Checks remain active in Release/NDEBUG builds.
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif

using namespace OpenDIA;
namespace
{
  template <typename T>
  T ready(std::future<T>& result)
  {
    // Abort rather than let a blocked future destructor hide a failed wakeup.
    if (result.wait_for(std::chrono::seconds(5)) != std::future_status::ready) { std::abort(); }
    return result.get();
  }

  bool rejects(OpenMS::Size capacity)
  {
    try { BoundedQueue<int> queue(capacity); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }

  struct ReentrantDelete
  {
    BoundedQueue<std::unique_ptr<int, ReentrantDelete>>* queue{};
    int* destroyed{};
    void operator()(int* value) const
    {
      // cancel() must not hold its mutex while releasing payload ownership.
      (void)queue->snapshot();
      ++*destroyed;
      delete value;
    }
  };
}

#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIABoundedQueue, "$Id$")
START_SECTION((bounded queue ownership, lifecycle, stop and contention))
#else
int main()
{
#endif
  CHECK(rejects(0))
  CHECK(rejects(std::numeric_limits<OpenMS::Size>::max()))
  static_assert(!std::is_copy_constructible_v<BoundedQueue<int>>);
  static_assert(!std::is_move_constructible_v<BoundedQueue<int>>);
  {
    BoundedQueue<int> queue(3);
    CHECK(queue.push(0))
    CHECK(queue.push(1))
    CHECK(queue.push(2))
    CHECK(queue.snapshot().size == 3)
    for (int i = 0; i < 100; ++i)
    {
      CHECK(queue.pop() == i)
      CHECK(queue.push(int(i + 3)))
    }
    queue.close();
    CHECK(queue.pop() == 100)
    CHECK(queue.pop() == 101)
    CHECK(queue.pop() == 102)
    CHECK(!queue.pop())
    const auto state = queue.snapshot();
    CHECK(state.capacity == 3 && state.size == 0 && state.peak_size == 3)
    CHECK(state.pushed == 103 && state.popped == 103 && state.discarded == 0)
    CHECK(state.closed && !state.cancelled)
  }
  {
    BoundedQueue<std::unique_ptr<int>> queue(1);
    auto first = std::make_unique<int>(7);
    CHECK(queue.push(std::move(first)))
    CHECK(!first)
    queue.close();
    auto rejected = std::make_unique<int>(8);
    CHECK(!queue.push(std::move(rejected)))
    CHECK(rejected && *rejected == 8)
    auto item = queue.pop();
    CHECK(item && **item == 7)
    CHECK(!queue.pop())
  }
  {
    BoundedQueue<std::unique_ptr<int, ReentrantDelete>> queue(2);
    int destroyed = 0;
    CHECK(queue.push({new int(1), ReentrantDelete{&queue, &destroyed}}))
    CHECK(queue.push({new int(2), ReentrantDelete{&queue, &destroyed}}))
    queue.cancel();
    queue.cancel();
    queue.close();
    const auto state = queue.snapshot();
    CHECK(destroyed == 2 && state.size == 0 && state.discarded == 2)
    CHECK(state.closed && state.cancelled && state.pushed == 2 && state.popped == 0)
    CHECK(!queue.pop())
    auto rejected = std::unique_ptr<int, ReentrantDelete>(new int(3), ReentrantDelete{&queue, &destroyed});
    CHECK(!queue.push(std::move(rejected)))
    CHECK(rejected && *rejected == 3 && destroyed == 2)
  }
  // A pre-requested stop takes precedence even when the queue is ready.
  {
    BoundedQueue<std::unique_ptr<int>> queue(2);
    std::stop_source stop;
    stop.request_stop();
    auto item = std::make_unique<int>(9);
    CHECK(!queue.push(std::move(item), stop.get_token()))
    CHECK(item && *item == 9)
    CHECK(queue.push(std::move(item)))
    CHECK(!queue.pop(stop.get_token()))
    CHECK(queue.snapshot().size == 1 && !queue.snapshot().cancelled)
    CHECK(**queue.pop() == 9)
  }
  // Full and empty predicates force waiting. Latches establish arrival without sleeps.
  for (int operation = 0; operation < 3; ++operation)
  {
    BoundedQueue<std::unique_ptr<int>> full(1);
    BoundedQueue<int> empty(1);
    CHECK(full.push(std::make_unique<int>(1)))
    std::stop_source stop;
    std::latch entered(2);
    auto push = std::async(std::launch::async, [&]
    {
      auto item = std::make_unique<int>(2);
      entered.count_down();
      const bool admitted = full.push(std::move(item), stop.get_token());
      return !admitted && item && *item == 2;
    });
    auto pop = std::async(std::launch::async, [&]
    {
      entered.count_down();
      return !empty.pop(stop.get_token());
    });
    entered.wait();
    CHECK(push.wait_for(std::chrono::seconds(0)) == std::future_status::timeout)
    CHECK(pop.wait_for(std::chrono::seconds(0)) == std::future_status::timeout)
    if (operation == 0) { full.close(); empty.close(); }
    else if (operation == 1) { full.cancel(); empty.cancel(); }
    else { stop.request_stop(); }
    CHECK(ready(push))
    CHECK(ready(pop))
    if (operation != 1)
    {
      CHECK(**full.pop() == 1)
      CHECK(full.snapshot().discarded == 0)
    }
    if (operation == 2)
    {
      CHECK(!full.snapshot().closed && !empty.snapshot().cancelled)
      CHECK(empty.push(3))
      CHECK(empty.pop() == 3)
    }
  }
  {
    // Stall the consumer behind an explicit gate while producers contend on full.
    BoundedQueue<int> queue(2);
    CHECK(queue.push(-2))
    CHECK(queue.push(-1))
    std::latch entered(4);
    std::promise<void> drain;
    auto gate = drain.get_future();
    std::array<std::future<bool>, 4> producers;
    for (int i = 0; i < 4; ++i)
    {
      producers[i] = std::async(std::launch::async, [&, i]
      {
        entered.count_down();
        return queue.push(int(i));
      });
    }
    auto consumer = std::async(std::launch::async, [&]
    {
      gate.wait();
      std::array<int, 6> values{};
      for (auto& value : values) { value = *queue.pop(); }
      return values;
    });
    entered.wait();
    for (auto& producer : producers)
    {
      CHECK(producer.wait_for(std::chrono::seconds(0)) == std::future_status::timeout)
    }
    CHECK(queue.snapshot().size == 2 && queue.snapshot().peak_size == 2)
    drain.set_value();
    for (auto& producer : producers) { CHECK(ready(producer)) }
    auto values = ready(consumer);
    CHECK(values[0] == -2 && values[1] == -1)
    std::sort(values.begin() + 2, values.end());
    CHECK((values == std::array<int, 6>{-2, -1, 0, 1, 2, 3}))
    const auto state = queue.snapshot();
    CHECK(state.size == 0 && state.peak_size == 2 && state.pushed == 6 && state.popped == 6)
  }
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
