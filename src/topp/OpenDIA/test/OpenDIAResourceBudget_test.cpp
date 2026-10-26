// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIAResourceBudget.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <atomic>
#include <chrono>
#include <future>
#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

#ifdef OPENDIA_CLASS_TEST
#include <OpenMS/CONCEPT/ClassTest.h>
#define CHECK(...) TEST_TRUE((__VA_ARGS__))
#else
#define CHECK(...) do { if (!(__VA_ARGS__)) throw std::runtime_error("Check failed at line " + std::to_string(__LINE__)); } while (false);
#endif

using namespace OpenDIA;
using namespace std::chrono_literals;

namespace
{
  template <typename F> bool invalid(F operation)
  {
    try { operation(); }
    catch (const OpenMS::Exception::InvalidValue&) { return true; }
    return false;
  }
  // Synchronize on observable admission state, with a deadline rather than timing sleeps.
  bool waiters(ByteBudget& budget, ByteCount count)
  {
    const auto deadline = std::chrono::steady_clock::now() + 3s;
    while (budget.snapshot().waiters != count)
    {
      if (std::chrono::steady_clock::now() >= deadline) return false;
      std::this_thread::yield();
    }
    return true;
  }
}

#ifdef OPENDIA_CLASS_TEST
START_TEST(OpenDIAResourceBudget, "$Id$")
START_SECTION((byte accounting and cancellable admission))
#else
int main()
{
#endif
  static_assert(!std::is_copy_constructible_v<ByteBudget>);
  static_assert(!std::is_copy_constructible_v<ByteReservation>);
  static_assert(!std::is_copy_constructible_v<MemoryLease>);
  static_assert(std::is_nothrow_move_constructible_v<ByteReservation>);
  static_assert(std::is_nothrow_move_constructible_v<MemoryLease>);

  const auto maximum = std::numeric_limits<ByteCount>::max();
  CHECK(checkedByteAdd(maximum, 0) == maximum)
  CHECK(checkedByteAdd(maximum - 1, 1) == maximum)
  CHECK(invalid([&] { checkedByteAdd(maximum, 1); }))
  CHECK(checkedByteMultiply(maximum, 1) == maximum)
  CHECK(checkedByteMultiply(maximum, 0) == 0)
  CHECK(checkedByteMultiply(0, maximum) == 0)
  CHECK(invalid([&] { checkedByteMultiply(maximum, 2); }))
  CHECK(invalid([] { ByteBudget b(0, 0); }))
  CHECK(invalid([] { ByteBudget b(5, 5); }))
  CHECK(invalid([] { ByteBudget b(5, 6); }))
  CHECK(invalid([] { ByteBudget b(5, 0, 0); }))
  CHECK(invalid([] { ByteBudget b(5, 2, 4); }))
  CHECK(invalid([&] { ByteBudget b(1, maximum); }))
  ByteBudget smallest(1, 0, 1);
  CHECK(smallest.tryReserve(1).has_value())
  CHECK(smallest.snapshot().peak_total_bytes == 1)
  CHECK(invalid([&] { smallest.tryReserve(2); }))

  ByteBudget budget(100, 20, 10);
  auto initial = budget.snapshot();
  CHECK(initial.total_bytes == 100 && initial.native_reserve == 20 && initial.usable_bytes == 80)
  CHECK(initial.reserved_bytes == 0 && initial.committed_bytes == 0 && initial.waiters == 0 && !initial.cancelled)
  CHECK(invalid([&] { budget.tryReserve(9); }))
  CHECK(invalid([&] { budget.reserve(81); }))
  CHECK(invalid([&] { budget.reserve(0); }))
  auto job = budget.reserve(80);
  CHECK(job && job->bytes() == 80)
  CHECK(!budget.tryReserve(10))
  CHECK(invalid([&] { job->split(0); }))
  CHECK(invalid([&] { job->split(81); }))
  CHECK(invalid([&] { job->commit(81); }))
  CHECK(job->bytes() == 80 && budget.snapshot().reserved_bytes == 80 && budget.snapshot().committed_bytes == 0)
  auto tiny = job->split(1);
  auto tiny_lease = tiny.commit(1);
  CHECK(tiny.bytes() == 0 && tiny_lease.bytes() == 1 && job->bytes() == 79)
  auto lease = job->commit(59);
  auto counts = budget.snapshot();
  CHECK(counts.reserved_bytes == 0 && counts.committed_bytes == 60)
  CHECK(counts.peak_reserved_bytes == 80 && counts.peak_committed_bytes == 60 && counts.peak_total_bytes == 80)
  auto remainder = budget.tryReserve(20);
  CHECK(remainder && !budget.tryReserve(10))
  CHECK(remainder->commit(0).bytes() == 0)
  MemoryLease moved(std::move(lease));
  CHECK(lease.bytes() == 0 && moved.bytes() == 59)
  tiny_lease = std::move(moved); // Replaces and releases the old committed byte.
  CHECK(moved.bytes() == 0 && budget.snapshot().committed_bytes == 59)
  tiny_lease.reset();
  tiny_lease.reset();
  CHECK(budget.snapshot().committed_bytes == 0)
  auto first = budget.tryReserve(20);
  auto second = budget.tryReserve(30);
  ByteReservation pending(std::move(*first));
  CHECK(first->bytes() == 0 && pending.bytes() == 20)
  pending = std::move(*second);
  CHECK(second->bytes() == 0 && budget.snapshot().reserved_bytes == 30)
  auto entire = pending.split(30);
  CHECK(pending.bytes() == 0 && entire.bytes() == 30)
  entire.reset();
  ByteReservation empty;
  CHECK(empty.commit(0).bytes() == 0)
  CHECK(invalid([&] { empty.commit(1); }))
  CHECK(invalid([&] { empty.split(1); }))
  try { auto exception_job = budget.reserve(80); throw std::runtime_error("unwind"); }
  catch (const std::runtime_error&) {}
  CHECK(budget.snapshot().reserved_bytes == 0)
  try
  {
    auto exception_job = budget.reserve(80);
    auto exception_lease = exception_job->commit(80);
    throw std::runtime_error("committed unwind");
  }
  catch (const std::runtime_error&) {}
  CHECK(budget.snapshot().committed_bytes == 0)

  // Huge accounting values must not allocate proportional storage or overflow peaks.
  ByteBudget giant(maximum, 0);
  auto huge = giant.reserve(maximum);
  auto one = huge->split(1).commit(1);
  auto huge_lease = huge->commit(maximum - 1);
  CHECK(giant.snapshot().committed_bytes == maximum && giant.snapshot().peak_total_bytes == maximum)
  CHECK(!giant.tryReserve(1))
  one.reset();
  CHECK(giant.tryReserve(1).has_value())
  huge_lease.reset();
  ByteReservation surviving;
  MemoryLease surviving_lease;
  {
    ByteBudget temporary(10, 0);
    auto ticket = temporary.reserve(10);
    surviving = ticket->split(3);
    surviving_lease = ticket->commit(7);
  }
  CHECK(surviving.bytes() == 3 && surviving_lease.bytes() == 7)
  auto late_commit = surviving.commit(2);
  late_commit.reset();
  surviving_lease.reset();

  // Stop cancellation affects only the stopped waiter, not other jobs or tokens.
  ByteBudget blocking(10, 0);
  auto held = blocking.reserve(10);
  std::promise<bool> stopped_result, admitted_result;
  auto stopped_future = stopped_result.get_future();
  auto admitted_future = admitted_result.get_future();
  std::jthread stopped([&](std::stop_token stop) { stopped_result.set_value(!blocking.reserve(10, stop)); });
  std::jthread admitted([&](std::stop_token stop) { admitted_result.set_value(blocking.reserve(10, stop).has_value()); });
  const bool two_waiting = waiters(blocking, 2);
  CHECK(two_waiting)
  stopped.request_stop();
  const bool stop_ready = stopped_future.wait_for(3s) == std::future_status::ready;
  if (!stop_ready) blocking.cancel();
  CHECK(stop_ready)
  CHECK(stopped_future.get())
  CHECK(waiters(blocking, 1))
  CHECK(blocking.snapshot().reserved_bytes == 10 && !blocking.snapshot().cancelled)
  held->reset();
  const bool admit_ready = admitted_future.wait_for(3s) == std::future_status::ready;
  if (!admit_ready) blocking.cancel();
  CHECK(admit_ready)
  CHECK(admitted_future.get())
  stopped.join();
  admitted.join();
  CHECK(blocking.snapshot().reserved_bytes == 0 && blocking.snapshot().waiters == 0)
  std::stop_source pre_stopped;
  pre_stopped.request_stop();
  CHECK(!blocking.reserve(1, pre_stopped.get_token()))

  ByteBudget cancelled(10, 0);
  auto active = cancelled.reserve(10);
  std::vector<std::jthread> cancelled_threads;
  std::atomic<unsigned> rejected = 0;
  for (unsigned i = 0; i < 4; ++i)
  {
    cancelled_threads.emplace_back([&](std::stop_token stop) { if (!cancelled.reserve(1, stop)) ++rejected; });
  }
  CHECK(waiters(cancelled, 4))
  cancelled.cancel();
  const bool woke = waiters(cancelled, 0);
  CHECK(woke)
  cancelled_threads.clear();
  CHECK(rejected.load() == 4 && active->bytes() == 10)
  auto after_cancel = active->commit(10);
  CHECK(after_cancel.bytes() == 10 && cancelled.snapshot().committed_bytes == 10)
  after_cancel.reset();
  CHECK(!cancelled.tryReserve(1) && !cancelled.reserve(1))
  CHECK(invalid([&] { cancelled.reserve(11); }))

  ByteBudget contended(64, 8);
  std::atomic<bool> good = true;
  std::atomic<unsigned> completed = 0;
  std::vector<std::jthread> workers;
  for (unsigned i = 0; i < 8; ++i)
  {
    workers.emplace_back([&](std::stop_token stop)
    {
      for (unsigned j = 0; j < 200; ++j)
      {
        auto ticket = contended.reserve(14, stop);
        if (!ticket) { good = false; ++completed; return; }
        auto output = ticket->split(7).commit(6);
        auto scratch = ticket->commit(7);
        const auto s = contended.snapshot();
        if (s.reserved_bytes + s.committed_bytes > s.usable_bytes) good = false;
      }
      ++completed;
    });
  }
  const auto contention_deadline = std::chrono::steady_clock::now() + 5s;
  while (completed.load() != 8 && std::chrono::steady_clock::now() < contention_deadline)
  {
    std::this_thread::yield();
  }
  const bool contention_finished = completed.load() == 8;
  if (!contention_finished) contended.cancel();
  for (auto& worker : workers) worker.join();
  CHECK(contention_finished)
  CHECK(good.load())
  const auto final_counts = contended.snapshot();
  CHECK(final_counts.reserved_bytes == 0 && final_counts.committed_bytes == 0 && final_counts.waiters == 0)
  CHECK(final_counts.peak_total_bytes <= 56 && final_counts.peak_total_bytes >= 14)
#ifdef OPENDIA_CLASS_TEST
END_SECTION
END_TEST
#else
  return 0;
}
#endif
