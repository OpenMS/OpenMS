// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIAResourceBudget.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <utility>

namespace OpenDIA
{
  namespace
  {
    [[noreturn]] void invalid(const char* message)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "OpenDIA byte budget");
    }
  }

  ByteCount checkedByteAdd(ByteCount a, ByteCount b)
  {
    if (b > std::numeric_limits<ByteCount>::max() - a) { invalid("Byte addition overflow"); }
    return a + b;
  }

  ByteCount checkedByteMultiply(ByteCount a, ByteCount b)
  {
    if (a != 0 && b > std::numeric_limits<ByteCount>::max() / a) { invalid("Byte multiplication overflow"); }
    return a * b;
  }

  namespace Detail
  {
    struct BudgetState
    {
      std::mutex mutex;
      std::condition_variable_any changed;
      BudgetSnapshot counts;
      ByteCount minimum;

      void validate(ByteCount bytes) const
      {
        if (bytes < minimum || bytes > counts.usable_bytes) { invalid("Request outside byte budget bounds"); }
      }
      bool fits(ByteCount bytes) const
      {
        // Subtraction avoids overflow even for an admissible SIZE_MAX-sized job.
        return bytes <= counts.usable_bytes - counts.reserved_bytes - counts.committed_bytes;
      }
      void peaks()
      {
        counts.peak_reserved_bytes = std::max(counts.peak_reserved_bytes, counts.reserved_bytes);
        counts.peak_committed_bytes = std::max(counts.peak_committed_bytes, counts.committed_bytes);
        counts.peak_total_bytes = std::max(counts.peak_total_bytes, counts.reserved_bytes + counts.committed_bytes);
      }
    };
  }

  ByteBudget::ByteBudget(ByteCount total_bytes, ByteCount native_reserve, ByteCount minimum_request)
  {
    if (native_reserve >= total_bytes || minimum_request == 0 || minimum_request > total_bytes - native_reserve)
    {
      invalid("Invalid byte budget configuration");
    }
    state_ = std::make_shared<Detail::BudgetState>();
    state_->counts.total_bytes = total_bytes;
    state_->counts.native_reserve = native_reserve;
    state_->counts.usable_bytes = total_bytes - native_reserve;
    state_->minimum = minimum_request;
  }
  ByteBudget::~ByteBudget() = default;

  std::optional<ByteReservation> ByteBudget::tryReserve(ByteCount bytes)
  {
    std::lock_guard lock(state_->mutex);
    state_->validate(bytes);
    if (state_->counts.cancelled || !state_->fits(bytes)) { return std::nullopt; }
    state_->counts.reserved_bytes += bytes;
    state_->peaks();
    return ByteReservation(state_, bytes);
  }

  std::optional<ByteReservation> ByteBudget::reserve(ByteCount bytes, std::stop_token stop)
  {
    std::unique_lock lock(state_->mutex);
    state_->validate(bytes);
    if (state_->counts.cancelled || stop.stop_requested()) { return std::nullopt; }
    if (!state_->fits(bytes))
    {
      // Count only blocked admissions and restore on every exit, including wait exceptions.
      struct Waiter
      {
        ByteCount& count;
        explicit Waiter(ByteCount& value) : count(value) { ++count; }
        ~Waiter() { --count; }
      } waiter(state_->counts.waiters);
      state_->changed.wait(lock, stop, [&] { return state_->counts.cancelled || state_->fits(bytes); });
    }
    if (state_->counts.cancelled || stop.stop_requested()) { return std::nullopt; }
    state_->counts.reserved_bytes += bytes;
    state_->peaks();
    return ByteReservation(state_, bytes);
  }

  void ByteBudget::cancel()
  {
    { std::lock_guard lock(state_->mutex); state_->counts.cancelled = true; }
    state_->changed.notify_all();
  }
  BudgetSnapshot ByteBudget::snapshot() const
  {
    std::lock_guard lock(state_->mutex);
    return state_->counts;
  }

  ByteReservation::ByteReservation(std::shared_ptr<Detail::BudgetState> state, ByteCount bytes) noexcept : state_(std::move(state)), bytes_(bytes) {}
  ByteReservation::~ByteReservation() { reset(); }
  ByteReservation::ByteReservation(ByteReservation&& other) noexcept : state_(std::move(other.state_)), bytes_(std::exchange(other.bytes_, 0)) {}
  ByteReservation& ByteReservation::operator=(ByteReservation&& other) noexcept
  {
    if (this != &other) { reset(); state_ = std::move(other.state_); bytes_ = std::exchange(other.bytes_, 0); }
    return *this;
  }
  void ByteReservation::reset() noexcept
  {
    if (!state_) { return; }
    { std::lock_guard lock(state_->mutex); state_->counts.reserved_bytes -= bytes_; }
    state_->changed.notify_all();
    bytes_ = 0;
    state_.reset();
  }
  ByteReservation ByteReservation::split(ByteCount nonzero_bytes)
  {
    if (nonzero_bytes == 0 || nonzero_bytes > bytes_) { invalid("Invalid reservation split"); }
    ByteReservation result(state_, nonzero_bytes);
    bytes_ -= nonzero_bytes;
    if (bytes_ == 0) { state_.reset(); }
    return result;
  }
  MemoryLease ByteReservation::commit(ByteCount actual_bytes)
  {
    if (actual_bytes > bytes_) { invalid("Commit exceeds reservation"); }
    if (actual_bytes == 0) { reset(); return {}; }
    MemoryLease result(state_, actual_bytes);
    {
      std::lock_guard lock(state_->mutex);
      state_->counts.reserved_bytes -= bytes_;
      state_->counts.committed_bytes += actual_bytes;
      state_->peaks();
    }
    state_->changed.notify_all();
    bytes_ = 0;
    state_.reset();
    return result;
  }

  MemoryLease::MemoryLease(std::shared_ptr<Detail::BudgetState> state, ByteCount bytes) noexcept : state_(std::move(state)), bytes_(bytes) {}
  MemoryLease::~MemoryLease() { reset(); }
  MemoryLease::MemoryLease(MemoryLease&& other) noexcept : state_(std::move(other.state_)), bytes_(std::exchange(other.bytes_, 0)) {}
  MemoryLease& MemoryLease::operator=(MemoryLease&& other) noexcept
  {
    if (this != &other) { reset(); state_ = std::move(other.state_); bytes_ = std::exchange(other.bytes_, 0); }
    return *this;
  }
  void MemoryLease::reset() noexcept
  {
    if (!state_) { return; }
    { std::lock_guard lock(state_->mutex); state_->counts.committed_bytes -= bytes_; }
    state_->changed.notify_all();
    bytes_ = 0;
    state_.reset();
  }
}
