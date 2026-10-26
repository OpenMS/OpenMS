// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <memory>
#include <optional>
#include <stop_token>

namespace OpenDIA
{
  using ByteCount = OpenMS::Size;

  /// Coherent accounting snapshot. Reserved and committed bytes are disjoint.
  struct BudgetSnapshot
  {
    ByteCount total_bytes = 0, native_reserve = 0, usable_bytes = 0;
    ByteCount reserved_bytes = 0, committed_bytes = 0;
    ByteCount peak_reserved_bytes = 0, peak_committed_bytes = 0, peak_total_bytes = 0;
    ByteCount waiters = 0;
    bool cancelled = false;
  };

  /** @brief Checked byte addition, throwing InvalidValue on overflow.
    @param[in] a First byte count.
    @param[in] b Second byte count.
  */
  ByteCount checkedByteAdd(ByteCount a, ByteCount b);
  /** @brief Checked byte multiplication, throwing InvalidValue on overflow.
    @param[in] a Element count.
    @param[in] b Bytes per element.
  */
  ByteCount checkedByteMultiply(ByteCount a, ByteCount b);

  namespace Detail { struct BudgetState; }
  class ByteReservation;

  /**
    @brief Unique accounting ownership for a committed allocation.

    The shared accounting state outlives its ByteBudget if necessary. Declare
    this token before the allocations it accounts for, so it is destroyed last.
    Moving a token transfers ownership without changing the accounting total.
  */
  class MemoryLease
  {
  public:
    MemoryLease() noexcept = default;
    ~MemoryLease();
    MemoryLease(const MemoryLease&) = delete;
    MemoryLease& operator=(const MemoryLease&) = delete;
    MemoryLease(MemoryLease&& other) noexcept;
    MemoryLease& operator=(MemoryLease&& other) noexcept;
    ByteCount bytes() const noexcept { return bytes_; }
    /// Release this accounting only after its associated allocations are gone.
    void reset() noexcept;

  private:
    friend class ByteReservation;
    MemoryLease(std::shared_ptr<Detail::BudgetState> state, ByteCount bytes) noexcept;
    std::shared_ptr<Detail::BudgetState> state_;
    ByteCount bytes_ = 0;
  };

  /**
    @brief Pending peak-job bytes, released automatically on destruction.

    Splits transfer ownership without another admission. Minimum request size
    applies to initial admission, not to individual components of an admitted job.
  */
  class ByteReservation
  {
  public:
    ByteReservation() noexcept = default;
    ~ByteReservation();
    ByteReservation(const ByteReservation&) = delete;
    ByteReservation& operator=(const ByteReservation&) = delete;
    ByteReservation(ByteReservation&& other) noexcept;
    ByteReservation& operator=(ByteReservation&& other) noexcept;
    ByteCount bytes() const noexcept { return bytes_; }
    void reset() noexcept;
    /** @brief Transfer part or all of this reservation without blocking.
      @param[in] nonzero_bytes Positive amount no greater than the remaining reservation.
      @throws OpenMS::Exception::InvalidValue for zero or oversized splits.
    */
    ByteReservation split(ByteCount nonzero_bytes);
    /** @brief Convert pending bytes to a unique committed charge and release unused bytes.
      @param[in] actual_bytes Charge no greater than this reservation. Zero releases all bytes.
      @throws OpenMS::Exception::InvalidValue if the charge exceeds the reservation.
    */
    MemoryLease commit(ByteCount actual_bytes);

  private:
    friend class ByteBudget;
    ByteReservation(std::shared_ptr<Detail::BudgetState> state, ByteCount bytes) noexcept;
    std::shared_ptr<Detail::BudgetState> state_;
    ByteCount bytes_ = 0;
  };

  /**
    @brief One fixed native allowance plus disjoint pending and committed bytes.

    Caller-declared charges are accounting, not interception of allocations or
    an RSS limit. Reserve bounded downstream drain, fixed workspaces and queue
    nodes before producers. Atomically admit each whole peak job before allocating
    its input/output/scratch. Split that reservation without further blocking
    admission while holding the job. Consumers drain with pre-admitted workspace.
    Acyclic forward queues and this acquisition order are caller responsibilities.

    Cancellation permanently rejects new admission and wakes all waiters, but
    never revokes live reservations or charges. Join all callers before destroying
    the budget. Tokens may outlive it. No fairness or admission deadline is promised.
  */
  class ByteBudget
  {
  public:
    /** @brief Configure admission bounds without allocating proportional payload storage.
      @param[in] total_bytes Positive total including native allowance.
      @param[in] native_reserve Runtime/native allowance, strictly less than total.
      @param[in] minimum_request Minimum initial admission, within the usable remainder.
      @throws OpenMS::Exception::InvalidValue for impossible configuration.
    */
    ByteBudget(ByteCount total_bytes, ByteCount native_reserve, ByteCount minimum_request = 1);
    ~ByteBudget();
    ByteBudget(const ByteBudget&) = delete;
    ByteBudget& operator=(const ByteBudget&) = delete;
    /** @brief Attempt immediate admission, returning nullopt if full or cancelled.
      @param[in] bytes Initial request within configured minimum and usable bounds.
      @throws OpenMS::Exception::InvalidValue for invalid requests, even after cancellation.
    */
    std::optional<ByteReservation> tryReserve(ByteCount bytes);
    /** @brief Wait for admission, returning nullopt on stop or permanent cancellation.
      @param[in] bytes Initial request within configured minimum and usable bounds.
      @param[in] stop Optional per-call cancellation, which does not cancel other callers.
      @throws OpenMS::Exception::InvalidValue for invalid requests, without waiting.
    */
    std::optional<ByteReservation> reserve(ByteCount bytes, std::stop_token stop = {});
    void cancel();
    BudgetSnapshot snapshot() const;

  private:
    std::shared_ptr<Detail::BudgetState> state_;
  };
}
