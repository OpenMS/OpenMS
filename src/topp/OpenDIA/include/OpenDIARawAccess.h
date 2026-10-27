// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include "OpenDIARecords.h"
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOnDisc.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <array>
#include <list>
#include <mutex>

namespace OpenDIA
{
  enum class RawPhase { UNCLASSIFIED, CALIBRATION, MS1_EXTRACTION, MS2_EXTRACTION, SCORING, COUNT };

  /// Successful read/copy timings are summed per call and may overlap across threads.
  /// Numeric bytes count materialized double values, not compressed input or native workspace.
  struct RawWorkCounters
  {
    ByteCount hits = 0, misses = 0, evictions = 0, decodes = 0, copies = 0;
    ByteCount decoded_numeric_bytes = 0, copied_numeric_bytes = 0;
    ByteCount read_ns = 0, decode_ns = 0, copy_ns = 0;
  };

  struct RawCacheSnapshot
  {
    BudgetSnapshot budget;
    OpenMS::Size cache_entries = 0, hits = 0, misses = 0, evictions = 0;
    RawWorkCounters work;
    std::array<RawWorkCounters, static_cast<size_t>(RawPhase::COUNT)> phases{};
  };

  /// Read-only scan. The lease stays valid after cache removal.
  class RawScan final
  {
  public:
    /// Native metadata has no peak or per-point arrays.
    const OpenMS::MSSpectrum& metadata() const { return metadata_; }
    /// Count all numeric columns, including m/z and intensity.
    OpenMS::Size arrayCount() const { return columns_->getDataArrays().size(); }
    /** @brief Read one numeric column.
      @param[in] index Column index, less than arrayCount().
    */
    const OpenSwath::BinaryDataArray& array(OpenMS::Size index) const;
  private:
    friend class RawAccess;
    RawScan() = default;
    MemoryLease lease_;
    OpenMS::MSSpectrum metadata_;
    OpenSwath::SpectrumPtr columns_;
  };

  /** @brief Bounded numeric scan storage.
    Native decode, metadata, index and control allocations use the declared native
    allowance. This allowance is not an intercepted allocation or RSS bound.
    Numeric charges use vector and string capacities. Admission never waits.
    Charges do not track later growth of mutable legacy copies.
  */
  class RawCache final
  {
  public:
    /** @brief Set the admission bounds.
      @param[in] total_bytes Total bytes, including native allowance.
      @param[in] native_reserve Fixed native allowance, below total_bytes.
      @param[in] max_scan_bytes Positive scan admission bound, within usable bytes.
    */
    RawCache(ByteCount total_bytes, ByteCount native_reserve, ByteCount max_scan_bytes);
    /// Read counters and budget state under one lock.
    RawCacheSnapshot snapshot() const;
    /// Change only between joined workflow stages, never while readers are active.
    void setCalibration(bool active);
    /// Remove cached scans. Live leases stay valid. This does not undo cancel().
    void clear();
    /// Reject new admission and remove cached scans. Live leases stay valid.
    void cancel();
  private:
    friend class RawAccess;
    struct Entry
    {
      std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source;
      int id;
      std::shared_ptr<const RawScan> scan;
    };
    std::optional<ByteReservation> reserve_(ByteCount bytes, RawPhase phase);
    mutable std::mutex mutex_;
    ByteBudget budget_;
    ByteCount max_scan_bytes_;
    std::list<Entry> entries_;
    RawWorkCounters work_;
    std::array<RawWorkCounters, static_cast<size_t>(RawPhase::COUNT)> phases_{};
    bool calibrating_ = false;
  };

  /// Private scan adapter. Clones share the cache and native handle lock.
  class RawAccess final : public OpenSwath::ISpectrumAccess
  {
  public:
    /** @brief Share one native source and cache.
      All access to this source must use this cache lock. Do not use the supplied
      source directly while an adapter call is active.
      Source scan contents and indexing must remain unchanged for the adapter lifetime.
      @param[in] source Non-null native source.
      @param[in] cache Non-null cache shared by all clones.
    */
    RawAccess(std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source, std::shared_ptr<RawCache> cache);
    /** @brief Lease one complete scan. Fail if admission is not possible.
      @param[in] id Visible scan index.
    */
    std::shared_ptr<const RawScan> acquire(int id);
    /** @brief Select scans by RT and full acquisition isolation overlap.
      Source isolation intervals are closed. Query edges follow their policy.
      @param[in] rt RT bounds in seconds.
      @param[in] isolation Optional isolation bounds in Th.
      @param[in] faims Optional exact FAIMS voltage in volts.
    */
    std::vector<std::size_t> select(const QueryInterval& rt,
      std::optional<QueryInterval> isolation = std::nullopt, std::optional<double> faims = std::nullopt);
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;
    /// Alias immutable cached arrays through the scan lease, without a detached numeric copy.
    SpectrumView getSpectrumView(int id) override;
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;
    std::vector<std::size_t> getSpectraByRT(double rt, double delta_rt) const override;
    size_t getNrSpectra() const override;
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;
    size_t getNrChromatograms() const override;
    std::string getChromatogramNativeID(int id) const override;
  private:
    RawPhase phase_(int id) const; // cache lock held
    std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source_;
    std::shared_ptr<RawCache> cache_;
    std::shared_ptr<OpenMS::IMSExperimentRandomAccess> metadata_source_;
    // Immutable-source sizing metadata, shared by clones under the cache lock.
    // One ByteCount per visible scan uses the native metadata/control allowance.
    std::shared_ptr<std::vector<ByteCount>> scan_bounds_;
  };
}
