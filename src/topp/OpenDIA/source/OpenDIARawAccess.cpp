// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARawAccess.h"
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/Precursor.h>
#include <algorithm>
#include <chrono>
#include <cmath>

namespace OpenDIA
{
  namespace
  {
    void invalid(const char* text)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, text, "raw scan cache");
    }

    bool lowerClosed(const QueryInterval& q)
    {
      return q.edges == EdgePolicy::CLOSED || q.edges == EdgePolicy::LOWER_CLOSED;
    }

    bool upperClosed(const QueryInterval& q)
    {
      return q.edges == EdgePolicy::CLOSED || q.edges == EdgePolicy::UPPER_CLOSED;
    }

    void validate(const QueryInterval& q)
    {
      if (!std::isfinite(q.lower) || !std::isfinite(q.upper) || q.lower > q.upper) { invalid("Invalid interval"); }
      switch (q.edges)
      {
        case EdgePolicy::OPEN: case EdgePolicy::CLOSED:
        case EdgePolicy::LOWER_CLOSED: case EdgePolicy::UPPER_CLOSED: break;
        default: invalid("Invalid edge policy");
      }
    }

    bool contains(const QueryInterval& q, double value)
    {
      return (value > q.lower || (lowerClosed(q) && value == q.lower)) &&
        (value < q.upper || (upperClosed(q) && value == q.upper));
    }

    bool overlaps(const QueryInterval& q, double lower, double upper)
    {
      if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper) { return false; }
      if (q.lower == q.upper && !(lowerClosed(q) && upperClosed(q))) { return false; }
      return (upper > q.lower || (upper == q.lower && lowerClosed(q))) &&
        (lower < q.upper || (lower == q.upper && upperClosed(q)));
    }

    using Clock = std::chrono::steady_clock;
    ByteCount elapsed(Clock::time_point start)
    {
      return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
    }

    ByteCount numericBytes(const OpenSwath::Spectrum& spectrum)
    {
      ByteCount bytes = 0;
      for (const auto& array : spectrum.getDataArrays())
      {
        bytes = checkedByteAdd(bytes, checkedByteMultiply(array->data.size(), sizeof(double)));
      }
      return bytes;
    }

    ByteCount charge(const OpenSwath::Spectrum& spectrum)
    {
      // Count capacity, not size. Include the string terminator.
      const auto& arrays = spectrum.getDataArrays();
      ByteCount bytes = checkedByteAdd(sizeof(RawScan), sizeof(OpenSwath::Spectrum));
      bytes = checkedByteAdd(bytes, checkedByteMultiply(arrays.capacity(), sizeof(OpenSwath::BinaryDataArrayPtr)));
      for (const auto& a : arrays)
      {
        if (!a) { invalid("Null numeric column"); }
        bytes = checkedByteAdd(bytes, sizeof(OpenSwath::BinaryDataArray));
        bytes = checkedByteAdd(bytes, checkedByteMultiply(a->data.capacity(), sizeof(double)));
        bytes = checkedByteAdd(bytes, checkedByteAdd(a->description.capacity(), 1));
      }
      return bytes;
    }

    ByteCount scanBytes(const OpenMS::MSSpectrum& spectrum)
    {
      const auto count = checkedByteAdd(2, checkedByteAdd(spectrum.getFloatDataArrays().size(), spectrum.getIntegerDataArrays().size()));
      ByteCount bytes = checkedByteAdd(sizeof(RawScan), sizeof(OpenSwath::Spectrum));
      bytes = checkedByteAdd(bytes, checkedByteMultiply(count, sizeof(OpenSwath::BinaryDataArrayPtr)));
      const auto add = [&](OpenMS::Size size, const std::string& name)
      {
        bytes = checkedByteAdd(bytes, sizeof(OpenSwath::BinaryDataArray));
        bytes = checkedByteAdd(bytes, checkedByteMultiply(size, sizeof(double)));
        bytes = checkedByteAdd(bytes, checkedByteAdd(std::max(std::string{}.capacity(), name.size()), 1));
      };
      add(spectrum.size(), "");
      add(spectrum.size(), "");
      for (const auto& array : spectrum.getFloatDataArrays()) { add(array.size(), array.getName()); }
      for (const auto& array : spectrum.getIntegerDataArrays()) { add(array.size(), array.getName()); }
      return bytes;
    }

    template <typename Arrays>
    void appendArrays(const Arrays& arrays, OpenSwath::Spectrum& output)
    {
      for (const auto& array : arrays)
      {
        auto column = std::make_shared<OpenSwath::BinaryDataArray>();
        column->description = array.getName();
        column->data.assign(array.begin(), array.end());
        output.getDataArrays().push_back(std::move(column));
      }
    }

    // Each array alias keeps all columns and their charge alive.
    struct MutableColumns
    {
      MemoryLease lease;
      OpenSwath::SpectrumPtr spectrum;
    };

    struct MutableSpectrum
    {
      MemoryLease lease;
      OpenSwath::Spectrum spectrum;
    };
  }

  const OpenSwath::BinaryDataArray& RawScan::array(OpenMS::Size index) const
  {
    if (index >= arrayCount())
    {
      throw OpenMS::Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, arrayCount());
    }
    return *columns_->getDataArrays()[index];
  }

  RawCache::RawCache(ByteCount total_bytes, ByteCount native_reserve, ByteCount max_scan_bytes) :
    budget_(total_bytes, native_reserve), max_scan_bytes_(max_scan_bytes)
  {
    if (max_scan_bytes == 0 || max_scan_bytes > total_bytes - native_reserve) { invalid("Invalid maximum scan size"); }
  }

  std::optional<ByteReservation> RawCache::reserve_(ByteCount bytes, RawPhase phase)
  {
    // Evict on pressure. Pinned scans keep leases, so admission can still fail.
    auto reservation = budget_.tryReserve(bytes);
    while (!reservation && !entries_.empty() && !budget_.snapshot().cancelled)
    {
      entries_.pop_back();
      ++work_.evictions;
      ++phases_[static_cast<size_t>(phase)].evictions;
      reservation = budget_.tryReserve(bytes);
    }
    return reservation;
  }

  RawCacheSnapshot RawCache::snapshot() const
  {
    std::lock_guard lock(mutex_);
    return {budget_.snapshot(), entries_.size(), work_.hits, work_.misses, work_.evictions, work_, phases_};
  }

  void RawCache::setCalibration(bool active)
  {
    std::lock_guard lock(mutex_);
    calibrating_ = active;
  }

  void RawCache::clear()
  {
    std::lock_guard lock(mutex_);
    entries_.clear();
  }

  void RawCache::cancel()
  {
    std::lock_guard lock(mutex_);
    budget_.cancel();
    entries_.clear();
  }

  RawAccess::RawAccess(std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source, std::shared_ptr<RawCache> cache) :
    source_(std::move(source)), cache_(std::move(cache))
  {
    if (!source_ || !cache_) { invalid("Null source or cache"); }
    std::lock_guard lock(cache_->mutex_);
    // One metadata handle is enough. The cache lock also protects this handle.
    metadata_source_ = source_->getExperiment().clone();
    scan_bounds_ = std::make_shared<std::vector<ByteCount>>(source_->getNrSpectra(), 0);
    region_bounds_ = std::make_shared<std::vector<ByteCount>>(source_->getNrSpectra(), 0);
  }

  RawPhase RawAccess::phase_(int id) const
  {
    if (cache_->calibrating_) return RawPhase::CALIBRATION;
    using Scope = OpenSwath::SpectrumAccessScope;
    if (Scope::current() == Scope::Purpose::SCORING) return RawPhase::SCORING;
    if (Scope::current() == Scope::Purpose::EXTRACTION)
    {
      switch (source_->getSpectrumMetaById(id).ms_level)
      {
        case 1: return RawPhase::MS1_EXTRACTION;
        case 2: return RawPhase::MS2_EXTRACTION;
      }
    }
    return RawPhase::UNCLASSIFIED;
  }

  std::shared_ptr<const RawScan> RawAccess::acquire(int id)
  {
    const auto started = Clock::now();
    std::lock_guard lock(cache_->mutex_);
    const auto index = source_->getExperimentIndex(id);
    const auto phase_id = phase_(id);
    auto& work = cache_->work_;
    auto& phase = cache_->phases_[static_cast<size_t>(phase_id)];
    const auto finish = [&]
    {
      const auto ns = elapsed(started);
      work.read_ns += ns;
      phase.read_ns += ns;
    };
    if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
    auto& entries = cache_->entries_;
    // ponytail: one lock and linear LRU lookup, add an index only if profiling requires it.
    auto found = std::find_if(entries.begin(), entries.end(), [&](const auto& e) { return e.scan && e.source == source_ && e.id == id; });
    if (found != entries.end())
    {
      ++work.hits;
      ++phase.hits;
      entries.splice(entries.begin(), entries, found);
      finish();
      return entries.front().scan;
    }
    ++work.misses;
    ++phase.misses;
    // First reads reserve the maximum. Immutable revisits reuse their proven
    // capacity bound instead of evicting hot scans for unused maximum headroom.
    auto& known_bound = scan_bounds_->at(static_cast<size_t>(id));
    const auto scan_bound = known_bound ? known_bound : cache_->max_scan_bytes_;
    auto reservation = cache_->reserve_(scan_bound, phase_id);
    if (!reservation) { invalid("Scan admission failed: pinned views or cancellation"); }
    // Native decode is outside private capacity enforcement. Check before conversion.
    const auto decode_started = Clock::now();
    const auto native = metadata_source_->getSpectrumPeaks(index);
    const auto decode_ns = elapsed(decode_started);
    ++work.decodes;
    ++phase.decodes;
    work.decode_ns += decode_ns;
    phase.decode_ns += decode_ns;
    const auto native_bound = scanBytes(native);
    if (native_bound > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (native_bound > scan_bound) { invalid("Scan exceeds previously admitted scan bytes"); }
    auto scan = std::shared_ptr<RawScan>(new RawScan);
    scan->metadata_ = metadata_source_->getSpectrumMetadata(index);
    scan->columns_ = OpenMS::OpenSwathDataAccessHelper::convertToSpectrumPtr(native);
    scan->columns_->getDataArrays().reserve(2 + native.getFloatDataArrays().size() + native.getIntegerDataArrays().size());
    appendArrays(native.getFloatDataArrays(), *scan->columns_);
    appendArrays(native.getIntegerDataArrays(), *scan->columns_);
    const auto bytes = charge(*scan->columns_);
    if (bytes > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (bytes > scan_bound) { invalid("Scan exceeds previously admitted scan bytes"); }
    scan->lease_ = reservation->commit(bytes);
    // Publish only a complete scan whose retained capacities fit the bound.
    entries.push_front({source_, id, scan});
    known_bound = std::max(native_bound, bytes);
    const auto numeric_bytes = numericBytes(*scan->columns_);
    work.decoded_numeric_bytes = checkedByteAdd(work.decoded_numeric_bytes, numeric_bytes);
    phase.decoded_numeric_bytes += numeric_bytes;
    finish();
    return scan;
  }

  RawAccess::SpectrumView RawAccess::getSpectrumView(int id)
  {
    const auto scan = acquire(id);
    const auto mobility = scan->columns_->getDriftTimeArray();
    return {std::shared_ptr<const OpenSwath::BinaryDataArray>(scan, scan->columns_->getMZArray().get()),
      std::shared_ptr<const OpenSwath::BinaryDataArray>(scan, scan->columns_->getIntensityArray().get()),
      mobility ? std::shared_ptr<const OpenSwath::BinaryDataArray>(scan, mobility.get()) : nullptr};
  }

  RawRegion RawAccess::acquireRegion(int id, RecordContext context, std::optional<MZCalibration> calibration)
  {
    if (context.schema_version != 0
        || (calibration && (! std::isfinite(calibration->a) || ! std::isfinite(calibration->b) || ! std::isfinite(calibration->c))))
    {
      invalid("Invalid region context or calibration");
    }
    const auto started = Clock::now();
    std::lock_guard lock(cache_->mutex_);
    const auto index = source_->getExperimentIndex(id);
    const auto phase_id = phase_(id);
    auto& work = cache_->work_;
    auto& phase = cache_->phases_[static_cast<size_t>(phase_id)];
    if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
    auto& entries = cache_->entries_;
    const auto found = std::find_if(entries.begin(), entries.end(), [&](const auto& e) {
      return e.region && e.source == source_ && e.id == id && e.context == context && e.calibration == calibration;
    });
    if (found != entries.end())
    {
      ++work.hits;
      ++phase.hits;
      entries.splice(entries.begin(), entries, found);
      const auto ns = elapsed(started);
      work.read_ns += ns;
      phase.read_ns += ns;
      return *entries.front().region;
    }
    ++work.misses;
    ++phase.misses;
    auto& known_bound = region_bounds_->at(static_cast<size_t>(id));
    const auto bound = known_bound ? known_bound : cache_->max_scan_bytes_;
    auto reservation = cache_->reserve_(bound, phase_id);
    if (! reservation) { invalid("Region admission failed: pinned views or cancellation"); }
    const auto decode_started = Clock::now();
    const auto native = metadata_source_->getSpectrum(index);
    const auto decode_ns = elapsed(decode_started);
    ++work.decodes;
    ++phase.decodes;
    work.decode_ns += decode_ns;
    phase.decode_ns += decode_ns;
    const auto bytes = RawRegion::requiredBytes({1, native.size()}, native.getNativeID().size());
    if (bytes > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (bytes > bound) { invalid("Region exceeds previously admitted scan bytes"); }
    // The shared cache lock excludes other admissions during this handoff.
    // Native decode and sort workspace remain outside the retained-column cap.
    reservation.reset();
    auto region = RawRegion::fromSpectrum(native, index, context, cache_->budget_, calibration);
    entries.push_front({source_, id, {}, region, context, calibration});
    known_bound = bytes;
    const auto numeric_bytes = checkedByteMultiply(native.size(), sizeof(double) + sizeof(float) * (native.containsIMData() ? 2 : 1));
    work.decoded_numeric_bytes = checkedByteAdd(work.decoded_numeric_bytes, numeric_bytes);
    phase.decoded_numeric_bytes += numeric_bytes;
    const auto ns = elapsed(started);
    work.read_ns += ns;
    phase.read_ns += ns;
    return region;
  }

  OpenSwath::SpectrumPtr RawAccess::getSpectrumById(int id)
  {
    const auto scan = acquire(id);
    std::lock_guard lock(cache_->mutex_);
    const auto started = Clock::now();
    const auto phase_id = phase_(id);
    auto& work = cache_->work_;
    auto& phase = cache_->phases_[static_cast<size_t>(phase_id)];
    // Bound the detached copy from the pinned columns, not the maximum scan size.
    const auto copy_bound = checkedByteAdd(checkedByteAdd(charge(*scan->columns_), sizeof(MutableColumns)),
      checkedByteAdd(sizeof(MutableSpectrum), checkedByteMultiply(scan->arrayCount(), sizeof(OpenSwath::BinaryDataArrayPtr))));
    if (copy_bound > cache_->max_scan_bytes_) { invalid("Mutable scan exceeds maximum scan bytes"); }
    auto reservation = cache_->reserve_(copy_bound, phase_id);
    if (!reservation) { invalid("Mutable scan admission failed"); }
    auto owner = std::make_shared<MutableColumns>();
    owner->spectrum = std::make_shared<OpenSwath::Spectrum>();
    auto& out = owner->spectrum->getDataArrays();
    out.clear();
    out.reserve(scan->arrayCount());
    for (OpenMS::Size i = 0; i < scan->arrayCount(); ++i)
    {
      out.push_back(std::make_shared<OpenSwath::BinaryDataArray>(scan->array(i)));
    }
    // Return detached columns. An array-only owner must retain the lease.
    auto result = std::make_shared<MutableSpectrum>();
    auto& aliases = result->spectrum.getDataArrays();
    aliases.clear();
    aliases.reserve(out.size());
    for (const auto& array : out) { aliases.emplace_back(owner, array.get()); }
    const auto column_bytes = checkedByteAdd(charge(*owner->spectrum), sizeof(MutableColumns));
    const auto result_bytes = checkedByteAdd(sizeof(MutableSpectrum),
      checkedByteMultiply(aliases.capacity(), sizeof(OpenSwath::BinaryDataArrayPtr)));
    if (checkedByteAdd(column_bytes, result_bytes) > copy_bound)
    {
      invalid("Mutable scan exceeds reserved copy bytes");
    }
    // Split ownership. A caller can clear every array but retain the spectrum.
    auto columns_reservation = reservation->split(column_bytes);
    owner->lease = columns_reservation.commit(column_bytes);
    result->lease = reservation->commit(result_bytes);
    const auto numeric_bytes = numericBytes(*owner->spectrum);
    const auto ns = elapsed(started);
    ++work.copies;
    ++phase.copies;
    work.copied_numeric_bytes = checkedByteAdd(work.copied_numeric_bytes, numeric_bytes);
    phase.copied_numeric_bytes += numeric_bytes;
    work.copy_ns += ns;
    phase.copy_ns += ns;
    return OpenSwath::SpectrumPtr(result, &result->spectrum);
  }

  std::vector<std::size_t> RawAccess::select(const QueryInterval& rt,
    std::optional<QueryInterval> isolation, std::optional<double> faims)
  {
    validate(rt);
    if (isolation) { validate(*isolation); }
    if (faims && !std::isfinite(*faims)) { invalid("Invalid FAIMS voltage"); }
    std::lock_guard lock(cache_->mutex_);
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < source_->getNrSpectra(); ++i)
    {
      if (!contains(rt, source_->getSpectrumMetaById(static_cast<int>(i)).RT)) { continue; }
      const auto metadata = metadata_source_->getSpectrumMetadata(source_->getExperimentIndex(static_cast<int>(i)));
      if (faims && (metadata.getDriftTimeUnit() != OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE || metadata.getDriftTime() != *faims)) { continue; }
      // Test full isolation geometry, not the locations of observed peaks.
      if (isolation)
      {
        const auto& precursors = metadata.getPrecursors();
        if (!std::any_of(precursors.begin(), precursors.end(), [&](const auto& p)
          { return overlaps(*isolation, p.getMZ() - p.getIsolationWindowLowerOffset(), p.getMZ() + p.getIsolationWindowUpperOffset()); })) { continue; }
      }
      result.push_back(i);
    }
    return result;
  }

  std::shared_ptr<OpenSwath::ISpectrumAccess> RawAccess::lightClone() const
  {
    // Share storage and both native handles. Do not clone the cache.
    return std::make_shared<RawAccess>(*this);
  }

  OpenSwath::SpectrumMeta RawAccess::getSpectrumMetaById(int id) const
  {
    std::lock_guard lock(cache_->mutex_);
    return source_->getSpectrumMetaById(id);
  }

  std::vector<std::size_t> RawAccess::getSpectraByRT(double rt, double delta_rt) const
  {
    if (!std::isfinite(rt) || !std::isfinite(delta_rt) || delta_rt < 0) { invalid("Invalid RT lookup"); }
    std::lock_guard lock(cache_->mutex_);
    return source_->getSpectraByRT(rt, delta_rt);
  }

  size_t RawAccess::getNrSpectra() const
  {
    std::lock_guard lock(cache_->mutex_);
    return source_->getNrSpectra();
  }

  OpenSwath::ChromatogramPtr RawAccess::getChromatogramById(int id)
  {
    std::lock_guard lock(cache_->mutex_);
    if (id < 0 || static_cast<size_t>(id) >= source_->getNrChromatograms()) { invalid("Invalid chromatogram index"); }
    return source_->getChromatogramById(id);
  }

  size_t RawAccess::getNrChromatograms() const
  {
    std::lock_guard lock(cache_->mutex_);
    return source_->getNrChromatograms();
  }

  std::string RawAccess::getChromatogramNativeID(int id) const
  {
    std::lock_guard lock(cache_->mutex_);
    if (id < 0 || static_cast<size_t>(id) >= source_->getNrChromatograms()) { invalid("Invalid chromatogram index"); }
    return source_->getChromatogramNativeID(id);
  }
}
