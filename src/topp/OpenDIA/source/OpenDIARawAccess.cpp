// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
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
#include <limits>
#include <string>
#include <unordered_map>

namespace OpenDIA
{
  namespace
  {
    void invalid(const std::string& text)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, text, "raw scan cache");
    }

    std::string admissionFailure(const char* resource, ByteCount requested, const BudgetSnapshot& budget)
    {
      const ByteCount used = checkedByteAdd(budget.reserved_bytes, budget.committed_bytes);
      const ByteCount available = budget.usable_bytes > used ? budget.usable_bytes - used : 0;
      return std::string(resource) + " admission failed: requested " + std::to_string(requested) +
        " bytes, " + std::to_string(available) + " bytes available after cache eviction. " +
        (budget.cancelled ? "The cache was cancelled." :
          "Live views may pin memory. Reduce -threads or increase "
          "TargetedDataExtraction:raw_memory_mb.");
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
    struct BatchAdmissionFailure
    {
    };
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

  struct RawCache::BatchFrame
  {
    struct Member
    {
      std::weak_ptr<RawAccess> access;
      int id = 0;
      OpenMS::Size index = 0;
    };
    std::shared_ptr<OpenMS::IMSExperimentRandomAccess> source;
    std::vector<OpenMS::Size> indices;
    std::unordered_map<OpenMS::Size, std::vector<Member>> members_by_index;
  };

  const OpenSwath::BinaryDataArray& RawScan::array(OpenMS::Size index) const
  {
    if (index >= arrayCount())
    {
      throw OpenMS::Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, arrayCount());
    }
    return *columns_->getDataArrays()[index];
  }

  bool RawScan::validateForExtraction() const
  {
    if (!std::isfinite(metadata_.getRT()) || metadata_.getMSLevel() == 0 || !native_points_valid_)
    {
      invalid("Invalid acquisition metadata or native points");
    }
    const auto scalar_unit = metadata_.getDriftTimeUnit();
    if (scalar_unit != OpenMS::DriftTimeUnit::NONE)
    {
      if (!std::isfinite(metadata_.getDriftTime()) ||
          (scalar_unit != OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE &&
           scalar_unit != OpenMS::DriftTimeUnit::MILLISECOND && scalar_unit != OpenMS::DriftTimeUnit::VSSC))
      {
        invalid("Unsupported or invalid scalar mobility");
      }
    }
    if (point_mobility_unit_ && *point_mobility_unit_ != OpenMS::DriftTimeUnit::MILLISECOND &&
        *point_mobility_unit_ != OpenMS::DriftTimeUnit::VSSC)
    {
      invalid("Unsupported point mobility unit");
    }
    return masses_sorted_;
  }

  RawCache::RawCache(ByteCount total_bytes, ByteCount native_reserve, ByteCount max_scan_bytes) :
    budget_(total_bytes, native_reserve), max_scan_bytes_(max_scan_bytes)
  {
    if (max_scan_bytes == 0 || max_scan_bytes > total_bytes - native_reserve) { invalid("Invalid maximum scan size"); }
  }

  RawCache::~RawCache()
  {
    clear();
  }

  bool RawCache::enableNativeFrames(std::shared_ptr<OpenMS::IMSExperimentRandomAccess> source)
  {
    if (! source) { invalid("Null native frame source"); }
    std::lock_guard lock(mutex_);
    if (native_source_) { return false; }
    if (! entries_.empty()) { invalid("Configure native frames before readers start"); }
    const auto usable = budget_.snapshot().usable_bytes;
    // ponytail: one source and at most one quarter of usable bytes. Revisit only with measured multi-source demand.
    const auto bytes = std::min(max_scan_bytes_, usable / 4);
    if (! bytes || usable - bytes < max_scan_bytes_) { return false; }
    auto reservation = budget_.tryReserve(bytes);
    if (! reservation || ! source->setNativeFrameCacheSize(bytes)) { return false; }
    native_frame_lease_ = reservation->commit(bytes);
    native_source_ = std::move(source);
    return true;
  }

  void RawCache::registerBatchGroup(const std::vector<std::shared_ptr<RawAccess>>& accessors)
  {
    std::lock_guard lock(mutex_);
    if (!entries_.empty()) { invalid("Configure batch groups before readers start"); }
    if (budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
    if (accessors.empty()) { return; }

    std::unordered_map<OpenMS::Size, std::vector<BatchFrame::Member>> members_by_acquisition;
    std::shared_ptr<OpenMS::IMSExperimentRandomAccess> source;
    for (const auto& access : accessors)
    {
      if (!access || access->cache_.get() != this) { invalid("Batch lanes must share this cache"); }
      if (!source) { source = access->metadata_source_; }
      const auto count = access->source_->getNrSpectra();
      if (count > static_cast<OpenMS::Size>(std::numeric_limits<int>::max())) { invalid("Too many batch lane spectra"); }
      access->batch_frames_by_id_.assign(count, {});
      for (OpenMS::Size i = 0; i < count; ++i)
      {
        const auto native_index = access->source_->getExperimentIndex(static_cast<int>(i));
        const auto acquisition = access->metadata_source_->getAcquisitionBatchId(native_index);
        if (acquisition)
        {
          members_by_acquisition[*acquisition].push_back({access, static_cast<int>(i), native_index});
        }
      }
    }

    if (members_by_acquisition.empty()) { return; }
    std::shared_ptr<OpenMS::IMSExperimentRandomAccess> batch_source(source->clone());
    batch_sources_.push_back(batch_source);
    for (auto& [acquisition, members] : members_by_acquisition)
    {
      static_cast<void>(acquisition);
      if (members.size() < 2) { continue; }
      auto frame = std::make_shared<BatchFrame>();
      frame->source = batch_source;
      for (const auto& member : members)
      {
        frame->members_by_index[member.index].push_back(member);
      }
      frame->indices.reserve(frame->members_by_index.size());
      for (const auto& [index, aliases] : frame->members_by_index)
      {
        static_cast<void>(aliases);
        frame->indices.push_back(index);
      }
      std::sort(frame->indices.begin(), frame->indices.end());
      for (const auto& member : members)
      {
        if (const auto access = member.access.lock())
        {
          access->batch_frames_by_id_[static_cast<OpenMS::Size>(member.id)] = frame;
        }
      }
      batch_frames_.push_back(std::move(frame));
    }
    if (batch_frames_.empty()) { batch_sources_.pop_back(); }
  }

  void RawCache::clearNativeFrames_()
  {
    if (! native_source_) { return; }
    native_frame_stats_ = native_source_->getNativeFrameCacheStats();
    native_source_->setNativeFrameCacheSize(0);
    native_frame_stats_.retained_bytes = 0;
    native_source_.reset();
    native_frame_lease_.reset();
  }

  std::optional<ByteReservation> RawCache::reserve_(ByteCount bytes, RawPhase phase)
  {
    // Optional decoded frames must not displace complete scans. Once pressure
    // exhausts headroom, release their reservation before evicting scan entries.
    auto reservation = budget_.tryReserve(bytes);
    if (!reservation && native_source_ && !budget_.snapshot().cancelled)
    {
      clearNativeFrames_();
      reservation = budget_.tryReserve(bytes);
    }
    // Pinned scans keep leases, so admission can still fail after reclamation.
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
    OpenMS::NativeFrameCacheStats frame_stats = native_frame_stats_;
    if (!batch_sources_.empty())
    {
      frame_stats = {};
      for (const auto& source : batch_sources_)
      {
        const auto stats = source->getNativeFrameCacheStats();
        frame_stats.limit_bytes = checkedByteAdd(frame_stats.limit_bytes, stats.limit_bytes);
        frame_stats.retained_bytes = checkedByteAdd(frame_stats.retained_bytes, stats.retained_bytes);
        frame_stats.peak_retained_bytes = checkedByteAdd(frame_stats.peak_retained_bytes, stats.peak_retained_bytes);
        frame_stats.hits = checkedByteAdd(frame_stats.hits, stats.hits);
        frame_stats.misses = checkedByteAdd(frame_stats.misses, stats.misses);
        frame_stats.evictions = checkedByteAdd(frame_stats.evictions, stats.evictions);
        frame_stats.decoded_frames = checkedByteAdd(frame_stats.decoded_frames, stats.decoded_frames);
        frame_stats.sorted_frames = checkedByteAdd(frame_stats.sorted_frames, stats.sorted_frames);
      }
    }
    else if (native_source_)
    {
      frame_stats = native_source_->getNativeFrameCacheStats();
    }
    return {budget_.snapshot(), entries_.size(), work_.hits, work_.misses,
            work_.evictions,    work_,           phases_,    frame_stats};
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
    clearNativeFrames_();
  }

  void RawCache::cancel()
  {
    std::lock_guard lock(mutex_);
    budget_.cancel();
    entries_.clear();
    clearNativeFrames_();
  }

  RawAccess::RawAccess(std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source, std::shared_ptr<RawCache> cache) :
    source_(std::move(source)), cache_(std::move(cache)), source_mutex_(std::make_shared<std::mutex>())
  {
    if (!source_ || !cache_) { invalid("Null source or cache"); }
    std::lock_guard lock(*source_mutex_);
    // One metadata handle is enough. The source lock protects both native handles.
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

  std::shared_ptr<const RawScan> RawAccess::acquire(int id, std::optional<MZCalibration> calibration)
  {
    if (calibration && (!std::isfinite(calibration->a) || !std::isfinite(calibration->b) || !std::isfinite(calibration->c)))
    {
      invalid("Invalid scan calibration");
    }
    const auto started = Clock::now();
    std::lock_guard source_lock(*source_mutex_);
    std::unique_lock cache_lock(cache_->mutex_);
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
    auto found = std::find_if(entries.begin(), entries.end(), [&](const auto& e) {
      return e.scan && e.source == source_ && e.id == id && e.calibration == calibration;
    });
    if (found != entries.end())
    {
      ++work.hits;
      ++phase.hits;
      entries.splice(entries.begin(), entries, found);
      finish();
      return entries.front().scan;
    }
    std::shared_ptr<RawCache::BatchFrame> frame;
    if (static_cast<OpenMS::Size>(id) < batch_frames_by_id_.size())
    {
      frame = batch_frames_by_id_[static_cast<OpenMS::Size>(id)].lock();
    }
    // Extraction streams sibling lanes and benefits from physical-frame batches.
    // Point scoring/calibration reads must not populate every sibling on a miss.
    if (frame && phase_id != RawPhase::SCORING && phase_id != RawPhase::CALIBRATION)
    {
      struct StagedScan
      {
        OpenMS::Size index;
        ByteCount bound;
        std::shared_ptr<const RawScan> scan;
      };
      std::vector<StagedScan> staged;
      auto admission = cache_->reserve_(cache_->max_scan_bytes_, phase_id);
      const auto decode_started = Clock::now();
      cache_lock.unlock();
      try
      {
        if (!admission) { throw BatchAdmissionFailure {}; }
        auto remaining = frame->indices.size();
        frame->source->forEachSpectrum(frame->indices, [&](OpenMS::Size native_index, OpenMS::MSSpectrum&& native)
        {
          if (!frame->members_by_index.contains(native_index)) { invalid("Unregistered batch index"); }
          if (!admission) { throw BatchAdmissionFailure {}; }
          const auto native_bound = scanBytes(native);
          auto scan = materialize_(native, native_index, calibration, std::move(*admission));
          admission.reset();
          const auto bytes = charge(*scan->columns_);
          staged.push_back({native_index, std::max(native_bound, bytes), std::move(scan)});
          if (--remaining)
          {
            std::lock_guard lock(cache_->mutex_);
            admission = cache_->reserve_(cache_->max_scan_bytes_, phase_id);
            if (!admission) { throw BatchAdmissionFailure {}; }
          }
        });
      }
      catch (const BatchAdmissionFailure&)
      {
        // Complete batches only. Retained staged leases are bounded by admission,
        // and released before the ordinary one-scan fallback attempts admission.
        staged.clear();
      }
      cache_lock.lock();
      if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
      if (!staged.empty())
      {
        std::shared_ptr<const RawScan> current;
        for (const auto& decoded : staged)
        {
          const auto numeric_bytes = numericBytes(*decoded.scan->columns_);
          ++work.decodes;
          ++phase.decodes;
          work.decoded_numeric_bytes = checkedByteAdd(work.decoded_numeric_bytes, numeric_bytes);
          phase.decoded_numeric_bytes = checkedByteAdd(phase.decoded_numeric_bytes, numeric_bytes);
          bool published_current = false;
          for (const auto& member : frame->members_by_index.at(decoded.index))
          {
            const auto access = member.access.lock();
            if (!access) { continue; }
            ++work.misses;
            ++cache_->phases_[static_cast<size_t>(access->phase_(member.id))].misses;
            access->scan_bounds_->at(static_cast<size_t>(member.id)) = decoded.bound;
            entries.push_front({access->source_, member.id, decoded.scan, {}, {}, calibration});
            if (access->source_ == source_ && member.id == id)
            {
              current = decoded.scan;
              published_current = true;
            }
          }
          // A surviving light clone shares source identity and bounds even after
          // the registered adapter itself is destroyed. Retain its requested row.
          if (decoded.index == index && !published_current)
          {
            ++work.misses;
            ++phase.misses;
            scan_bounds_->at(static_cast<size_t>(id)) = decoded.bound;
            entries.push_front({source_, id, decoded.scan, {}, {}, calibration});
            current = decoded.scan;
          }
        }
        const auto ns = elapsed(decode_started);
        work.decode_ns += ns;
        phase.decode_ns += ns;
        finish();
        if (!current) { invalid("Batch omitted requested scan"); }
        return current;
      }
    }
    ++work.misses;
    ++phase.misses;
    // First reads reserve the maximum. Immutable revisits reuse their proven
    // capacity bound instead of evicting hot scans for unused maximum headroom.
    auto& known_bound = scan_bounds_->at(static_cast<size_t>(id));
    const auto scan_bound = known_bound ? known_bound : cache_->max_scan_bytes_;
    auto reservation = cache_->reserve_(scan_bound, phase_id);
    if (!reservation) { invalid(admissionFailure("Scan", scan_bound, cache_->budget_.snapshot())); }
    // Native decode is outside private capacity enforcement. Check before conversion.
    const auto decode_started = Clock::now();
    cache_lock.unlock();
    const auto native = metadata_source_->getSpectrumPeaks(index);
    const auto decode_ns = elapsed(decode_started);
    const auto native_bound = scanBytes(native);
    if (native_bound > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (native_bound > scan_bound) { invalid("Scan exceeds previously admitted scan bytes"); }
    auto scan = materialize_(native, index, calibration, std::move(*reservation));
    const auto bytes = charge(*scan->columns_);
    cache_lock.lock();
    if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
    ++work.decodes;
    ++phase.decodes;
    work.decode_ns += decode_ns;
    phase.decode_ns += decode_ns;
    // Publish only a complete scan whose retained capacities fit the bound.
    entries.push_front({source_, id, scan, {}, {}, calibration});
    known_bound = std::max(native_bound, bytes);
    const auto numeric_bytes = numericBytes(*scan->columns_);
    work.decoded_numeric_bytes = checkedByteAdd(work.decoded_numeric_bytes, numeric_bytes);
    phase.decoded_numeric_bytes += numeric_bytes;
    finish();
    return scan;
  }

  std::shared_ptr<const RawScan> RawAccess::materialize_(const OpenMS::MSSpectrum& native,
    OpenMS::Size index, std::optional<MZCalibration> calibration, ByteReservation reservation)
  {
    const auto native_bound = scanBytes(native);
    if (native_bound > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (native_bound > reservation.bytes()) { invalid("Scan exceeds admitted scan bytes"); }
    auto scan = std::shared_ptr<RawScan>(new RawScan);
    scan->metadata_ = metadata_source_->getSpectrumMetadata(index);
    scan->native_points_valid_ = std::all_of(native.begin(), native.end(), [](const auto& point)
    {
      return std::isfinite(point.getMZ()) && point.getMZ() >= 0 && std::isfinite(point.getIntensity());
    });
    if (native.containsIMData())
    {
      const auto [column, unit] = native.getIMData();
      scan->point_mobility_unit_ = unit;
      const auto& mobility = native.getFloatDataArrays().at(column);
      scan->native_points_valid_ = scan->native_points_valid_ && mobility.size() == native.size()
        && std::all_of(mobility.begin(), mobility.end(), [](float value) { return std::isfinite(value); });
    }
    scan->columns_ = OpenMS::OpenSwathDataAccessHelper::convertToSpectrumPtr(native);
    scan->columns_->getDataArrays().reserve(2 + native.getFloatDataArrays().size() + native.getIntegerDataArrays().size());
    appendArrays(native.getFloatDataArrays(), *scan->columns_);
    appendArrays(native.getIntegerDataArrays(), *scan->columns_);
    // Match QuadMZTransforming exactly. Never reorder or mutate a published scan.
    if (calibration)
    {
      for (double& value : scan->columns_->getMZArray()->data)
      {
        const double predict = calibration->a + calibration->b * value + calibration->c * value * value;
        value = calibration->ppm ? value - predict * value / 1000000 : predict;
      }
    }
    // Compute extraction compatibility once while materializing immutable data.
    // Auxiliary columns keep their legacy semantics and need not be peak-aligned.
    const auto& masses = scan->columns_->getMZArray()->data;
    scan->native_points_valid_ = scan->native_points_valid_ && std::all_of(masses.begin(), masses.end(),
      [](double value) { return std::isfinite(value) && value >= 0; });
    scan->masses_sorted_ = std::is_sorted(masses.begin(), masses.end());
    const auto bytes = charge(*scan->columns_);
    if (bytes > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
    if (bytes > reservation.bytes()) { invalid("Scan exceeds admitted scan bytes"); }
    scan->lease_ = reservation.commit(bytes);
    return scan;
  }

  OpenSwath::SpectrumPtr RawAccess::spectrumForScoring(int id, std::optional<MZCalibration> calibration)
  {
    const auto scan = acquire(id, calibration);
    // Do not alias the stored spectrum itself: its arrays would then own only
    // their vectors, not the scan lease, when retained without the spectrum.
    struct BorrowedSpectrum
    {
      std::shared_ptr<const RawScan> scan;
      OpenSwath::Spectrum spectrum;
    };
    auto owner = std::make_shared<BorrowedSpectrum>();
    owner->scan = scan;
    auto& arrays = owner->spectrum.getDataArrays();
    arrays.clear();
    arrays.reserve(scan->arrayCount());
    for (const auto& array : scan->columns_->getDataArrays())
    {
      arrays.emplace_back(scan, array.get());
    }
    return OpenSwath::SpectrumPtr(owner, &owner->spectrum);
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
    std::lock_guard source_lock(*source_mutex_);
    std::unique_lock cache_lock(cache_->mutex_);
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

    std::shared_ptr<RawCache::BatchFrame> batch_frame;
    if (static_cast<OpenMS::Size>(id) < batch_frames_by_id_.size())
    {
      batch_frame = batch_frames_by_id_[static_cast<OpenMS::Size>(id)].lock();
    }
    if (batch_frame)
    {
      struct StagedRegion
      {
        OpenMS::Size index = 0;
        ByteCount bytes = 0;
        RawRegion region;
        std::vector<RawCache::BatchFrame::Member> members;
      };
      std::vector<StagedRegion> staged;
      const auto decode_started = Clock::now();
      // Admit the first output before the decoder starts. Reserve each next
      // output before returning to forEachSpectrum, which decodes sequentially.
      auto admission = cache_->reserve_(cache_->max_scan_bytes_, phase_id);
      cache_lock.unlock();
      try
      {
        if (!admission) { throw BatchAdmissionFailure {}; }
        OpenMS::Size remaining = batch_frame->indices.size();
        batch_frame->source->forEachSpectrum(batch_frame->indices, [&](OpenMS::Size native_index, OpenMS::MSSpectrum&& spectrum)
        {
          const auto found_members = batch_frame->members_by_index.find(native_index);
          if (found_members == batch_frame->members_by_index.end()) { invalid("Batch source returned an unregistered index"); }
          const auto bytes = RawRegion::requiredBytes({1, spectrum.size()}, spectrum.getNativeID().size());
          if (bytes > cache_->max_scan_bytes_) { invalid("Scan exceeds maximum scan bytes"); }
          if (!admission) { throw BatchAdmissionFailure {}; }
          auto slot = std::move(*admission);
          admission.reset();
          auto region = RawRegion::fromSpectrum(spectrum, native_index, context, cache_->budget_, calibration, std::move(slot));
          staged.push_back({native_index, bytes, std::move(region), found_members->second});
          if (remaining > 0) { --remaining; }
          if (remaining > 0)
          {
            std::lock_guard lock(cache_->mutex_);
            admission = cache_->reserve_(cache_->max_scan_bytes_, phase_id);
            if (!admission)
            {
              if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
              throw BatchAdmissionFailure {};
            }
          }
        });
      }
      catch (const BatchAdmissionFailure&)
      {
        // A batch may exceed the remaining region budget. Discard staged leases,
        // then preserve the established one-scan admission behavior.
        staged.clear();
      }

      cache_lock.lock();
      if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
      if (!staged.empty())
      {
        bool current_found = false;
        RawRegion current_region;
        for (const auto& decoded : staged)
        {
          const auto scan = decoded.region.scan(0);
          const ByteCount values = scan.mz().size();
          const ByteCount numeric_bytes = checkedByteMultiply(values,
            sizeof(double) + sizeof(float) * (scan.mobility().empty() ? 1 : 2));
          ++work.decodes;
          work.decoded_numeric_bytes = checkedByteAdd(work.decoded_numeric_bytes, numeric_bytes);
          phase.decoded_numeric_bytes = checkedByteAdd(phase.decoded_numeric_bytes, numeric_bytes);
          bool published_current = false;
          for (const auto& member : decoded.members)
          {
            const auto access = member.access.lock();
            if (!access) { continue; }
            auto& member_phase = cache_->phases_[static_cast<size_t>(access->phase_(member.id))];
            ++work.misses;
            ++member_phase.misses;
            access->region_bounds_->at(static_cast<OpenMS::Size>(member.id)) = decoded.bytes;
            entries.push_front({access->source_, member.id, {}, decoded.region, context, calibration});
            if (access->source_ == source_ && member.id == id)
            {
              current_region = decoded.region;
              current_found = true;
              published_current = true;
            }
          }
          if (decoded.index == index && !published_current)
          {
            ++work.misses;
            ++phase.misses;
            region_bounds_->at(static_cast<size_t>(id)) = decoded.bytes;
            entries.push_front({source_, id, {}, decoded.region, context, calibration});
            current_region = decoded.region;
            current_found = true;
          }
        }
        const auto ns = elapsed(started);
        work.read_ns += ns;
        phase.read_ns += ns;
        const auto decode_ns = elapsed(decode_started);
        work.decode_ns += decode_ns;
        phase.decode_ns += decode_ns;
        phase.decodes += staged.size();
        if (current_found) { return current_region; }
        invalid("Batch did not produce the requested source index");
      }
    }

    ++work.misses;
    ++phase.misses;
    auto& known_bound = region_bounds_->at(static_cast<size_t>(id));
    const auto bound = known_bound ? known_bound : cache_->max_scan_bytes_;
    auto reservation = cache_->reserve_(bound, phase_id);
    if (! reservation) { invalid("Region admission failed: pinned views or cancellation"); }
    const auto decode_started = Clock::now();
    cache_lock.unlock();
    const auto native = metadata_source_->getSpectrum(index);
    const auto decode_ns = elapsed(decode_started);
    cache_lock.lock();
    if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
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
    RawPhase phase_id;
    {
      std::lock_guard source_lock(*source_mutex_);
      std::lock_guard lock(cache_->mutex_);
      phase_id = phase_(id);
    }
    std::unique_lock cache_lock(cache_->mutex_);
    const auto started = Clock::now();
    auto& work = cache_->work_;
    auto& phase = cache_->phases_[static_cast<size_t>(phase_id)];
    // Bound the detached copy from the pinned columns, not the maximum scan size.
    const auto copy_bound = checkedByteAdd(checkedByteAdd(charge(*scan->columns_), sizeof(MutableColumns)),
      checkedByteAdd(sizeof(MutableSpectrum), checkedByteMultiply(scan->arrayCount(), sizeof(OpenSwath::BinaryDataArrayPtr))));
    if (copy_bound > cache_->max_scan_bytes_) { invalid("Mutable scan exceeds maximum scan bytes"); }
    auto reservation = cache_->reserve_(copy_bound, phase_id);
    if (!reservation) { invalid(admissionFailure("Mutable scan", copy_bound, cache_->budget_.snapshot())); }
    cache_lock.unlock();
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
    cache_lock.lock();
    if (cache_->budget_.snapshot().cancelled) { invalid("Cache is cancelled"); }
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
