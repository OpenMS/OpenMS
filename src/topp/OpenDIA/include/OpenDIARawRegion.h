// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#pragma once
#include "OpenDIARecords.h"

#include <span>
#include <string_view>

namespace OpenMS
{
class IMSExperimentRandomAccess;
class MSSpectrum;
} // namespace OpenMS
namespace OpenDIA
{
/// Tool-private quadratic calibration, matching the OpenSWATH transform.
struct MZCalibration
{
  double a = 0, b = 1, c = 0;
  bool ppm = false;
  bool operator==(const MZCalibration&) const = default;
};
struct RegionAcquisition
{
  RecordId source_index = 0;
  double rt_seconds = 0;
  unsigned ms_level = 0;
  std::optional<MobilityUnit> point_mobility;
  /// Exact scalar acquisition annotation, never substituted for point IM.
  std::optional<Mobility> acquisition_mobility;
  std::optional<double> faims_voltage;
  /// CSR offsets in the full region, including for a narrowed Scan view.
  OpenMS::Size begin = 0, end = 0, native_id_size = 0;
};
struct RegionSum
{
  /// Region row index, not a native source index. Source identity is never renumbered.
  OpenMS::Size acquisition = 0;
  QueryInterval mz;
  std::optional<MobilityInterval> mobility;
};

/** @brief Immutable, bounded native acquisition/point columns, independent of OpenSWATH.
  Caller selects source indices and lanes using full native isolation geometry.
  This is NOT a generic isolation/region selector or a native-frame backend.
  Each requested index produces one row, in request order (including repeats).
  Unsorted native peaks are stably sorted, retaining duplicates and point IM.
  Scalar drift-time annotations are retained separately from point mobility.
  Unsupported scalar units (e.g. CCS) are rejected, never silently discarded.
  Other native metadata and isolation geometry remain on the source.
  Decode/sort workspace and allocator/control overhead use the native allowance,
  not an RSS guarantee. Fixed arrays allocate exactly limits.records rows,
  limits.values doubles plus two floats per point, and records * id_capacity
  chars. All capacities are admitted before allocation, never grown.
*/
class RawRegion
{
  struct Storage;

public:
  /// Owning scan/range view. Spans and native ID borrow this view's lifetime.
  class Scan
  {
  public:
    const RegionAcquisition& acquisition() const;
    RecordContext context() const;
    std::string_view nativeID() const;
    std::span<const double> mz() const;
    std::span<const float> intensity() const;
    /// Empty when mobility is absent, even for nonempty peaks.
    std::span<const float> mobility() const;
    /// Binary range lookup, no allocation. Returned view retains the whole lease.
    Scan range(QueryInterval interval) const;

  private:
    friend class RawRegion;
    std::shared_ptr<const Storage> storage_;
    OpenMS::Size row_ = 0, begin_ = 0, end_ = 0;
  };
  /** @brief Decode explicit native indices, failing before publishing any storage.
    @param[in] source Exclusive, unchanged native reader for the duration of this call.
    @param[in] indices Native source indices in deterministic caller order.
    @param[in] context Caller namespace, schema version must be zero.
    @param[in] limits Fixed maximum allocated acquisition and point capacities.
    @param[in] id_capacity Maximum bytes per native ID, no terminator needed.
    @param[in] budget Retained allocation admission (nonblocking).
    @param[in] stop Cancellation checked during decode and before publication.
  */
  static RawRegion read(OpenMS::IMSExperimentRandomAccess& source,
                        std::span<const OpenMS::Size> indices,
                        RecordContext context,
                        RecordLimits limits,
                        OpenMS::Size id_capacity,
                        ByteBudget& budget,
                        std::stop_token stop = {});
  /// Copy one decoded spectrum into exact-size admitted storage, without changing it.
  static RawRegion fromSpectrum(const OpenMS::MSSpectrum& spectrum,
                                OpenMS::Size source_index,
                                RecordContext context,
                                ByteBudget& budget,
                                std::optional<MZCalibration> calibration = {});
  static ByteCount requiredBytes(RecordLimits limits, OpenMS::Size id_capacity);
  OpenMS::Size size() const;
  Scan scan(OpenMS::Size index) const;
  /// Validate all requests before writing. Output size must equal request size.
  /// Missing point IM and incompatible units throw before any output is written.
  /// Requests need not be sorted. No matrix, cache, or hot-path allocations.
  void sum(std::span<const RegionSum> requests, std::span<double> output) const;

private:
  std::shared_ptr<const Storage> storage_;
};
} // namespace OpenDIA
