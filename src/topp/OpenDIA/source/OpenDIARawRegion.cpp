// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include "OpenDIARawRegion.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/IMSExperimentRandomAccess.h>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace OpenDIA
{
namespace
{
  [[noreturn]] void invalid(const char* message)
  {
    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "native region");
  }
  void validate(QueryInterval q)
  {
    if (! std::isfinite(q.lower) || ! std::isfinite(q.upper) || q.lower > q.upper
        || (q.edges != EdgePolicy::OPEN && q.edges != EdgePolicy::CLOSED && q.edges != EdgePolicy::LOWER_CLOSED
            && q.edges != EdgePolicy::UPPER_CLOSED))
    {
      invalid("Invalid interval");
    }
  }
  bool lowerClosed(QueryInterval q)
  {
    return q.edges == EdgePolicy::CLOSED || q.edges == EdgePolicy::LOWER_CLOSED;
  }
  bool upperClosed(QueryInterval q)
  {
    return q.edges == EdgePolicy::CLOSED || q.edges == EdgePolicy::UPPER_CLOSED;
  }
  bool contains(QueryInterval q, double value)
  {
    return (lowerClosed(q) ? value >= q.lower : value > q.lower) && (upperClosed(q) ? value <= q.upper : value < q.upper);
  }
  std::pair<OpenMS::Size, OpenMS::Size> rangeOffsets(std::span<const double> values, QueryInterval q)
  {
    auto first = lowerClosed(q) ? std::lower_bound(values.begin(), values.end(), q.lower) : std::upper_bound(values.begin(), values.end(), q.lower);
    auto last = upperClosed(q) ? std::upper_bound(first, values.end(), q.upper) : std::lower_bound(first, values.end(), q.upper);
    return {first - values.begin(), last - values.begin()};
  }
  void unit(MobilityUnit u)
  {
    if (u != MobilityUnit::MILLISECOND && u != MobilityUnit::INVERSE_REDUCED) { invalid("Invalid mobility unit"); }
  }
} // namespace
struct RawRegion::Storage
{
  MemoryLease lease; // Destroyed after every retained allocation.
  std::unique_ptr<RegionAcquisition[]> rows;
  std::unique_ptr<double[]> mz;
  std::unique_ptr<float[]> intensity, im;
  std::unique_ptr<char[]> ids;
  RecordContext context;
  OpenMS::Size count = 0, id_capacity = 0;
  static std::shared_ptr<Storage> allocate(RecordContext context, RecordLimits limits, OpenMS::Size id_capacity, ByteBudget& budget);
  void append(const OpenMS::MSSpectrum& spectrum, OpenMS::Size index, RecordLimits limits, std::optional<MZCalibration> calibration = {});
};
ByteCount RawRegion::requiredBytes(RecordLimits limits, OpenMS::Size id_capacity)
{
  auto bytes = checkedByteAdd(sizeof(Storage), checkedByteMultiply(limits.records, sizeof(RegionAcquisition)));
  bytes = checkedByteAdd(bytes, checkedByteMultiply(limits.values, sizeof(double) + 2 * sizeof(float)));
  return checkedByteAdd(bytes, checkedByteMultiply(limits.records, id_capacity));
}
RawRegion RawRegion::read(OpenMS::IMSExperimentRandomAccess& source,
                          std::span<const OpenMS::Size> indices,
                          RecordContext context,
                          RecordLimits limits,
                          OpenMS::Size id_capacity,
                          ByteBudget& budget,
                          std::stop_token stop)
{
  if (context.schema_version != 0 || indices.size() > limits.records) { invalid("Invalid context or acquisition limit"); }
  for (auto index : indices)
  {
    if (index >= source.getNrSpectra()) { invalid("Invalid source index"); }
  }
  if (stop.stop_requested()) { invalid("Cancelled region"); }
  auto storage = Storage::allocate(context, limits, id_capacity, budget);
  for (auto index : indices)
  {
    if (stop.stop_requested() || budget.snapshot().cancelled) { invalid("Cancelled region"); }
    auto spectrum = source.getSpectrum(index); // Decode each native acquisition exactly once.
    storage->append(spectrum, index, limits);
  }
  if (stop.stop_requested() || budget.snapshot().cancelled) { invalid("Cancelled region"); }
  RawRegion result;
  result.storage_ = std::move(storage);
  return result;
}
std::shared_ptr<RawRegion::Storage>
RawRegion::Storage::allocate(RecordContext context, RecordLimits limits, OpenMS::Size id_capacity, ByteBudget& budget)
{
  if (context.schema_version != 0) { invalid("Invalid context"); }
  auto reservation = budget.tryReserve(requiredBytes(limits, id_capacity));
  if (! reservation) { invalid("Region admission failed"); }
  auto storage = std::make_shared<Storage>();
  storage->lease = reservation->commit(reservation->bytes());
  storage->rows = std::make_unique<RegionAcquisition[]>(limits.records);
  storage->mz = std::make_unique<double[]>(limits.values);
  storage->intensity = std::make_unique<float[]>(limits.values);
  storage->im = std::make_unique<float[]>(limits.values);
  storage->ids = std::make_unique<char[]>(checkedByteMultiply(limits.records, id_capacity));
  storage->context = context;
  storage->id_capacity = id_capacity;
  return storage;
}
void RawRegion::Storage::append(const OpenMS::MSSpectrum& spectrum, OpenMS::Size index, RecordLimits limits, std::optional<MZCalibration> calibration)
{
  auto* storage = this;
  OpenMS::Size offset = count ? rows[count - 1].end : 0;
  if (count >= limits.records) { invalid("Region acquisition capacity exceeded"); }
  const auto calibrated = [&](OpenMS::Size p) {
    const double mz = spectrum[p].getMZ();
    if (! calibration) { return mz; }
    const auto& [a, b, c, ppm] = *calibration;
    const double predict = a + b * mz + c * mz * mz;
    return ppm ? mz - predict * mz / 1000000 : predict;
  };
  if (spectrum.size() > limits.values - offset || spectrum.getNativeID().size() > id_capacity) { invalid("Region capacity exceeded"); }
  auto& row = storage->rows[storage->count];
  row.source_index = index;
  row.rt_seconds = spectrum.getRT();
  row.ms_level = spectrum.getMSLevel();
  if (! std::isfinite(row.rt_seconds) || row.ms_level == 0) { invalid("Invalid acquisition metadata"); }
  row.native_id_size = spectrum.getNativeID().size();
  std::copy(spectrum.getNativeID().begin(), spectrum.getNativeID().end(), storage->ids.get() + storage->count * id_capacity);
  if (spectrum.getDriftTimeUnit() == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
  {
    if (! std::isfinite(spectrum.getDriftTime())) { invalid("Invalid FAIMS"); }
    row.faims_voltage = spectrum.getDriftTime();
  }
  else if (spectrum.getDriftTimeUnit() != OpenMS::DriftTimeUnit::NONE)
  {
    if (! std::isfinite(spectrum.getDriftTime())) { invalid("Invalid scalar mobility"); }
    if (spectrum.getDriftTimeUnit() == OpenMS::DriftTimeUnit::MILLISECOND)
    {
      row.acquisition_mobility = Mobility {spectrum.getDriftTime(), MobilityUnit::MILLISECOND};
    }
    else if (spectrum.getDriftTimeUnit() == OpenMS::DriftTimeUnit::VSSC)
    {
      row.acquisition_mobility = Mobility {spectrum.getDriftTime(), MobilityUnit::INVERSE_REDUCED};
    }
    else { invalid("Unsupported scalar mobility unit"); }
  }
  const OpenMS::MSSpectrum::FloatDataArray* im = nullptr;
  if (spectrum.containsIMData())
  {
    const auto [column, native_unit] = spectrum.getIMData();
    if (native_unit == OpenMS::DriftTimeUnit::MILLISECOND) { row.point_mobility = MobilityUnit::MILLISECOND; }
    else if (native_unit == OpenMS::DriftTimeUnit::VSSC) { row.point_mobility = MobilityUnit::INVERSE_REDUCED; }
    else { invalid("Unsupported point mobility unit"); }
    im = &spectrum.getFloatDataArrays().at(column);
    if (im->size() != spectrum.size()) { invalid("Mismatched point mobility column"); }
  }
  // Validate calibrated order before sorting. Sorted results need no permutation workspace.
  bool sorted = true;
  for (OpenMS::Size p = 0; p < spectrum.size(); ++p)
  {
    const double mz = calibrated(p);
    if (p && calibrated(p - 1) > mz) { sorted = false; }
    if (! std::isfinite(mz) || mz < 0 || ! std::isfinite(spectrum[p].getMZ()) || spectrum[p].getMZ() < 0
        || ! std::isfinite(spectrum[p].getIntensity()) || (im && ! std::isfinite((*im)[p])))
    {
      invalid("Invalid native point");
    }
  }
  row.begin = offset;
  const auto append = [&](OpenMS::Size p) {
    storage->mz[offset] = calibrated(p);
    storage->intensity[offset] = spectrum[p].getIntensity();
    if (im) { storage->im[offset] = (*im)[p]; }
    ++offset;
  };
  if (sorted)
  {
    for (OpenMS::Size p = 0; p < spectrum.size(); ++p)
    {
      append(p);
    }
  }
  else
  {
    // Stable-sort/permutation workspace is transient native allowance, outside the retained cap.
    std::vector<OpenMS::Size> order(spectrum.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return calibrated(a) < calibrated(b); });
    for (auto p : order)
    {
      append(p);
    }
  }
  row.end = offset;
  ++storage->count;
}
RawRegion RawRegion::fromSpectrum(const OpenMS::MSSpectrum& spectrum,
                                  OpenMS::Size source_index,
                                  RecordContext context,
                                  ByteBudget& budget,
                                  std::optional<MZCalibration> calibration)
{
  if (calibration && (! std::isfinite(calibration->a) || ! std::isfinite(calibration->b) || ! std::isfinite(calibration->c)))
  {
    invalid("Invalid calibration coefficients");
  }
  const RecordLimits limits {1, spectrum.size()};
  auto storage = Storage::allocate(context, limits, spectrum.getNativeID().size(), budget);
  storage->append(spectrum, source_index, limits, calibration);
  if (budget.snapshot().cancelled) { invalid("Cancelled region"); }
  RawRegion result;
  result.storage_ = std::move(storage);
  return result;
}
OpenMS::Size RawRegion::size() const
{
  if (! storage_) { invalid("Empty region handle"); }
  return storage_->count;
}
RawRegion::Scan RawRegion::scan(OpenMS::Size index) const
{
  if (index >= size()) { invalid("Invalid acquisition row"); }
  Scan result;
  result.storage_ = storage_;
  result.row_ = index;
  result.begin_ = storage_->rows[index].begin;
  result.end_ = storage_->rows[index].end;
  return result;
}
const RegionAcquisition& RawRegion::Scan::acquisition() const
{
  if (! storage_) { invalid("Empty scan handle"); }
  return storage_->rows[row_];
}
RecordContext RawRegion::Scan::context() const
{
  acquisition();
  return storage_->context;
}
std::string_view RawRegion::Scan::nativeID() const
{
  const auto length = acquisition().native_id_size;
  return {storage_->ids.get() + row_ * storage_->id_capacity, length};
}
std::span<const double> RawRegion::Scan::mz() const
{
  acquisition();
  return {storage_->mz.get() + begin_, end_ - begin_};
}
std::span<const float> RawRegion::Scan::intensity() const
{
  acquisition();
  return {storage_->intensity.get() + begin_, end_ - begin_};
}
std::span<const float> RawRegion::Scan::mobility() const
{
  if (! acquisition().point_mobility) { return {}; }
  return {storage_->im.get() + begin_, end_ - begin_};
}
RawRegion::Scan RawRegion::Scan::range(QueryInterval q) const
{
  validate(q);
  const auto [first, last] = rangeOffsets(mz(), q);
  Scan result = *this;
  result.begin_ += first;
  result.end_ = begin_ + last;
  return result;
}
void RawRegion::sum(std::span<const RegionSum> requests, std::span<double> output) const
{
  const auto count = size();
  if (requests.size() != output.size()) { invalid("Output span size mismatch"); }
  for (const auto& request : requests)
  {
    if (request.acquisition >= count) { invalid("Invalid acquisition row"); }
    const auto& row = storage_->rows[request.acquisition];
    validate(request.mz);
    if (request.mobility)
    {
      validate(request.mobility->interval);
      unit(request.mobility->unit);
      if (! row.point_mobility) { invalid("Point mobility is absent"); }
      if (*row.point_mobility != request.mobility->unit) { invalid("Mobility unit mismatch"); }
    }
  }
  // Direct-row kernel: the region already pins storage, no per-query lease copies.
  for (OpenMS::Size q = 0; q < requests.size(); ++q)
  {
    const auto& request = requests[q];
    const auto& row = storage_->rows[request.acquisition];
    const auto [first, last] = rangeOffsets({storage_->mz.get() + row.begin, row.end - row.begin}, request.mz);
    double total = 0;
    for (OpenMS::Size p = row.begin + first; p < row.begin + last; ++p)
    {
      if (! request.mobility || contains(request.mobility->interval, storage_->im[p])) { total += storage_->intensity[p]; }
    }
    output[q] = total;
  }
}
} // namespace OpenDIA
