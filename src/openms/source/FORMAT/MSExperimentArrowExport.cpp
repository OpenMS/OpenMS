// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/c/bridge.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <cstdint>

namespace OpenMS
{

namespace // anonymous namespace for helper functions
{

/// Check if a column should be included based on user selection
bool shouldIncludeColumn(const std::string& col_name,
                         const std::vector<std::string>& requested_columns)
{
  if (requested_columns.empty())
  {
    return true; // Include all if no filter specified
  }
  return std::find(requested_columns.begin(), requested_columns.end(), col_name)
         != requested_columns.end();
}

/// Check if spectrum passes MS level filter only (RT filtering is done via binary search)
bool passesMSLevelFilter(const MSSpectrum& spec,
                         const std::unordered_set<UInt>& ms_levels_set)
{
  if (ms_levels_set.empty()) return true;
  return ms_levels_set.find(spec.getMSLevel()) != ms_levels_set.end();
}

/// Get iterator range for RT-filtered spectra using binary search
/// Returns pair of (begin, end) iterators for spectra in RT range
std::pair<MSExperiment::ConstIterator, MSExperiment::ConstIterator>
getRTFilteredRange(const MSExperiment& exp, const ArrowSpectraExportConfig& config)
{
  MSExperiment::ConstIterator begin_it = exp.begin();
  MSExperiment::ConstIterator end_it = exp.end();

  // Use binary search for RT filtering if specified
  if (config.min_rt != 0)
  {
    begin_it = exp.RTBegin(config.min_rt);
  }
  if (config.max_rt != 0)
  {
    end_it = exp.RTEnd(config.max_rt);
  }

  return {begin_it, end_it};
}

/// Get iterator range for m/z-filtered peaks using binary search
/// Returns pair of (begin, end) iterators for peaks in m/z range
std::pair<MSSpectrum::ConstIterator, MSSpectrum::ConstIterator>
getMZFilteredRange(const MSSpectrum& spec, double min_mz, double max_mz)
{
  MSSpectrum::ConstIterator begin_it = spec.begin();
  MSSpectrum::ConstIterator end_it = spec.end();

  // Use binary search for m/z filtering if specified
  if (min_mz != 0)
  {
    begin_it = spec.MZBegin(min_mz);
  }
  if (max_mz != 0)
  {
    end_it = spec.MZEnd(max_mz);
  }

  return {begin_it, end_it};
}

/// Check if any spectrum in the experiment has ion mobility data
bool experimentHasIMData(const MSExperiment& exp,
                         const std::unordered_set<UInt>& ms_levels_set,
                         const ArrowSpectraExportConfig& config)
{
  auto [begin_it, end_it] = getRTFilteredRange(exp, config);
  for (auto it = begin_it; it != end_it; ++it)
  {
    if (!passesMSLevelFilter(*it, ms_levels_set)) continue;
    if (it->containsIMData()) return true;
  }
  return false;
}

/// Count total peaks that will be exported (for capacity reservation)
Size countTotalPeaks(const MSExperiment& exp,
                     const ArrowSpectraExportConfig& config,
                     const std::unordered_set<UInt>& ms_levels_set)
{
  Size total = 0;
  bool filter_mz = (config.min_mz != 0 || config.max_mz != 0);

  auto [begin_it, end_it] = getRTFilteredRange(exp, config);
  for (auto it = begin_it; it != end_it; ++it)
  {
    if (!passesMSLevelFilter(*it, ms_levels_set)) continue;

    if (filter_mz)
    {
      // Use binary search to count peaks in m/z range
      auto [mz_begin, mz_end] = getMZFilteredRange(*it, config.min_mz, config.max_mz);
      total += std::distance(mz_begin, mz_end);
    }
    else
    {
      total += it->size();
    }
  }
  return total;
}

/// Count spectra that will be exported
Size countFilteredSpectra(const MSExperiment& exp,
                          const ArrowSpectraExportConfig& config,
                          const std::unordered_set<UInt>& ms_levels_set)
{
  Size count = 0;
  auto [begin_it, end_it] = getRTFilteredRange(exp, config);
  for (auto it = begin_it; it != end_it; ++it)
  {
    if (passesMSLevelFilter(*it, ms_levels_set)) ++count;
  }
  return count;
}

/// Build long format Arrow table
std::shared_ptr<arrow::Table> buildLongFormatTable(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  const std::unordered_set<UInt>& ms_levels_set)
{
  // Determine which columns to include
  const auto& cols = config.columns;
  bool inc_mz = shouldIncludeColumn(SpectraLongSchema::MZ, cols);
  bool inc_intensity = shouldIncludeColumn(SpectraLongSchema::INTENSITY, cols);
  bool inc_rt = shouldIncludeColumn(SpectraLongSchema::RT, cols);
  bool inc_spectrum_index = shouldIncludeColumn(SpectraLongSchema::SPECTRUM_INDEX, cols);
  bool inc_ms_level = shouldIncludeColumn(SpectraLongSchema::MS_LEVEL, cols);
  bool inc_native_id = shouldIncludeColumn(SpectraLongSchema::NATIVE_ID, cols);
  bool inc_ion_mobility = config.include_ion_mobility && shouldIncludeColumn(SpectraLongSchema::ION_MOBILITY, cols);
  bool inc_precursor_mz = config.include_precursor_info && shouldIncludeColumn(SpectraLongSchema::PRECURSOR_MZ, cols);
  bool inc_precursor_charge = config.include_precursor_info && shouldIncludeColumn(SpectraLongSchema::PRECURSOR_CHARGE, cols);
  bool inc_precursor_intensity = config.include_precursor_info && shouldIncludeColumn(SpectraLongSchema::PRECURSOR_INTENSITY, cols);
  bool inc_isolation_lower = config.include_precursor_info && shouldIncludeColumn(SpectraLongSchema::ISOLATION_LOWER, cols);
  bool inc_isolation_upper = config.include_precursor_info && shouldIncludeColumn(SpectraLongSchema::ISOLATION_UPPER, cols);

  // Check if experiment has IM data (only if we want to include it)
  bool has_im_data = inc_ion_mobility && experimentHasIMData(exp, ms_levels_set, config);

  // Pre-count total peaks for capacity reservation
  Size total_peaks = countTotalPeaks(exp, config, ms_levels_set);

  // Create Arrow memory pool
  arrow::MemoryPool* pool = arrow::default_memory_pool();

  // Create builders with reserved capacity
  arrow::DoubleBuilder mz_builder(pool);
  arrow::FloatBuilder intensity_builder(pool);
  arrow::FloatBuilder rt_builder(pool);
  arrow::FloatBuilder ion_mobility_builder(pool);
  arrow::UInt32Builder spectrum_index_builder(pool);
  arrow::UInt8Builder ms_level_builder(pool);
  arrow::StringBuilder native_id_builder(pool);
  arrow::DoubleBuilder precursor_mz_builder(pool);
  arrow::Int16Builder precursor_charge_builder(pool);
  arrow::FloatBuilder precursor_intensity_builder(pool);
  arrow::DoubleBuilder isolation_lower_builder(pool);
  arrow::DoubleBuilder isolation_upper_builder(pool);

  // Reserve capacity
  arrow::Status status;
  if (inc_mz) { status = mz_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_intensity) { status = intensity_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow intensity_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_rt) { status = rt_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_spectrum_index) { status = spectrum_index_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow spectrum_index_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_ms_level) { status = ms_level_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ms_level_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (has_im_data) { status = ion_mobility_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_mz) { status = precursor_mz_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_charge) { status = precursor_charge_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_charge_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_intensity) { status = precursor_intensity_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_intensity_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_isolation_lower) { status = isolation_lower_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_lower_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_isolation_upper) { status = isolation_upper_builder.Reserve(total_peaks); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_upper_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; } }

  // Use binary search for RT range filtering
  auto [rt_begin, rt_end] = getRTFilteredRange(exp, config);

  // Iterate through RT-filtered spectra
  UInt32 spectrum_idx = static_cast<UInt32>(std::distance(exp.begin(), rt_begin));
  for (auto spec_it = rt_begin; spec_it != rt_end; ++spec_it, ++spectrum_idx)
  {
    const auto& spec = *spec_it;

    // Check MS level filter only (RT already filtered by binary search)
    if (!passesMSLevelFilter(spec, ms_levels_set))
    {
      continue;
    }

    // Get spectrum-level data once per spectrum
    float rt = static_cast<float>(spec.getRT());
    uint8_t ms_level = static_cast<uint8_t>(spec.getMSLevel());
    const String& native_id = spec.getNativeID();

    // Precursor info (null for MS1)
    bool has_precursor = !spec.getPrecursors().empty();
    double prec_mz = has_precursor ? spec.getPrecursors()[0].getMZ() : 0.0;
    int16_t prec_charge = has_precursor ? static_cast<int16_t>(spec.getPrecursors()[0].getCharge()) : 0;
    float prec_intensity = has_precursor ? static_cast<float>(spec.getPrecursors()[0].getIntensity()) : 0.0f;
    double iso_lower = has_precursor ? spec.getPrecursors()[0].getIsolationWindowLowerOffset() : 0.0;
    double iso_upper = has_precursor ? spec.getPrecursors()[0].getIsolationWindowUpperOffset() : 0.0;

    // Get ion mobility data if present
    const float* im_data_ptr = nullptr;
    if (has_im_data && spec.containsIMData())
    {
      auto [im_idx, im_unit] = spec.getIMData();
      const auto& im_array = spec.getFloatDataArrays()[im_idx];
      im_data_ptr = im_array.data();
    }

    // Use binary search for m/z range filtering
    auto [mz_begin, mz_end] = getMZFilteredRange(spec, config.min_mz, config.max_mz);

    // Iterate through m/z-filtered peaks
    for (auto peak_it = mz_begin; peak_it != mz_end; ++peak_it)
    {
      double mz = peak_it->getMZ();
      // Calculate peak index for ion mobility array access
      Size peak_idx = static_cast<Size>(std::distance(spec.begin(), peak_it));

      // Append peak data using UnsafeAppend (we already reserved capacity)
      if (inc_mz) mz_builder.UnsafeAppend(mz);
      if (inc_intensity) intensity_builder.UnsafeAppend(static_cast<float>(peak_it->getIntensity()));
      if (inc_rt) rt_builder.UnsafeAppend(rt);
      if (inc_spectrum_index) spectrum_index_builder.UnsafeAppend(spectrum_idx);
      if (inc_ms_level) ms_level_builder.UnsafeAppend(ms_level);

      // native_id requires safe append (string)
      if (inc_native_id)
      {
        status = native_id_builder.Append(native_id);
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow native_id_builder Append failed: " << status.ToString() << std::endl; return nullptr; }
      }

      // Ion mobility
      if (has_im_data)
      {
        if (im_data_ptr)
        {
          ion_mobility_builder.UnsafeAppend(im_data_ptr[peak_idx]);
        }
        else
        {
          status = ion_mobility_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ion_mobility_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      // Precursor info
      if (inc_precursor_mz)
      {
        if (has_precursor)
        {
          precursor_mz_builder.UnsafeAppend(prec_mz);
        }
        else
        {
          status = precursor_mz_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_mz_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      if (inc_precursor_charge)
      {
        if (has_precursor)
        {
          precursor_charge_builder.UnsafeAppend(prec_charge);
        }
        else
        {
          status = precursor_charge_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_charge_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      if (inc_precursor_intensity)
      {
        if (has_precursor)
        {
          precursor_intensity_builder.UnsafeAppend(prec_intensity);
        }
        else
        {
          status = precursor_intensity_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_intensity_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      if (inc_isolation_lower)
      {
        if (has_precursor)
        {
          isolation_lower_builder.UnsafeAppend(iso_lower);
        }
        else
        {
          status = isolation_lower_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_lower_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      if (inc_isolation_upper)
      {
        if (has_precursor)
        {
          isolation_upper_builder.UnsafeAppend(iso_upper);
        }
        else
        {
          status = isolation_upper_builder.AppendNull();
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_upper_builder AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }
    }
  }

  // Build schema and arrays
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  std::shared_ptr<arrow::Array> arr;

  if (inc_mz)
  {
    status = mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::MZ, arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_intensity)
  {
    status = intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow intensity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::INTENSITY, arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_rt)
  {
    status = rt_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::RT, arrow::float32()));
    arrays.push_back(arr);
  }

  if (has_im_data)
  {
    status = ion_mobility_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::ION_MOBILITY, arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_spectrum_index)
  {
    status = spectrum_index_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow spectrum_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::SPECTRUM_INDEX, arrow::uint32()));
    arrays.push_back(arr);
  }

  if (inc_ms_level)
  {
    status = ms_level_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ms_level_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::MS_LEVEL, arrow::uint8()));
    arrays.push_back(arr);
  }

  if (inc_native_id)
  {
    status = native_id_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow native_id_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::NATIVE_ID, arrow::utf8()));
    arrays.push_back(arr);
  }

  if (inc_precursor_mz)
  {
    status = precursor_mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::PRECURSOR_MZ, arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_precursor_charge)
  {
    status = precursor_charge_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::PRECURSOR_CHARGE, arrow::int16()));
    arrays.push_back(arr);
  }

  if (inc_precursor_intensity)
  {
    status = precursor_intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_intensity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::PRECURSOR_INTENSITY, arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_isolation_lower)
  {
    status = isolation_lower_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_lower_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::ISOLATION_LOWER, arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_isolation_upper)
  {
    status = isolation_upper_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_upper_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraLongSchema::ISOLATION_UPPER, arrow::float64()));
    arrays.push_back(arr);
  }

  auto schema = arrow::schema(fields);
  auto table = arrow::Table::Make(schema, arrays);

  // Validate table against registry schema (subset — dynamic columns are a subset of the superset)
  auto validation = ArrowSchemaValidation::validate(table, SpectraLongSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Spectra long format schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


/// Build semi-wide format Arrow table (one row per spectrum)
std::shared_ptr<arrow::Table> buildSemiWideFormatTable(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  const std::unordered_set<UInt>& ms_levels_set)
{
  // Determine which columns to include
  const auto& cols = config.columns;
  bool inc_mz = shouldIncludeColumn(SpectraSemiWideSchema::MZ, cols);
  bool inc_intensity = shouldIncludeColumn(SpectraSemiWideSchema::INTENSITY, cols);
  bool inc_rt = shouldIncludeColumn(SpectraSemiWideSchema::RT, cols);
  bool inc_spectrum_index = shouldIncludeColumn(SpectraSemiWideSchema::SPECTRUM_INDEX, cols);
  bool inc_ms_level = shouldIncludeColumn(SpectraSemiWideSchema::MS_LEVEL, cols);
  bool inc_native_id = shouldIncludeColumn(SpectraSemiWideSchema::NATIVE_ID, cols);
  bool inc_ion_mobility = config.include_ion_mobility && shouldIncludeColumn(SpectraSemiWideSchema::ION_MOBILITY, cols);
  bool inc_precursor_mz = config.include_precursor_info && shouldIncludeColumn(SpectraSemiWideSchema::PRECURSOR_MZ, cols);
  bool inc_precursor_charge = config.include_precursor_info && shouldIncludeColumn(SpectraSemiWideSchema::PRECURSOR_CHARGE, cols);
  bool inc_precursor_intensity = config.include_precursor_info && shouldIncludeColumn(SpectraSemiWideSchema::PRECURSOR_INTENSITY, cols);
  bool inc_isolation_lower = config.include_precursor_info && shouldIncludeColumn(SpectraSemiWideSchema::ISOLATION_LOWER, cols);
  bool inc_isolation_upper = config.include_precursor_info && shouldIncludeColumn(SpectraSemiWideSchema::ISOLATION_UPPER, cols);

  // Check if experiment has IM data
  bool has_im_data = inc_ion_mobility && experimentHasIMData(exp, ms_levels_set, config);

  // Count spectra for capacity reservation
  Size num_spectra = countFilteredSpectra(exp, config, ms_levels_set);

  // Create Arrow memory pool
  arrow::MemoryPool* pool = arrow::default_memory_pool();

  // Scalar column builders
  arrow::UInt32Builder spectrum_index_builder(pool);
  arrow::FloatBuilder rt_builder(pool);
  arrow::UInt8Builder ms_level_builder(pool);
  arrow::StringBuilder native_id_builder(pool);
  arrow::DoubleBuilder precursor_mz_builder(pool);
  arrow::Int16Builder precursor_charge_builder(pool);
  arrow::FloatBuilder precursor_intensity_builder(pool);
  arrow::DoubleBuilder isolation_lower_builder(pool);
  arrow::DoubleBuilder isolation_upper_builder(pool);

  // List builders for peak data
  auto mz_value_builder = std::make_shared<arrow::DoubleBuilder>(pool);
  arrow::ListBuilder mz_list_builder(pool, mz_value_builder);

  auto intensity_value_builder = std::make_shared<arrow::FloatBuilder>(pool);
  arrow::ListBuilder intensity_list_builder(pool, intensity_value_builder);

  auto im_value_builder = std::make_shared<arrow::FloatBuilder>(pool);
  arrow::ListBuilder im_list_builder(pool, im_value_builder);

  // Reserve capacity for scalar columns
  arrow::Status status;
  if (inc_spectrum_index) { status = spectrum_index_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_rt) { status = rt_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_ms_level) { status = ms_level_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_mz) { status = precursor_mz_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_charge) { status = precursor_charge_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_precursor_intensity) { status = precursor_intensity_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_isolation_lower) { status = isolation_lower_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
  if (inc_isolation_upper) { status = isolation_upper_builder.Reserve(num_spectra); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }

  // Use binary search for RT range filtering
  auto [rt_begin, rt_end] = getRTFilteredRange(exp, config);

  // Iterate through RT-filtered spectra
  UInt32 spectrum_idx = static_cast<UInt32>(std::distance(exp.begin(), rt_begin));
  for (auto spec_it = rt_begin; spec_it != rt_end; ++spec_it, ++spectrum_idx)
  {
    const auto& spec = *spec_it;

    // Check MS level filter only (RT already filtered by binary search)
    if (!passesMSLevelFilter(spec, ms_levels_set))
    {
      continue;
    }

    // Append scalar columns
    if (inc_spectrum_index) spectrum_index_builder.UnsafeAppend(spectrum_idx);
    if (inc_rt) rt_builder.UnsafeAppend(static_cast<float>(spec.getRT()));
    if (inc_ms_level) ms_level_builder.UnsafeAppend(static_cast<uint8_t>(spec.getMSLevel()));

    if (inc_native_id)
    {
      status = native_id_builder.Append(spec.getNativeID());
      if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Append failed: " << status.ToString() << std::endl; return nullptr; }
    }

    // Precursor info
    bool has_precursor = !spec.getPrecursors().empty();

    if (inc_precursor_mz)
    {
      if (has_precursor) precursor_mz_builder.UnsafeAppend(spec.getPrecursors()[0].getMZ());
      else { status = precursor_mz_builder.AppendNull(); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; } }
    }

    if (inc_precursor_charge)
    {
      if (has_precursor) precursor_charge_builder.UnsafeAppend(static_cast<int16_t>(spec.getPrecursors()[0].getCharge()));
      else { status = precursor_charge_builder.AppendNull(); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; } }
    }

    if (inc_precursor_intensity)
    {
      if (has_precursor) precursor_intensity_builder.UnsafeAppend(static_cast<float>(spec.getPrecursors()[0].getIntensity()));
      else { status = precursor_intensity_builder.AppendNull(); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; } }
    }

    if (inc_isolation_lower)
    {
      if (has_precursor) isolation_lower_builder.UnsafeAppend(spec.getPrecursors()[0].getIsolationWindowLowerOffset());
      else { status = isolation_lower_builder.AppendNull(); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; } }
    }

    if (inc_isolation_upper)
    {
      if (has_precursor) isolation_upper_builder.UnsafeAppend(spec.getPrecursors()[0].getIsolationWindowUpperOffset());
      else { status = isolation_upper_builder.AppendNull(); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; } }
    }

    // Get ion mobility data if present
    const float* im_data_ptr = nullptr;
    if (has_im_data && spec.containsIMData())
    {
      auto [im_idx, im_unit] = spec.getIMData();
      const auto& im_array = spec.getFloatDataArrays()[im_idx];
      im_data_ptr = im_array.data();
    }

    // Use binary search to get m/z-filtered peak range
    auto [mz_begin, mz_end] = getMZFilteredRange(spec, config.min_mz, config.max_mz);

    // Append list columns (mz, intensity, ion_mobility arrays)
    if (inc_mz)
    {
      status = mz_list_builder.Append();
      if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

      for (auto peak_it = mz_begin; peak_it != mz_end; ++peak_it)
      {
        status = mz_value_builder->Append(peak_it->getMZ());
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
      }
    }

    if (inc_intensity)
    {
      status = intensity_list_builder.Append();
      if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

      for (auto peak_it = mz_begin; peak_it != mz_end; ++peak_it)
      {
        status = intensity_value_builder->Append(static_cast<float>(peak_it->getIntensity()));
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
      }
    }

    if (has_im_data)
    {
      if (im_data_ptr)
      {
        status = im_list_builder.Append();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

        for (auto peak_it = mz_begin; peak_it != mz_end; ++peak_it)
        {
          // Calculate peak index for ion mobility array access
          Size peak_idx = static_cast<Size>(std::distance(spec.begin(), peak_it));
          status = im_value_builder->Append(im_data_ptr[peak_idx]);
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }
      else
      {
        status = im_list_builder.AppendNull();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
      }
    }
  }

  // Build schema and arrays
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  std::shared_ptr<arrow::Array> arr;

  if (inc_spectrum_index)
  {
    status = spectrum_index_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::SPECTRUM_INDEX, arrow::uint32()));
    arrays.push_back(arr);
  }

  if (inc_rt)
  {
    status = rt_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::RT, arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_ms_level)
  {
    status = ms_level_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::MS_LEVEL, arrow::uint8()));
    arrays.push_back(arr);
  }

  if (inc_native_id)
  {
    status = native_id_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::NATIVE_ID, arrow::utf8()));
    arrays.push_back(arr);
  }

  if (inc_mz)
  {
    status = mz_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::MZ, arrow::list(arrow::float64())));
    arrays.push_back(arr);
  }

  if (inc_intensity)
  {
    status = intensity_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::INTENSITY, arrow::list(arrow::float32())));
    arrays.push_back(arr);
  }

  if (has_im_data)
  {
    status = im_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::ION_MOBILITY, arrow::list(arrow::float32())));
    arrays.push_back(arr);
  }

  if (inc_precursor_mz)
  {
    status = precursor_mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::PRECURSOR_MZ, arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_precursor_charge)
  {
    status = precursor_charge_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::PRECURSOR_CHARGE, arrow::int16()));
    arrays.push_back(arr);
  }

  if (inc_precursor_intensity)
  {
    status = precursor_intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::PRECURSOR_INTENSITY, arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_isolation_lower)
  {
    status = isolation_lower_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::ISOLATION_LOWER, arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_isolation_upper)
  {
    status = isolation_upper_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field(SpectraSemiWideSchema::ISOLATION_UPPER, arrow::float64()));
    arrays.push_back(arr);
  }

  auto schema = arrow::schema(fields);
  auto table = arrow::Table::Make(schema, arrays);

  // Validate table against registry schema (subset — dynamic columns are a subset of the superset)
  auto validation = ArrowSchemaValidation::validate(table, SpectraSemiWideSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Spectra semi-wide format schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

// Internal function - not part of public API
std::shared_ptr<arrow::Table> exportSpectraToArrow(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config)
{
  // Convert ms_levels vector to set for O(1) lookup
  std::unordered_set<UInt> ms_levels_set(config.ms_levels.begin(), config.ms_levels.end());

  // Dispatch to appropriate format builder
  if (config.format == ArrowExportFormat::Long)
  {
    return buildLongFormatTable(exp, config, ms_levels_set);
  }
  else
  {
    return buildSemiWideFormatTable(exp, config, ms_levels_set);
  }
}

} // anonymous namespace


std::vector<std::string> MSExperimentArrowExport::getSpectraArrowColumnNames(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config)
{
  std::vector<std::string> columns;

  // Convert ms_levels vector to set for O(1) lookup
  std::unordered_set<UInt> ms_levels_set(config.ms_levels.begin(), config.ms_levels.end());

  // Check if experiment has IM data
  bool has_im_data = config.include_ion_mobility && experimentHasIMData(exp, ms_levels_set, config);

  if (config.format == ArrowExportFormat::Long)
  {
    columns = {SpectraLongSchema::MZ, SpectraLongSchema::INTENSITY, SpectraLongSchema::RT};
    if (has_im_data) columns.push_back(SpectraLongSchema::ION_MOBILITY);
    columns.push_back(SpectraLongSchema::SPECTRUM_INDEX);
    columns.push_back(SpectraLongSchema::MS_LEVEL);
    columns.push_back(SpectraLongSchema::NATIVE_ID);

    if (config.include_precursor_info)
    {
      columns.push_back(SpectraLongSchema::PRECURSOR_MZ);
      columns.push_back(SpectraLongSchema::PRECURSOR_CHARGE);
      columns.push_back(SpectraLongSchema::PRECURSOR_INTENSITY);
      columns.push_back(SpectraLongSchema::ISOLATION_LOWER);
      columns.push_back(SpectraLongSchema::ISOLATION_UPPER);
    }
  }
  else // SemiWide
  {
    columns = {SpectraSemiWideSchema::SPECTRUM_INDEX, SpectraSemiWideSchema::RT, SpectraSemiWideSchema::MS_LEVEL,
               SpectraSemiWideSchema::NATIVE_ID, SpectraSemiWideSchema::MZ, SpectraSemiWideSchema::INTENSITY};
    if (has_im_data) columns.push_back(SpectraSemiWideSchema::ION_MOBILITY);

    if (config.include_precursor_info)
    {
      columns.push_back(SpectraSemiWideSchema::PRECURSOR_MZ);
      columns.push_back(SpectraSemiWideSchema::PRECURSOR_CHARGE);
      columns.push_back(SpectraSemiWideSchema::PRECURSOR_INTENSITY);
      columns.push_back(SpectraSemiWideSchema::ISOLATION_LOWER);
      columns.push_back(SpectraSemiWideSchema::ISOLATION_UPPER);
    }
  }

  // Filter by requested columns if specified
  if (!config.columns.empty())
  {
    std::vector<std::string> filtered;
    for (const auto& col : columns)
    {
      if (std::find(config.columns.begin(), config.columns.end(), col) != config.columns.end())
      {
        filtered.push_back(col);
      }
    }
    return filtered;
  }

  return columns;
}


namespace
{

// Internal function - not part of public API
std::shared_ptr<arrow::Table> exportChromatogramsToArrow(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config)
{
  const auto& chroms = exp.getChromatograms();
  const auto& cols = config.columns;
  bool inc_rt = shouldIncludeColumn(ChromatogramSchema::RT, cols);
  bool inc_intensity = shouldIncludeColumn(ChromatogramSchema::INTENSITY, cols);
  bool inc_chrom_index = shouldIncludeColumn(ChromatogramSchema::CHROMATOGRAM_INDEX, cols);
  bool inc_native_id = shouldIncludeColumn(ChromatogramSchema::NATIVE_ID, cols);
  bool inc_precursor_mz = shouldIncludeColumn(ChromatogramSchema::PRECURSOR_MZ, cols);
  bool inc_product_mz = shouldIncludeColumn(ChromatogramSchema::PRODUCT_MZ, cols);

  arrow::MemoryPool* pool = arrow::default_memory_pool();
  arrow::Status status;

  if (config.format == ArrowExportFormat::Long)
  {
    // Count total data points
    Size total_points = 0;
    for (const auto& chrom : chroms)
    {
      total_points += chrom.size();
    }

    arrow::DoubleBuilder rt_builder(pool);
    arrow::FloatBuilder intensity_builder(pool);
    arrow::UInt32Builder chrom_index_builder(pool);
    arrow::StringBuilder native_id_builder(pool);
    arrow::DoubleBuilder precursor_mz_builder(pool);
    arrow::DoubleBuilder product_mz_builder(pool);

    if (inc_rt) { status = rt_builder.Reserve(total_points); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_intensity) { status = intensity_builder.Reserve(total_points); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_chrom_index) { status = chrom_index_builder.Reserve(total_points); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_precursor_mz) { status = precursor_mz_builder.Reserve(total_points); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_product_mz) { status = product_mz_builder.Reserve(total_points); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }

    UInt32 chrom_idx = 0;
    for (const auto& chrom : chroms)
    {
      double prec_mz = chrom.getPrecursor().getMZ();
      double prod_mz = chrom.getProduct().getMZ();
      const String& native_id = chrom.getNativeID();

      for (const auto& point : chrom)
      {
        double rt = point.getRT();

        // Apply RT filter
        if (config.min_rt != 0 && rt < config.min_rt) continue;
        if (config.max_rt != 0 && rt > config.max_rt) continue;

        if (inc_rt) rt_builder.UnsafeAppend(rt);
        if (inc_intensity) intensity_builder.UnsafeAppend(static_cast<float>(point.getIntensity()));
        if (inc_chrom_index) chrom_index_builder.UnsafeAppend(chrom_idx);

        if (inc_native_id)
        {
          status = native_id_builder.Append(native_id);
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Append failed: " << status.ToString() << std::endl; return nullptr; }
        }

        if (inc_precursor_mz) precursor_mz_builder.UnsafeAppend(prec_mz);
        if (inc_product_mz) product_mz_builder.UnsafeAppend(prod_mz);
      }
      ++chrom_idx;
    }

    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    std::shared_ptr<arrow::Array> arr;

    if (inc_rt) { status = rt_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::RT, arrow::float64())); arrays.push_back(arr); }
    if (inc_intensity) { status = intensity_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow intensity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::INTENSITY, arrow::float32())); arrays.push_back(arr); }
    if (inc_chrom_index) { status = chrom_index_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow chrom_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::CHROMATOGRAM_INDEX, arrow::uint32())); arrays.push_back(arr); }
    if (inc_native_id) { status = native_id_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow native_id_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::NATIVE_ID, arrow::utf8())); arrays.push_back(arr); }
    if (inc_precursor_mz) { status = precursor_mz_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::PRECURSOR_MZ, arrow::float64())); arrays.push_back(arr); }
    if (inc_product_mz) { status = product_mz_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow product_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSchema::PRODUCT_MZ, arrow::float64())); arrays.push_back(arr); }

    auto schema = arrow::schema(fields);
    auto table = arrow::Table::Make(schema, arrays);

    // Validate table against registry schema (subset — dynamic columns)
    auto validation = ArrowSchemaValidation::validate(table, ChromatogramSchema::schema(), ArrowSchemaValidation::Mode::Subset);
    if (!validation.valid)
    {
      OPENMS_LOG_ERROR << "MSExperimentArrowExport: Chromatogram long format schema validation failed: " << validation.toString() << "\n";
      return nullptr;
    }

    return table;
  }
  else // SemiWide
  {
    arrow::UInt32Builder chrom_index_builder(pool);
    arrow::StringBuilder native_id_builder(pool);
    arrow::DoubleBuilder precursor_mz_builder(pool);
    arrow::DoubleBuilder product_mz_builder(pool);

    auto rt_value_builder = std::make_shared<arrow::DoubleBuilder>(pool);
    arrow::ListBuilder rt_list_builder(pool, rt_value_builder);

    auto intensity_value_builder = std::make_shared<arrow::FloatBuilder>(pool);
    arrow::ListBuilder intensity_list_builder(pool, intensity_value_builder);

    Size num_chroms = chroms.size();
    if (inc_chrom_index) { status = chrom_index_builder.Reserve(num_chroms); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_precursor_mz) { status = precursor_mz_builder.Reserve(num_chroms); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }
    if (inc_product_mz) { status = product_mz_builder.Reserve(num_chroms); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Reserve failed: " << status.ToString() << std::endl; return nullptr; } }

    UInt32 chrom_idx = 0;
    for (const auto& chrom : chroms)
    {
      if (inc_chrom_index) chrom_index_builder.UnsafeAppend(chrom_idx);

      if (inc_native_id)
      {
        status = native_id_builder.Append(chrom.getNativeID());
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Append failed: " << status.ToString() << std::endl; return nullptr; }
      }

      if (inc_precursor_mz) precursor_mz_builder.UnsafeAppend(chrom.getPrecursor().getMZ());
      if (inc_product_mz) product_mz_builder.UnsafeAppend(chrom.getProduct().getMZ());

      if (inc_rt)
      {
        status = rt_list_builder.Append();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

        for (const auto& point : chrom)
        {
          double rt = point.getRT();
          if (config.min_rt != 0 && rt < config.min_rt) continue;
          if (config.max_rt != 0 && rt > config.max_rt) continue;

          status = rt_value_builder->Append(rt);
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      if (inc_intensity)
      {
        status = intensity_list_builder.Append();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

        for (const auto& point : chrom)
        {
          double rt = point.getRT();
          if (config.min_rt != 0 && rt < config.min_rt) continue;
          if (config.max_rt != 0 && rt > config.max_rt) continue;

          status = intensity_value_builder->Append(static_cast<float>(point.getIntensity()));
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }

      ++chrom_idx;
    }

    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    std::shared_ptr<arrow::Array> arr;

    if (inc_chrom_index) { status = chrom_index_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::CHROMATOGRAM_INDEX, arrow::uint32())); arrays.push_back(arr); }
    if (inc_native_id) { status = native_id_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::NATIVE_ID, arrow::utf8())); arrays.push_back(arr); }
    if (inc_rt) { status = rt_list_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::RT, arrow::list(arrow::float64()))); arrays.push_back(arr); }
    if (inc_intensity) { status = intensity_list_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::INTENSITY, arrow::list(arrow::float32()))); arrays.push_back(arr); }
    if (inc_precursor_mz) { status = precursor_mz_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::PRECURSOR_MZ, arrow::float64())); arrays.push_back(arr); }
    if (inc_product_mz) { status = product_mz_builder.Finish(&arr); if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; } fields.push_back(arrow::field(ChromatogramSemiWideSchema::PRODUCT_MZ, arrow::float64())); arrays.push_back(arr); }

    auto schema = arrow::schema(fields);
    auto table = arrow::Table::Make(schema, arrays);

    // Validate table against registry schema (subset — dynamic columns)
    auto validation = ArrowSchemaValidation::validate(table, ChromatogramSemiWideSchema::schema(), ArrowSchemaValidation::Mode::Subset);
    if (!validation.valid)
    {
      OPENMS_LOG_ERROR << "MSExperimentArrowExport: Chromatogram semi-wide format schema validation failed: " << validation.toString() << "\n";
      return nullptr;
    }

    return table;
  }
}

} // anonymous namespace


std::vector<std::string> MSExperimentArrowExport::getChromatogramArrowColumnNames(
  const MSExperiment& /* exp */,
  const ArrowChromatogramExportConfig& config)
{
  std::vector<std::string> columns;

  if (config.format == ArrowExportFormat::Long)
  {
    columns = {ChromatogramSchema::RT, ChromatogramSchema::INTENSITY, ChromatogramSchema::CHROMATOGRAM_INDEX,
               ChromatogramSchema::NATIVE_ID, ChromatogramSchema::PRECURSOR_MZ, ChromatogramSchema::PRODUCT_MZ};
  }
  else // SemiWide
  {
    columns = {ChromatogramSemiWideSchema::CHROMATOGRAM_INDEX, ChromatogramSemiWideSchema::NATIVE_ID,
               ChromatogramSemiWideSchema::RT, ChromatogramSemiWideSchema::INTENSITY,
               ChromatogramSemiWideSchema::PRECURSOR_MZ, ChromatogramSemiWideSchema::PRODUCT_MZ};
  }

  // Filter by requested columns if specified
  if (!config.columns.empty())
  {
    std::vector<std::string> filtered;
    for (const auto& col : columns)
    {
      if (std::find(config.columns.begin(), config.columns.end(), col) != config.columns.end())
      {
        filtered.push_back(col);
      }
    }
    return filtered;
  }

  return columns;
}


bool MSExperimentArrowExport::exportSpectraToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array)
{
  // Build the Arrow table
  auto table = exportSpectraToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to create Arrow table for spectra" << std::endl;
    return false;
  }

  // Convert to RecordBatch (required for C Data Interface export)
  auto batch_result = table->CombineChunksToBatch();
  if (!batch_result.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to combine table chunks: "
                     << batch_result.status().ToString() << std::endl;
    return false;
  }
  const auto& batch = batch_result.ValueOrDie();

  // Export schema
  auto schema_status = arrow::ExportSchema(*batch->schema(), out_schema);
  if (!schema_status.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to export schema: "
                     << schema_status.ToString() << std::endl;
    return false;
  }

  // Export record batch as struct array
  auto array_status = arrow::ExportRecordBatch(*batch, out_array);
  if (!array_status.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to export record batch: "
                     << array_status.ToString() << std::endl;
    // Release the schema on error
    if (out_schema->release)
    {
      out_schema->release(out_schema);
    }
    return false;
  }

  return true;
}


bool MSExperimentArrowExport::exportChromatogramsToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array)
{
  // Build the Arrow table
  auto table = exportChromatogramsToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to create Arrow table for chromatograms" << std::endl;
    return false;
  }

  // Convert to RecordBatch (required for C Data Interface export)
  auto batch_result = table->CombineChunksToBatch();
  if (!batch_result.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to combine table chunks: "
                     << batch_result.status().ToString() << std::endl;
    return false;
  }
  const auto& batch = batch_result.ValueOrDie();

  // Export schema
  auto schema_status = arrow::ExportSchema(*batch->schema(), out_schema);
  if (!schema_status.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to export schema: "
                     << schema_status.ToString() << std::endl;
    return false;
  }

  // Export record batch as struct array
  auto array_status = arrow::ExportRecordBatch(*batch, out_array);
  if (!array_status.ok())
  {
    OPENMS_LOG_ERROR << "MSExperimentArrowExport: Failed to export record batch: "
                     << array_status.ToString() << std::endl;
    // Release the schema on error
    if (out_schema->release)
    {
      out_schema->release(out_schema);
    }
    return false;
  }

  return true;
}


namespace // anonymous namespace for Parquet helpers
{

/// Convert our compression enum to Arrow compression type
arrow::Compression::type toArrowCompression(ParquetWriteConfig::Compression compression)
{
  switch (compression)
  {
    case ParquetWriteConfig::Compression::NONE:
      return arrow::Compression::UNCOMPRESSED;
    case ParquetWriteConfig::Compression::SNAPPY:
      return arrow::Compression::SNAPPY;
    case ParquetWriteConfig::Compression::GZIP:
      return arrow::Compression::GZIP;
    case ParquetWriteConfig::Compression::LZ4:
      return arrow::Compression::LZ4;
    case ParquetWriteConfig::Compression::ZSTD:
      return arrow::Compression::ZSTD;
    default:
      return arrow::Compression::ZSTD;
  }
}


/// Write an Arrow table to a Parquet file
bool writeTableToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const String& filename,
  const ParquetWriteConfig& config)
{
  // Open output file
  auto file_result = arrow::io::FileOutputStream::Open(filename);
  if (!file_result.ok())
  {
    OPENMS_LOG_ERROR << "ParquetExport: Failed to open file for writing: "
                     << filename << " - " << file_result.status().ToString() << std::endl;
    return false;
  }
  const auto& outfile = file_result.ValueOrDie();

  // Configure Parquet writer properties
  auto builder = parquet::WriterProperties::Builder();
  builder.compression(toArrowCompression(config.compression));

  // Set compression level for algorithms that support it
  if (config.compression == ParquetWriteConfig::Compression::ZSTD ||
      config.compression == ParquetWriteConfig::Compression::GZIP)
  {
    builder.compression_level(config.compression_level);
  }

  // Configure data page size
  builder.data_pagesize(config.data_page_size);

  // Enable/disable statistics
  if (config.write_statistics)
  {
    builder.enable_statistics();
  }
  else
  {
    builder.disable_statistics();
  }

  auto writer_properties = builder.build();

  // Configure Arrow writer properties (row group size)
  auto arrow_properties = parquet::ArrowWriterProperties::Builder()
    .store_schema()  // Store Arrow schema for better type fidelity
    ->build();

  // Write the table
  auto status = parquet::arrow::WriteTable(
    *table,
    arrow::default_memory_pool(),
    outfile,
    config.row_group_size,
    writer_properties,
    arrow_properties
  );

  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "ParquetExport: Failed to write Parquet file: "
                     << status.ToString() << std::endl;
    return false;
  }

  // Close the file
  auto close_status = outfile->Close();
  if (!close_status.ok())
  {
    OPENMS_LOG_ERROR << "ParquetExport: Failed to close file: "
                     << close_status.ToString() << std::endl;
    return false;
  }

  return true;
}

} // anonymous namespace


bool MSExperimentArrowExport::exportSpectraToParquet(
  const MSExperiment& exp,
  const String& filename,
  const ArrowSpectraExportConfig& config,
  const ParquetWriteConfig& parquet_config)
{
  // Build Arrow table
  auto table = exportSpectraToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ParquetExport: Failed to create Arrow table for spectra" << std::endl;
    return false;
  }

  // Write to Parquet file
  return writeTableToParquet(table, filename, parquet_config);
}


bool MSExperimentArrowExport::exportChromatogramsToParquet(
  const MSExperiment& exp,
  const String& filename,
  const ArrowChromatogramExportConfig& config,
  const ParquetWriteConfig& parquet_config)
{
  // Build Arrow table
  auto table = exportChromatogramsToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ParquetExport: Failed to create Arrow table for chromatograms" << std::endl;
    return false;
  }

  // Write to Parquet file
  return writeTableToParquet(table, filename, parquet_config);
}


} // namespace OpenMS
