// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowExport.h>

#ifdef WITH_PARQUET

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

/// Check if spectrum passes filter criteria
bool passesFilter(const MSSpectrum& spec,
                  const ArrowSpectraExportConfig& config,
                  const std::unordered_set<UInt>& ms_levels_set)
{
  // Check MS level filter
  if (!ms_levels_set.empty() && ms_levels_set.find(spec.getMSLevel()) == ms_levels_set.end())
  {
    return false;
  }

  // Check RT range filter
  double rt = spec.getRT();
  if (config.min_rt != 0 || config.max_rt != 0)
  {
    if (config.min_rt != 0 && rt < config.min_rt) return false;
    if (config.max_rt != 0 && rt > config.max_rt) return false;
  }

  return true;
}

/// Check if a peak passes m/z filter
inline bool peakPassesMzFilter(double mz, double min_mz, double max_mz)
{
  if (min_mz == 0 && max_mz == 0) return true;
  if (min_mz != 0 && mz < min_mz) return false;
  if (max_mz != 0 && mz > max_mz) return false;
  return true;
}

/// Check if any spectrum in the experiment has ion mobility data
bool experimentHasIMData(const MSExperiment& exp,
                         const std::unordered_set<UInt>& ms_levels_set,
                         const ArrowSpectraExportConfig& config)
{
  for (const auto& spec : exp)
  {
    if (!passesFilter(spec, config, ms_levels_set)) continue;
    if (spec.containsIMData()) return true;
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

  for (const auto& spec : exp)
  {
    if (!passesFilter(spec, config, ms_levels_set)) continue;

    if (filter_mz)
    {
      // Need to count peaks in m/z range
      for (const auto& peak : spec)
      {
        if (peakPassesMzFilter(peak.getMZ(), config.min_mz, config.max_mz))
        {
          ++total;
        }
      }
    }
    else
    {
      total += spec.size();
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
  for (const auto& spec : exp)
  {
    if (passesFilter(spec, config, ms_levels_set)) ++count;
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
  bool inc_mz = shouldIncludeColumn("mz", cols);
  bool inc_intensity = shouldIncludeColumn("intensity", cols);
  bool inc_rt = shouldIncludeColumn("rt", cols);
  bool inc_spectrum_index = shouldIncludeColumn("spectrum_index", cols);
  bool inc_ms_level = shouldIncludeColumn("ms_level", cols);
  bool inc_native_id = shouldIncludeColumn("native_id", cols);
  bool inc_ion_mobility = config.include_ion_mobility && shouldIncludeColumn("ion_mobility", cols);
  bool inc_precursor_mz = config.include_precursor_info && shouldIncludeColumn("precursor_mz", cols);
  bool inc_precursor_charge = config.include_precursor_info && shouldIncludeColumn("precursor_charge", cols);
  bool inc_precursor_intensity = config.include_precursor_info && shouldIncludeColumn("precursor_intensity", cols);
  bool inc_isolation_lower = config.include_precursor_info && shouldIncludeColumn("isolation_lower", cols);
  bool inc_isolation_upper = config.include_precursor_info && shouldIncludeColumn("isolation_upper", cols);

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

  // m/z filter flags
  bool filter_mz = (config.min_mz != 0 || config.max_mz != 0);

  // Iterate through spectra
  UInt32 spectrum_idx = 0;
  for (const auto& spec : exp)
  {
    if (!passesFilter(spec, config, ms_levels_set))
    {
      ++spectrum_idx;
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

    // Iterate through peaks
    Size peak_idx = 0;
    for (const auto& peak : spec)
    {
      double mz = peak.getMZ();

      // Apply m/z filter
      if (filter_mz && !peakPassesMzFilter(mz, config.min_mz, config.max_mz))
      {
        ++peak_idx;
        continue;
      }

      // Append peak data using UnsafeAppend (we already reserved capacity)
      if (inc_mz) mz_builder.UnsafeAppend(mz);
      if (inc_intensity) intensity_builder.UnsafeAppend(static_cast<float>(peak.getIntensity()));
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

      ++peak_idx;
    }

    ++spectrum_idx;
  }

  // Build schema and arrays
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  std::shared_ptr<arrow::Array> arr;

  if (inc_mz)
  {
    status = mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("mz", arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_intensity)
  {
    status = intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow intensity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("intensity", arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_rt)
  {
    status = rt_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("rt", arrow::float32()));
    arrays.push_back(arr);
  }

  if (has_im_data)
  {
    status = ion_mobility_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("ion_mobility", arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_spectrum_index)
  {
    status = spectrum_index_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow spectrum_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("spectrum_index", arrow::uint32()));
    arrays.push_back(arr);
  }

  if (inc_ms_level)
  {
    status = ms_level_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow ms_level_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("ms_level", arrow::uint8()));
    arrays.push_back(arr);
  }

  if (inc_native_id)
  {
    status = native_id_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow native_id_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("native_id", arrow::utf8()));
    arrays.push_back(arr);
  }

  if (inc_precursor_mz)
  {
    status = precursor_mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_mz", arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_precursor_charge)
  {
    status = precursor_charge_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_charge", arrow::int16()));
    arrays.push_back(arr);
  }

  if (inc_precursor_intensity)
  {
    status = precursor_intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow precursor_intensity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_intensity", arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_isolation_lower)
  {
    status = isolation_lower_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_lower_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("isolation_lower", arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_isolation_upper)
  {
    status = isolation_upper_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow isolation_upper_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("isolation_upper", arrow::float64()));
    arrays.push_back(arr);
  }

  auto schema = arrow::schema(fields);
  return arrow::Table::Make(schema, arrays);
}


/// Build semi-wide format Arrow table (one row per spectrum)
std::shared_ptr<arrow::Table> buildSemiWideFormatTable(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  const std::unordered_set<UInt>& ms_levels_set)
{
  // Determine which columns to include
  const auto& cols = config.columns;
  bool inc_mz = shouldIncludeColumn("mz", cols);
  bool inc_intensity = shouldIncludeColumn("intensity", cols);
  bool inc_rt = shouldIncludeColumn("rt", cols);
  bool inc_spectrum_index = shouldIncludeColumn("spectrum_index", cols);
  bool inc_ms_level = shouldIncludeColumn("ms_level", cols);
  bool inc_native_id = shouldIncludeColumn("native_id", cols);
  bool inc_ion_mobility = config.include_ion_mobility && shouldIncludeColumn("ion_mobility", cols);
  bool inc_precursor_mz = config.include_precursor_info && shouldIncludeColumn("precursor_mz", cols);
  bool inc_precursor_charge = config.include_precursor_info && shouldIncludeColumn("precursor_charge", cols);
  bool inc_precursor_intensity = config.include_precursor_info && shouldIncludeColumn("precursor_intensity", cols);
  bool inc_isolation_lower = config.include_precursor_info && shouldIncludeColumn("isolation_lower", cols);
  bool inc_isolation_upper = config.include_precursor_info && shouldIncludeColumn("isolation_upper", cols);

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

  // m/z filter flags
  bool filter_mz = (config.min_mz != 0 || config.max_mz != 0);

  // Iterate through spectra
  UInt32 spectrum_idx = 0;
  for (const auto& spec : exp)
  {
    if (!passesFilter(spec, config, ms_levels_set))
    {
      ++spectrum_idx;
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

    // Append list columns (mz, intensity, ion_mobility arrays)
    if (inc_mz)
    {
      status = mz_list_builder.Append();
      if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

      Size peak_idx = 0;
      for (const auto& peak : spec)
      {
        double mz = peak.getMZ();
        if (!filter_mz || peakPassesMzFilter(mz, config.min_mz, config.max_mz))
        {
          status = mz_value_builder->Append(mz);
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
        }
        ++peak_idx;
      }
    }

    if (inc_intensity)
    {
      status = intensity_list_builder.Append();
      if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

      for (const auto& peak : spec)
      {
        if (!filter_mz || peakPassesMzFilter(peak.getMZ(), config.min_mz, config.max_mz))
        {
          status = intensity_value_builder->Append(static_cast<float>(peak.getIntensity()));
          if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
        }
      }
    }

    if (has_im_data)
    {
      if (im_data_ptr)
      {
        status = im_list_builder.Append();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow list Append failed: " << status.ToString() << std::endl; return nullptr; }

        Size peak_idx = 0;
        for (const auto& peak : spec)
        {
          if (!filter_mz || peakPassesMzFilter(peak.getMZ(), config.min_mz, config.max_mz))
          {
            status = im_value_builder->Append(im_data_ptr[peak_idx]);
            if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow value Append failed: " << status.ToString() << std::endl; return nullptr; }
          }
          ++peak_idx;
        }
      }
      else
      {
        status = im_list_builder.AppendNull();
        if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow AppendNull failed: " << status.ToString() << std::endl; return nullptr; }
      }
    }

    ++spectrum_idx;
  }

  // Build schema and arrays
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;

  std::shared_ptr<arrow::Array> arr;

  if (inc_spectrum_index)
  {
    status = spectrum_index_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("spectrum_index", arrow::uint32()));
    arrays.push_back(arr);
  }

  if (inc_rt)
  {
    status = rt_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("rt", arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_ms_level)
  {
    status = ms_level_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("ms_level", arrow::uint8()));
    arrays.push_back(arr);
  }

  if (inc_native_id)
  {
    status = native_id_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("native_id", arrow::utf8()));
    arrays.push_back(arr);
  }

  if (inc_mz)
  {
    status = mz_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("mz", arrow::list(arrow::float64())));
    arrays.push_back(arr);
  }

  if (inc_intensity)
  {
    status = intensity_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("intensity", arrow::list(arrow::float32())));
    arrays.push_back(arr);
  }

  if (has_im_data)
  {
    status = im_list_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("ion_mobility", arrow::list(arrow::float32())));
    arrays.push_back(arr);
  }

  if (inc_precursor_mz)
  {
    status = precursor_mz_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_mz", arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_precursor_charge)
  {
    status = precursor_charge_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_charge", arrow::int16()));
    arrays.push_back(arr);
  }

  if (inc_precursor_intensity)
  {
    status = precursor_intensity_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("precursor_intensity", arrow::float32()));
    arrays.push_back(arr);
  }

  if (inc_isolation_lower)
  {
    status = isolation_lower_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("isolation_lower", arrow::float64()));
    arrays.push_back(arr);
  }

  if (inc_isolation_upper)
  {
    status = isolation_upper_builder.Finish(&arr);
    if (!status.ok()) { OPENMS_LOG_ERROR << "Arrow Finish failed: " << status.ToString() << std::endl; return nullptr; }
    fields.push_back(arrow::field("isolation_upper", arrow::float64()));
    arrays.push_back(arr);
  }

  auto schema = arrow::schema(fields);
  return arrow::Table::Make(schema, arrays);
}

} // anonymous namespace


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


std::vector<std::string> getSpectraArrowColumns(
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
    columns = {"mz", "intensity", "rt"};
    if (has_im_data) columns.push_back("ion_mobility");
    columns.push_back("spectrum_index");
    columns.push_back("ms_level");
    columns.push_back("native_id");

    if (config.include_precursor_info)
    {
      columns.push_back("precursor_mz");
      columns.push_back("precursor_charge");
      columns.push_back("precursor_intensity");
      columns.push_back("isolation_lower");
      columns.push_back("isolation_upper");
    }
  }
  else // SemiWide
  {
    columns = {"spectrum_index", "rt", "ms_level", "native_id", "mz", "intensity"};
    if (has_im_data) columns.push_back("ion_mobility");

    if (config.include_precursor_info)
    {
      columns.push_back("precursor_mz");
      columns.push_back("precursor_charge");
      columns.push_back("precursor_intensity");
      columns.push_back("isolation_lower");
      columns.push_back("isolation_upper");
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


std::shared_ptr<arrow::Table> exportChromatogramsToArrow(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config)
{
  const auto& chroms = exp.getChromatograms();

  if (chroms.empty())
  {
    // Return empty table with schema
    std::vector<std::shared_ptr<arrow::Field>> fields;
    if (config.format == ArrowExportFormat::Long)
    {
      fields = {
        arrow::field("rt", arrow::float64()),
        arrow::field("intensity", arrow::float32()),
        arrow::field("chromatogram_index", arrow::uint32()),
        arrow::field("native_id", arrow::utf8()),
        arrow::field("precursor_mz", arrow::float64()),
        arrow::field("product_mz", arrow::float64())
      };
    }
    else
    {
      fields = {
        arrow::field("chromatogram_index", arrow::uint32()),
        arrow::field("native_id", arrow::utf8()),
        arrow::field("rt", arrow::list(arrow::float64())),
        arrow::field("intensity", arrow::list(arrow::float32())),
        arrow::field("precursor_mz", arrow::float64()),
        arrow::field("product_mz", arrow::float64())
      };
    }
    auto schema = arrow::schema(fields);
    std::vector<std::shared_ptr<arrow::Array>> empty_arrays;
    for (size_t i = 0; i < fields.size(); ++i)
    {
      std::shared_ptr<arrow::Array> empty_arr;
      arrow::Status status;

      if (fields[i]->type()->id() == arrow::Type::FLOAT)
      {
        arrow::FloatBuilder builder;
        status = builder.Finish(&empty_arr);
      }
      else if (fields[i]->type()->id() == arrow::Type::DOUBLE)
      {
        arrow::DoubleBuilder builder;
        status = builder.Finish(&empty_arr);
      }
      else if (fields[i]->type()->id() == arrow::Type::UINT32)
      {
        arrow::UInt32Builder builder;
        status = builder.Finish(&empty_arr);
      }
      else if (fields[i]->type()->id() == arrow::Type::STRING)
      {
        arrow::StringBuilder builder;
        status = builder.Finish(&empty_arr);
      }
      else if (fields[i]->type()->id() == arrow::Type::LIST)
      {
        auto inner_type = std::static_pointer_cast<arrow::ListType>(fields[i]->type())->value_type();
        if (inner_type->id() == arrow::Type::DOUBLE)
        {
          auto val_builder = std::make_shared<arrow::DoubleBuilder>();
          arrow::ListBuilder list_builder(arrow::default_memory_pool(), val_builder);
          status = list_builder.Finish(&empty_arr);
        }
        else
        {
          auto val_builder = std::make_shared<arrow::FloatBuilder>();
          arrow::ListBuilder list_builder(arrow::default_memory_pool(), val_builder);
          status = list_builder.Finish(&empty_arr);
        }
      }

      if (!status.ok())
      {
        OPENMS_LOG_ERROR << "Arrow builder Finish failed: " << status.ToString() << std::endl;
        return nullptr;
      }
      empty_arrays.push_back(empty_arr);
    }
    return arrow::Table::Make(schema, empty_arrays);
  }

  const auto& cols = config.columns;
  bool inc_rt = shouldIncludeColumn("rt", cols);
  bool inc_intensity = shouldIncludeColumn("intensity", cols);
  bool inc_chrom_index = shouldIncludeColumn("chromatogram_index", cols);
  bool inc_native_id = shouldIncludeColumn("native_id", cols);
  bool inc_precursor_mz = shouldIncludeColumn("precursor_mz", cols);
  bool inc_product_mz = shouldIncludeColumn("product_mz", cols);

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

    if (inc_rt) { status = rt_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("rt", arrow::float64())); arrays.push_back(arr); }
    if (inc_intensity) { status = intensity_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("intensity", arrow::float32())); arrays.push_back(arr); }
    if (inc_chrom_index) { status = chrom_index_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("chromatogram_index", arrow::uint32())); arrays.push_back(arr); }
    if (inc_native_id) { status = native_id_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("native_id", arrow::utf8())); arrays.push_back(arr); }
    if (inc_precursor_mz) { status = precursor_mz_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("precursor_mz", arrow::float64())); arrays.push_back(arr); }
    if (inc_product_mz) { status = product_mz_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("product_mz", arrow::float64())); arrays.push_back(arr); }

    auto schema = arrow::schema(fields);
    return arrow::Table::Make(schema, arrays);
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
    if (inc_chrom_index) { status = chrom_index_builder.Reserve(num_chroms); if (!status.ok()) return nullptr; }
    if (inc_precursor_mz) { status = precursor_mz_builder.Reserve(num_chroms); if (!status.ok()) return nullptr; }
    if (inc_product_mz) { status = product_mz_builder.Reserve(num_chroms); if (!status.ok()) return nullptr; }

    UInt32 chrom_idx = 0;
    for (const auto& chrom : chroms)
    {
      if (inc_chrom_index) chrom_index_builder.UnsafeAppend(chrom_idx);

      if (inc_native_id)
      {
        status = native_id_builder.Append(chrom.getNativeID());
        if (!status.ok()) return nullptr;
      }

      if (inc_precursor_mz) precursor_mz_builder.UnsafeAppend(chrom.getPrecursor().getMZ());
      if (inc_product_mz) product_mz_builder.UnsafeAppend(chrom.getProduct().getMZ());

      if (inc_rt)
      {
        status = rt_list_builder.Append();
        if (!status.ok()) return nullptr;

        for (const auto& point : chrom)
        {
          double rt = point.getRT();
          if (config.min_rt != 0 && rt < config.min_rt) continue;
          if (config.max_rt != 0 && rt > config.max_rt) continue;

          status = rt_value_builder->Append(rt);
          if (!status.ok()) return nullptr;
        }
      }

      if (inc_intensity)
      {
        status = intensity_list_builder.Append();
        if (!status.ok()) return nullptr;

        for (const auto& point : chrom)
        {
          double rt = point.getRT();
          if (config.min_rt != 0 && rt < config.min_rt) continue;
          if (config.max_rt != 0 && rt > config.max_rt) continue;

          status = intensity_value_builder->Append(static_cast<float>(point.getIntensity()));
          if (!status.ok()) return nullptr;
        }
      }

      ++chrom_idx;
    }

    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    std::shared_ptr<arrow::Array> arr;

    if (inc_chrom_index) { status = chrom_index_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("chromatogram_index", arrow::uint32())); arrays.push_back(arr); }
    if (inc_native_id) { status = native_id_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("native_id", arrow::utf8())); arrays.push_back(arr); }
    if (inc_rt) { status = rt_list_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("rt", arrow::list(arrow::float64()))); arrays.push_back(arr); }
    if (inc_intensity) { status = intensity_list_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("intensity", arrow::list(arrow::float32()))); arrays.push_back(arr); }
    if (inc_precursor_mz) { status = precursor_mz_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("precursor_mz", arrow::float64())); arrays.push_back(arr); }
    if (inc_product_mz) { status = product_mz_builder.Finish(&arr); if (!status.ok()) return nullptr; fields.push_back(arrow::field("product_mz", arrow::float64())); arrays.push_back(arr); }

    auto schema = arrow::schema(fields);
    return arrow::Table::Make(schema, arrays);
  }
}


std::vector<std::string> getChromatogramArrowColumns(
  const MSExperiment& /* exp */,
  const ArrowChromatogramExportConfig& config)
{
  std::vector<std::string> columns;

  if (config.format == ArrowExportFormat::Long)
  {
    columns = {"rt", "intensity", "chromatogram_index", "native_id", "precursor_mz", "product_mz"};
  }
  else // SemiWide
  {
    columns = {"chromatogram_index", "native_id", "rt", "intensity", "precursor_mz", "product_mz"};
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


bool exportSpectraToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array)
{
  // Build the Arrow table
  auto table = exportSpectraToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to create Arrow table for spectra" << std::endl;
    return false;
  }

  // Convert to RecordBatch (required for C Data Interface export)
  auto batch_result = table->CombineChunksToBatch();
  if (!batch_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to combine table chunks: "
                     << batch_result.status().ToString() << std::endl;
    return false;
  }
  auto batch = batch_result.ValueOrDie();

  // Export schema
  auto schema_status = arrow::ExportSchema(*batch->schema(), out_schema);
  if (!schema_status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to export schema: "
                     << schema_status.ToString() << std::endl;
    return false;
  }

  // Export record batch as struct array
  auto array_status = arrow::ExportRecordBatch(*batch, out_array);
  if (!array_status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to export record batch: "
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


bool exportChromatogramsToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array)
{
  // Build the Arrow table
  auto table = exportChromatogramsToArrow(exp, config);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to create Arrow table for chromatograms" << std::endl;
    return false;
  }

  // Convert to RecordBatch (required for C Data Interface export)
  auto batch_result = table->CombineChunksToBatch();
  if (!batch_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to combine table chunks: "
                     << batch_result.status().ToString() << std::endl;
    return false;
  }
  auto batch = batch_result.ValueOrDie();

  // Export schema
  auto schema_status = arrow::ExportSchema(*batch->schema(), out_schema);
  if (!schema_status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to export schema: "
                     << schema_status.ToString() << std::endl;
    return false;
  }

  // Export record batch as struct array
  auto array_status = arrow::ExportRecordBatch(*batch, out_array);
  if (!array_status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowExport: Failed to export record batch: "
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
  auto outfile = file_result.ValueOrDie();

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


bool exportSpectraToParquet(
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


bool exportChromatogramsToParquet(
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

#endif // WITH_PARQUET
