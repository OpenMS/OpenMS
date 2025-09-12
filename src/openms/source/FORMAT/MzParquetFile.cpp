// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzParquetFile.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/compute/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/exception.h>

#include <unordered_map>
#include <vector>

using namespace std;

namespace OpenMS
{

  MzParquetFile::MzParquetFile() :
    ProgressLogger()
  {
  }

  MzParquetFile::~MzParquetFile() = default;

  void MzParquetFile::load(const String& filename, PeakMap& map)
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    startProgress(0, 0, "Loading mzParquet file");

    try
    {
      readParquetFile_(filename, map);
    }
    catch (const std::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        String("Error reading mzParquet file: ") + e.what(), filename);
    }

    endProgress();
  }

  void MzParquetFile::store(const String& /* filename */, const PeakMap& /* map */) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void MzParquetFile::readParquetFile_(const String& filename, PeakMap& map)
  {
    // Clear existing data
    map.clear(true);

    // Open parquet file
    std::shared_ptr<arrow::io::ReadableFile> infile;
    auto file_result = arrow::io::ReadableFile::Open(filename.c_str());
    if (!file_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to open parquet file: " + String(file_result.status().ToString()), filename);
    }
    infile = file_result.ValueOrDie();

    // Create parquet reader using simpler approach
    std::unique_ptr<parquet::arrow::FileReader> reader;
    auto result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!result.ok()) {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to create parquet reader: " + String(result.status().ToString()), filename);
    }
    reader = std::move(result.ValueOrDie());

    // Get file metadata to understand the structure
    auto parquet_metadata = reader->parquet_reader()->metadata();
    int num_row_groups = parquet_metadata->num_row_groups();

    // Get unique scan IDs efficiently using Arrow's compute functions
    std::shared_ptr<arrow::Table> scan_column_table;
    auto status = reader->ReadTable({0}, &scan_column_table); // Read only scan column (column 0)
    if (!status.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to read scan column: " + String(status.ToString()), filename);
    }
    
    auto scan_column = scan_column_table->column(0);
    auto unique_result = arrow::compute::CallFunction("unique", {scan_column});
    if (!unique_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to get unique scans: " + String(unique_result.status().ToString()), filename);
    }
    
    auto unique_array = unique_result.ValueOrDie().array_as<arrow::UInt32Array>();

    // Pre-allocate spectra map with exact scan IDs
    std::unordered_map<UInt32, MSSpectrum> spectra_map;
    spectra_map.reserve(unique_array->length()); // Reserve space for better performance
    for (int64_t i = 0; i < unique_array->length(); ++i)
    {
      if (!unique_array->IsNull(i))
      {
        spectra_map[unique_array->Value(i)] = MSSpectrum();
      }
    }

    // Second pass: process row groups in parallel if beneficial
    #pragma omp parallel for if(num_row_groups > 2)
    for (int rg = 0; rg < num_row_groups; ++rg)
    {
      // Read entire row group
      std::shared_ptr<arrow::Table> table;
      auto status = reader->ReadRowGroup(rg, &table);
      if (!status.ok())
      {
        #pragma omp critical
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Failed to read row group: " + String(status.ToString()), filename);
        }
      }

      // Process this row group's data
      processRowGroup_(table, spectra_map);
    }

    // Convert map to vector and add to MSExperiment
    map.reserveSpaceSpectra(spectra_map.size());
    for (auto& [scan_id, spectrum] : spectra_map)
    {
      // Sort peaks by m/z within each spectrum
      spectrum.sortByPosition();
      map.addSpectrum(std::move(spectrum));
    }

    // Sort spectra by RT
    map.sortSpectra(true);
  }

  void MzParquetFile::processRowGroup_(const std::shared_ptr<arrow::Table>& table, 
                                       std::unordered_map<UInt32, MSSpectrum>& spectra_map)
  {
    // Get all columns
    auto scan_column = table->GetColumnByName("scan");
    auto level_column = table->GetColumnByName("level");
    auto rt_column = table->GetColumnByName("rt");
    auto mz_column = table->GetColumnByName("mz");
    auto intensity_column = table->GetColumnByName("intensity");

    // Optional columns
    auto collision_energy_column = table->GetColumnByName("collision_energy");
    auto ion_mobility_column = table->GetColumnByName("ion_mobility");
    auto isolation_lower_column = table->GetColumnByName("isolation_lower");
    auto isolation_upper_column = table->GetColumnByName("isolation_upper");
    auto precursor_scan_column = table->GetColumnByName("precursor_scan");
    auto precursor_mz_column = table->GetColumnByName("precursor_mz");
    auto precursor_charge_column = table->GetColumnByName("precursor_charge");

    // Process all chunks in each column
    for (int chunk_idx = 0; chunk_idx < scan_column->num_chunks(); ++chunk_idx)
    {
      // Get arrays for this chunk
      auto scan_array = std::static_pointer_cast<arrow::UInt32Array>(scan_column->chunk(chunk_idx));
      auto level_array = std::static_pointer_cast<arrow::UInt32Array>(level_column->chunk(chunk_idx));
      auto rt_array = std::static_pointer_cast<arrow::FloatArray>(rt_column->chunk(chunk_idx));
      auto mz_array = std::static_pointer_cast<arrow::FloatArray>(mz_column->chunk(chunk_idx));
      auto intensity_array = std::static_pointer_cast<arrow::UInt32Array>(intensity_column->chunk(chunk_idx));

      // Optional arrays for this chunk
      auto collision_energy_array = collision_energy_column ? 
        std::static_pointer_cast<arrow::FloatArray>(collision_energy_column->chunk(chunk_idx)) : nullptr;
      auto ion_mobility_array = ion_mobility_column ? 
        std::static_pointer_cast<arrow::FloatArray>(ion_mobility_column->chunk(chunk_idx)) : nullptr;
      auto isolation_lower_array = isolation_lower_column ? 
        std::static_pointer_cast<arrow::FloatArray>(isolation_lower_column->chunk(chunk_idx)) : nullptr;
      auto isolation_upper_array = isolation_upper_column ? 
        std::static_pointer_cast<arrow::FloatArray>(isolation_upper_column->chunk(chunk_idx)) : nullptr;
      auto precursor_scan_array = precursor_scan_column ? 
        std::static_pointer_cast<arrow::UInt32Array>(precursor_scan_column->chunk(chunk_idx)) : nullptr;
      auto precursor_mz_array = precursor_mz_column ? 
        std::static_pointer_cast<arrow::FloatArray>(precursor_mz_column->chunk(chunk_idx)) : nullptr;
      auto precursor_charge_array = precursor_charge_column ? 
        std::static_pointer_cast<arrow::UInt32Array>(precursor_charge_column->chunk(chunk_idx)) : nullptr;

      // Process each row in the chunk
      for (int64_t i = 0; i < scan_array->length(); ++i)
      {
        UInt32 scan_id = scan_array->Value(i);
        
        // Get or create spectrum (thread-safe access needed)
        #pragma omp critical(spectra_access)
        {
          auto& spectrum = spectra_map[scan_id];
          
          // Set spectrum metadata on first peak of this scan
          if (spectrum.empty())
          {
            spectrum.setRT(rt_array->Value(i));
            spectrum.setMSLevel(level_array->Value(i));
            spectrum.setNativeID(String("scan=") + scan_id);

            // Set ion mobility if available
            if (ion_mobility_array && !ion_mobility_array->IsNull(i))
            {
              spectrum.setDriftTime(ion_mobility_array->Value(i));
            }

            // Set precursor information for MS2+ spectra
            if (level_array->Value(i) > 1)
            {
              Precursor precursor;
              
              if (precursor_mz_array && !precursor_mz_array->IsNull(i))
              {
                precursor.setMZ(precursor_mz_array->Value(i));
              }
              
              if (precursor_charge_array && !precursor_charge_array->IsNull(i))
              {
                precursor.setCharge(static_cast<Int>(precursor_charge_array->Value(i)));
              }

              // Set collision energy if available
              if (collision_energy_array && !collision_energy_array->IsNull(i))
              {
                precursor.setActivationEnergy(collision_energy_array->Value(i));
              }

              // Set isolation window
              if (isolation_lower_array && isolation_upper_array && 
                  !isolation_lower_array->IsNull(i) && !isolation_upper_array->IsNull(i))
              {
                float precursor_mz_val = precursor_mz_array ? precursor_mz_array->Value(i) : 0.0f;
                precursor.setIsolationWindowLowerOffset(precursor_mz_val - isolation_lower_array->Value(i));
                precursor.setIsolationWindowUpperOffset(isolation_upper_array->Value(i) - precursor_mz_val);
              }

              spectrum.getPrecursors().push_back(precursor);
            }
          }

          // Add peak to spectrum
          Peak1D peak;
          peak.setMZ(mz_array->Value(i));
          peak.setIntensity(static_cast<float>(intensity_array->Value(i)));
          spectrum.push_back(peak);
        }
      }
    }
  }

} // namespace OpenMS

#endif // WITH_PARQUET
