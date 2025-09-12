// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#ifdef WITH_PARQUET

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <unordered_map>
#include <memory>

// Forward declarations for Arrow types
namespace arrow
{
  class Table;
  class RecordBatch;
}

namespace OpenMS
{
  /**
    @brief File adapter for mzParquet files

    This implementation supports reading mzParquet files into MSExperiment objects.
    The mzParquet format is a columnar storage format for mass spectrometry data
    based on Apache Parquet, designed for efficient storage and processing of large
    MS datasets.

    The expected schema contains the following columns:
    - scan: scan number (required, uint32)
    - level: MS level (required, uint32) 
    - rt: retention time in seconds (required, float)
    - mz: mass-to-charge ratio (required, float)
    - intensity: peak intensity (required, uint32)
    - collision_energy: collision energy (optional, float)
    - ion_mobility: ion mobility value (optional, float)
    - isolation_lower: lower isolation window boundary (optional, float)
    - isolation_upper: upper isolation window boundary (optional, float)
    - precursor_scan: precursor scan number (optional, uint32)
    - precursor_mz: precursor m/z (optional, float)
    - precursor_charge: precursor charge (optional, uint32)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzParquetFile :
    public ProgressLogger
  {
public:
    /// Default constructor
    MzParquetFile();
    
    /// Destructor
    ~MzParquetFile() override;

    /**
      @brief Loads a map from an mzParquet file

      @param filename The filename with the data
      @param map MSExperiment to store the loaded data

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, PeakMap& map);

    /**
      @brief Stores a map in an mzParquet file

      @param filename The filename to write to
      @param map MSExperiment containing the data to store

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
      @exception Exception::NotImplemented is thrown as storing is not yet implemented
    */
    void store(const String& filename, const PeakMap& map) const;

  private:
    /// Copy constructor (not implemented)
    MzParquetFile(const MzParquetFile& rhs);
    
    /// Assignment operator (not implemented)
    MzParquetFile& operator=(const MzParquetFile& rhs);

    /**
      @brief Helper function to read parquet file and extract data

      @param filename The parquet file to read
      @param map MSExperiment to populate with data
    */
    void readParquetFile_(const String& filename, PeakMap& map);

    /**
      @brief Process a single row group from the parquet file

      @param table Arrow table containing the row group data
      @param spectra_map Map of scan IDs to MSSpectrum objects to populate
    */
    void processRowGroup_(const std::shared_ptr<arrow::Table>& table, 
                          std::unordered_map<UInt32, MSSpectrum>& spectra_map);
  };

} // namespace OpenMS

#endif // WITH_PARQUET
