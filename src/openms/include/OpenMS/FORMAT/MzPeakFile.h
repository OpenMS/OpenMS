// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{
  /**
    @brief File adapter for mzPeak files

    The mzPeak format stores mass spectrometry data in a ZIP archive containing
    Parquet files. The format supports:
    - spectra_metadata.mzpeak: Spectrum level metadata
    - spectra_data.mzpeak: Spectrum signal data (profile or centroid)
    - chromatograms_metadata.mzpeak: Chromatogram metadata (optional)
    - chromatograms_data.mzpeak: Chromatogram signal data (optional)

    For format specification, see:
    https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00435

    @note This file format requires Parquet support (WITH_PARQUET=ON during build)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzPeakFile :
    public ProgressLogger
  {
public:
    /// Type alias for map type
    typedef MSExperiment MapType;

    /**
      @brief Configuration class for MzPeakFile

      Contains configuration options for mzPeak file reading/writing
    */
    struct OPENMS_DLLAPI MzPeakConfig 
    {
      bool use_f32_mz{false};           ///< Use float32 for m/z values (default: float64)
      bool use_f32_intensity{true};     ///< Use float32 for intensity values
      bool use_f32_ion_mobility{false}; ///< Use float32 for ion mobility values
      bool write_chromatograms{true};   ///< Include chromatograms in output
      bool use_chunked_layout{false};   ///< Use chunked layout for data arrays (enables Numpress)
      Size buffer_size{5000};           ///< Number of spectra to buffer between writes
    };

    /// @name Constructors and Destructor
    //@{
    /// Default constructor
    MzPeakFile();

    /// Default destructor
    ~MzPeakFile() override;
    //@}

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// Set options for loading/storing
    void setOptions(const PeakFileOptions& options);

    /**
      @brief Loads a map from an mzPeak file.

      @param filename The filename with the data
      @param map Is an MSExperiment to store the data

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, MapType& map);

    /**
      @brief Stores a map in an mzPeak file.

      @param filename The filename to store to
      @param map The MSExperiment to store

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const MapType& map) const;

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer.

      The result will not be stored directly but is available through the
      events triggered by the parser and caught by the provided IMSDataConsumer
      object.

      @param filename_in Filename of input mzPeak file to transform
      @param consumer Consumer class to operate on the input filename
      @param skip_full_count Whether to skip computing the correct number of spectra/chromatograms
      @param skip_first_pass Skip first file parsing pass (meta-data only)
    */
    void transform(const String& filename_in, 
                   Interfaces::IMSDataConsumer* consumer, 
                   bool skip_full_count = false, 
                   bool skip_first_pass = false);

    /**
      @brief Set the configuration for reading/writing

      @param config The configuration to use
    */
    void setConfig(const MzPeakConfig& config);

    /**
      @brief Get the current configuration

      @return The current configuration
    */
    const MzPeakConfig& getConfig() const;

protected:
    /// Configuration for mzPeak operations
    MzPeakConfig config_;
    
    /// Options for peak file handling
    PeakFileOptions options_;

    /**
      @brief Read spectrum metadata from Parquet file within ZIP archive

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment to populate with spectrum metadata
    */
    void readSpectrumMetadata_(const String& archive_path, MapType& map);

    /**
      @brief Read spectrum signal data from Parquet file within ZIP archive

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment to populate with spectrum data
    */
    void readSpectrumData_(const String& archive_path, MapType& map);

    /**
      @brief Read chromatogram metadata and data from Parquet files

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment to populate with chromatogram data
    */
    void readChromatograms_(const String& archive_path, MapType& map);

    /**
      @brief Write spectrum metadata to Parquet file within ZIP archive

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment containing spectrum metadata to write
    */
    void writeSpectrumMetadata_(const String& archive_path, const MapType& map) const;

    /**
      @brief Write spectrum signal data to Parquet file within ZIP archive

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment containing spectrum data to write
    */
    void writeSpectrumData_(const String& archive_path, const MapType& map) const;

    /**
      @brief Write chromatogram data to Parquet files

      @param archive_path Path to the mzPeak archive
      @param map MSExperiment containing chromatogram data to write
    */
    void writeChromatograms_(const String& archive_path, const MapType& map) const;
  };

} // namespace OpenMS

#endif // WITH_PARQUET

