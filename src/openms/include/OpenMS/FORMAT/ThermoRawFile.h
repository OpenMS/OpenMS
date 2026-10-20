// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <limits>
#include <vector>

namespace OpenMS
{
  namespace Interfaces
  {
    class IMSDataConsumer;
  }

  /**
    @brief Reader for Thermo RAW files via the openms-thermo-bridge.

    Uses the openms-thermo-bridge library to access Thermo Fisher
    RawFileReader through the .NET host runtime. Reads spectra (MS1 and MSn)
    including retention times, precursor information, and instrument metadata.

    By default, scans are returned as they were acquired (profile or centroided);
    Options::centroid applies Thermo's peak picking. FileHandler::loadExperiment(),
    through which most TOPP tools read .raw input, sets it.

    RawFileReader reading tool. Copyright &copy; 2016 by Thermo Fisher Scientific, Inc.
    All rights reserved.

    Requires a .NET 8 (or newer) runtime on the machine and the
    openms-thermo-bridge managed runtime files (ThermoWrapperManaged.dll, its
    runtimeconfig.json and the Thermo CommonCore assemblies). These are looked
    up in this order:
      1. the directory named by the OPENMS_THERMO_MANAGED_DIR environment variable
         (skipped with a warning if it does not hold the files),
      2. <OpenMS share dir>/openms_thermo_bridge/managed (installed by OpenMS,
         also shipped inside pyOpenMS wheels),
      3. the bridge's own default, i.e. a 'managed' directory next to the
         openms_thermo_bridge shared library.
    Set DOTNET_ROOT if the .NET runtime is installed in a non-standard location.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ThermoRawFile : public ProgressLogger
  {
  public:
    /// Optional data exported by ThermoRawFileParser, plus lossless vendor metadata.
    struct Options
    {
      bool centroid = false; ///< Preserve acquired representation by default; true matches TRFP peak picking.
      bool charge_data = false; ///< Instrument-assigned centroid peak charge array.
      bool noise_data = false; ///< Noise/baseline arrays and their independent mass coordinates.
      bool all_detectors = false; ///< UV/PDA/Analog/MSAnalog traces and PDA spectra.
      bool preserve_trailers = true; ///< Retain all original scan trailer label/value pairs.
      bool instrument_methods = true; ///< Retain embedded method text once per run.
      bool checksum = true; ///< Compute source SHA-1 for mzML provenance (reads the whole file once more).

      /** @name Selection of MS scans
          Scans that are not selected are not decoded. Precursors of the selected scans
          (including the precursor intensity estimated from the parent scan) are the same as
          in a full load, also if the parent scan itself is not selected.
      */
      //@{
      int first_scan = 0; ///< First scan number to load (0: from the first scan).
      int last_scan = 0; ///< Last scan number to load (0: up to the last scan).
      double rt_min = -std::numeric_limits<double>::infinity(); ///< Load scans with rt_min <= RT (seconds) ...
      double rt_max = std::numeric_limits<double>::infinity(); ///< ... and RT < rt_max.
      std::vector<int> ms_levels; ///< MS levels to load (empty: all). PDA spectra (all_detectors) have level 0.
      /// Skip MSn scans whose first precursor m/z lies outside [precursor_mz_min, precursor_mz_max) (MS1 is kept).
      double precursor_mz_min = -std::numeric_limits<double>::infinity();
      double precursor_mz_max = std::numeric_limits<double>::infinity(); ///< see precursor_mz_min
      //@}

      /**
        @brief Load only spectrum metadata: no peaks, charge or noise arrays are read.

        Values that are derived from the peaks are not available then: the precursor intensity,
        "lowest/highest observed m/z" and, for centroided data, the base peak of the centroids
        (the vendor scan statistics base peak is reported instead).
      */
      bool metadata_only = false;
      bool skip_chromatograms = false; ///< Do not load the TIC and detector chromatograms.
    };

    /// Returns current loading options.
    const Options& getOptions() const { return options_; }
    /// Set loading options.
    void setOptions(const Options& options) { options_ = options; }
    /**
      @brief Load a Thermo RAW file into an MSExperiment.

      Reads all scans from the RAW file and populates the experiment with
      spectra, retention times, MS levels, precursor information (for MSn),
      and source file metadata.

      @param[in] path Path to the .raw file, or to a ZIP archive ('.raw.zip') that holds exactly one .raw file.
                      RawFileReader reads files on disk only, so that file is extracted into a temporary
                      directory, which is deleted when the call returns.
      @param[out] exp The experiment to populate

      @throws Exception::FileNotFound if the file does not exist
      @throws Exception::ParseError if the file cannot be read by the thermo bridge
    */
    void load(const std::string& path, MSExperiment& exp);

    /**
      @brief Stream a Thermo RAW file into a consumer, scan by scan, with constant memory.

      Gives the same spectra and chromatograms as load() with the same Options, without holding
      the run in memory. A first pass reads only the scan headers (no peaks) to count the
      selected spectra and to collect the instrument configurations, so that @p consumer gets,
      in this order:
        1. setExpectedSize() with the exact number of spectra and chromatograms,
        2. setExperimentalSettings() with the run-level settings, including all instrument configurations,
        3. consumeSpectrum() for every selected MS scan in acquisition (RT) order,
        4. consumeSpectrum() for the PDA spectra (Options::all_detectors, MS level 0),
        5. consumeChromatogram() for the TIC and then the detector chromatograms.

      Differences to load(): PDA spectra follow the MS scans instead of being sorted into them
      by RT, and the detector data (Options::all_detectors) is read into memory at once. The
      header pass also resolves the precursors only if a precursor m/z range is set.

      The precursor intensity of an MSn scan is estimated from its parent scan as in load(): the
      most recent converted parent candidates are kept, and a parent that was not converted or
      is no longer kept is decoded again.

      @param[in] path Path to the .raw file, or to a ZIP archive ('.raw.zip') that holds exactly one .raw file.
                      RawFileReader reads files on disk only, so that file is extracted into a temporary
                      directory, which is deleted when the call returns.
      @param[in] consumer Receives the run; must not be null

      @throws Exception::FileNotFound if the file does not exist
      @throws Exception::ParseError if the file cannot be read by the thermo bridge
    */
    void transform(const std::string& path, Interfaces::IMSDataConsumer* consumer);

  private:
    Options options_;
  };

} // namespace OpenMS

#endif // WITH_THERMO_RAW
