// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  class DataProcessing;
  class ExperimentalSettings;
  class MSChromatogram;
  class MSSpectrum;

namespace Internal
{
  /**
    @brief Scan-by-scan conversion of a Thermo RAW file into OpenMS spectra (private).

    Shared by ThermoRawFile::load() and by later streaming and random access readers. A reader
    opens the RAW file once; all calls must come from one thread.

    Usage:
      1. fillSettings() writes the run-level settings (source file, sample, instrument, data
         processing).
      2. header() must be called for every MS scan in acquisition order, from firstScan() on,
         because the precursor of an MSn scan is reconstructed from the scans before it.
         It reads only the scan header, unless the caller asks for the peaks as well because
         it will convert the scan anyway (one vendor read instead of two).
      3. convert() turns the scan whose header was read last into an MSSpectrum. Scans that
         are not needed can simply be skipped (no peaks are read for them).

    The precursor intensity estimate of an MSn scan sums peaks of its parent scan. convert()
    asks the ParentLookup for the parent spectrum first (e.g. the already loaded experiment);
    if it returns nullptr, the reader decodes the parent peaks itself, so filtered loads give
    the same precursors as a full load.
  */
  class OPENMS_DLLAPI ThermoRawScanReader
  {
  public:
    using Options = ThermoRawFile::Options;

    /// Scan header values needed to select scans before converting them.
    struct ScanInfo
    {
      int scan = 0;          ///< scan number
      int ms_level = 0;      ///< MS level (0 for unknown or non-MS controllers)
      double rt = 0.0;       ///< retention time in seconds
      double tic = 0.0;      ///< total ion current (NaN if the vendor reports none)
      double precursor_mz = 0.0; ///< m/z of the first precursor (NaN if there is none)
    };

    /**
      @brief Everything of a scan that depends on the scans before it (see header()).

      With a record, Accessor::convert() converts the scan in any order and on any thread.
    */
    struct ScanRecord
    {
      ScanInfo info;
      std::vector<ThermoPrecursor> precursors; ///< resolved in acquisition order
      std::string configuration_id;            ///< instrument configuration id
      std::string native_id;
    };

    /**
      @brief Converts scans independently of each other (random access).

      Created by ThermoRawScanReader::accessor() once the records of the scans have been taken.
      Each accessor has its own handle on the file: one accessor must only be used by one thread
      at a time, but distinct accessors (e.g. from clone()) may be used concurrently. Accessors
      may outlive the reader that created them.
    */
    class OPENMS_DLLAPI Accessor
    {
    public:
      ~Accessor();
      Accessor(const Accessor&) = delete;
      Accessor& operator=(const Accessor&) = delete;

      /**
        @brief Convert the scan of @p record, as ThermoRawScanReader::convert() would.

        The parent spectra for the precursor intensity estimate are decoded on demand and
        the most recent ones are kept in a small cache.

        @param[in] record Record of the scan from ThermoRawScanReader::record()
        @param[out] spectrum The converted spectrum
        @param[in] peaks If false, no peaks or auxiliary arrays are read (metadata only)
      */
      void convert(const ScanRecord& record, MSSpectrum& spectrum, bool peaks);

      /// Only the peaks (and the charge array if requested by the options) of @p scan; @p spectrum is reset.
      void readPeaks(int scan, MSSpectrum& spectrum);

      /// A new accessor on the same file; safe to call concurrently from several threads.
      std::unique_ptr<Accessor> clone() const;

    private:
      friend class ThermoRawScanReader;
      struct Impl;
      explicit Accessor(std::unique_ptr<Impl> impl);
      std::unique_ptr<Impl> impl_;
    };

    /// Returns the converted (peak-sorted) spectrum of a scan or nullptr if it is not available.
    using ParentLookup = std::function<const MSSpectrum*(int scan)>;

    /**
      @brief Open @p path.

      @param[in] path Thermo RAW file (symbolic links are resolved for the vendor library)
      @param[in] options Conversion options; ms_levels / rt / scan ranges are not applied here

      @throws Exception::FileNotFound if @p path does not exist
      @throws std::exception if the bridge cannot open the file
    */
    ThermoRawScanReader(const std::string& path, const Options& options);
    ~ThermoRawScanReader();
    ThermoRawScanReader(const ThermoRawScanReader&) = delete;
    ThermoRawScanReader& operator=(const ThermoRawScanReader&) = delete;

    /// True if the file contains a mass spectrometer controller.
    bool hasMSData() const;
    /// First MS scan number (only valid if hasMSData()).
    int firstScan() const;
    /// Last MS scan number (only valid if hasMSData()).
    int lastScan() const;

    /// Write run-level settings into @p settings (instrument configurations are added by convert() and addInstrumentConfigurations()).
    void fillSettings(ExperimentalSettings& settings) const;

    /**
      @brief Read the header of @p scan; must be called for every scan in increasing order.

      @param[in] scan Scan number, the successor of the previous call's scan
      @param[in] with_peaks Also read the peaks for a following convert() call. Saves a second
                 vendor read of the scan when the caller already knows it will convert it.
    */
    const ScanInfo& header(int scan, bool with_peaks = false);

    /**
      @brief Read the following scans up to @p last_scan on worker threads ahead of header().

      Uses up to three threads besides the calling one, limited by the OpenMP thread count
      (none if it is 1). Output is unchanged; header() takes the scans from the read-ahead
      buffer, which holds a bounded number of scans.

      @param[in] last_scan Last scan to read ahead
      @param[in] peaks_from Scans from this number on are read with peaks (as header(scan, true))
    */
    void readAhead(int last_scan, int peaks_from);

    /**
      @brief Visit the headers of scans @p first to @p last ahead of header(), without peaks.

      A cheap pass, e.g. to count the scans a filtered conversion selects before streaming them
      or to index a file for random access. Instrument configurations are registered in the
      same order as header() would register them (so ids match a single pass), and the state
      of header() and convert() is not changed otherwise.

      @param[in] first First scan number; use firstScan() if precursors are resolved
      @param[in] last Last scan number
      @param[in] precursors Also set ScanInfo::precursor_mz. Reads the trailers and resolves
                 the precursors from @p first on (slower); otherwise precursor_mz is NaN.
      @param[in] visit Called for every scan in order; returning false stops the pass
    */
    void scanHeaders(int first, int last, bool precursors, const std::function<bool(const ScanInfo&)>& visit);

    /**
      @brief Convert the scan of the last header() call.

      @param[out] spectrum The converted spectrum (peaks sorted by m/z)
      @param[in,out] settings Receives instrument configurations as they appear
      @param[in] parent Lookup of already converted parent spectra (may be empty)
      @param[in] peaks If false, no peaks or auxiliary arrays are read (metadata only)
    */
    void convert(MSSpectrum& spectrum, ExperimentalSettings& settings, const ParentLookup& parent, bool peaks);

    /**
      @brief Add the instrument configurations of all scans read by header() so far to @p settings.

      Configuration ids are assigned in the order of all scans read, independent of which scans
      are converted. Call this after the last header() so that configurations of scans that were
      not converted are also listed, as in a load without scan selection.
    */
    void addInstrumentConfigurations(ExperimentalSettings& settings);

    /// Record of the scan whose header was read last (by header()).
    ScanRecord record() const;

    /**
      @brief Random access converter for scans whose records were taken.

      Call it after the header() pass: precursor intensity estimates use parents from the first
      scan visited by header() on, as a sequential conversion does.
    */
    std::unique_ptr<Accessor> accessor() const;

    /// Data processing attached to every converted spectrum and chromatogram.
    std::shared_ptr<DataProcessing> dataProcessing() const;

    /// TIC chromatogram as reported by the vendor library (reads every scan header).
    MSChromatogram ticChromatogram() const;

    /// Detector (UV/PDA/analog) chromatograms and PDA spectra (Options::all_detectors); may be called at any time.
    void readDetectors(std::vector<MSChromatogram>& chromatograms, std::vector<MSSpectrum>& spectra);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
  /// Scan selection by ThermoRawFile::Options, shared by load(), transform() and OnDiscThermoRawExperiment.
  struct ThermoRawScanSelection
  {
    const ThermoRawFile::Options& options;

    bool level(int level) const
    {
      return options.ms_levels.empty() || std::find(options.ms_levels.begin(), options.ms_levels.end(), level) != options.ms_levels.end();
    }
    bool rt(double rt) const { return rt >= options.rt_min && rt < options.rt_max; }
    bool precursorFiltered() const
    {
      return options.precursor_mz_min > -std::numeric_limits<double>::infinity() || options.precursor_mz_max < std::numeric_limits<double>::infinity();
    }
    /// True if scans are selected by values of their header.
    bool headerFiltered() const
    {
      return ! options.ms_levels.empty() || precursorFiltered() || options.rt_min > -std::numeric_limits<double>::infinity()
             || options.rt_max < std::numeric_limits<double>::infinity();
    }
    /// Last scan to read: scans of one controller are in acquisition order, so no later scan can be selected.
    int lastScan(int last) const { return options.last_scan > 0 ? std::min(last, options.last_scan) : last; }
    /// True if no scan from @p info on can be selected (checked after reading its header).
    bool beyondRT(const ThermoRawScanReader::ScanInfo& info) const { return info.rt >= options.rt_max; }
    bool selected(const ThermoRawScanReader::ScanInfo& info) const
    {
      if (info.scan < options.first_scan || ! rt(info.rt) || ! level(info.ms_level)) { return false; }
      return info.ms_level <= 1 || std::isnan(info.precursor_mz)
             || (info.precursor_mz >= options.precursor_mz_min && info.precursor_mz < options.precursor_mz_max);
    }
  };

} // namespace Internal
} // namespace OpenMS
