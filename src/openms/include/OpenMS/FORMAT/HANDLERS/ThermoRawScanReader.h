// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/ThermoRawFile.h>

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

    /// Data processing attached to every converted spectrum and chromatogram.
    std::shared_ptr<DataProcessing> dataProcessing() const;

    /// TIC chromatogram as reported by the vendor library (reads every scan header).
    MSChromatogram ticChromatogram() const;

    /// Detector (UV/PDA/analog) chromatograms and PDA spectra (Options::all_detectors).
    void readDetectors(std::vector<MSChromatogram>& chromatograms, std::vector<MSSpectrum>& spectra);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
} // namespace Internal
} // namespace OpenMS
