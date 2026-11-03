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

#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/KERNEL/IMSExperimentRandomAccess.h>

#include <memory>
#include <string>

namespace OpenMS
{
  /**
    @brief Random access to the spectra of a Thermo RAW file without loading it.

    openFile() reads the scan headers once (no peaks) and keeps a small index: retention time,
    MS level and native ID of every selected scan, plus the precursors, which can only be
    reconstructed in acquisition order. getSpectrum() then decodes one scan on request, in
    any order, and gives the same spectrum as ThermoRawFile::load() with the same options at
    that index, including the precursor intensity estimate (the parent scan is decoded on
    demand; recent parents are cached).

    The options select the spectra as for load(): ms_levels, the RT, precursor m/z and scan
    ranges define which scans are indexed. Spectrum indices refer to the selected spectra.
    The source file checksum (Options::checksum) is never computed, as it reads the whole file.

    With Options::all_detectors, the detector chromatograms and the PDA spectra are read by
    openFile() and kept in memory (they are small); PDA spectra are merged into the spectrum
    order by retention time as in load().

    @code
    OnDiscThermoRawExperiment exp;
    ThermoRawFile::Options options;
    options.ms_levels = {2};
    exp.openFile("run.raw", options);
    Size first = exp.findFirstRT(600.0);
    for (Size i = first; i < exp.size() && exp.getSpectrumRT(i) < 660.0; ++i)
    {
      MSSpectrum spectrum = exp.getSpectrum(i);
    }
    @endcode

    @note An instance holds one handle on the file and is not thread-safe. Copies (and
    clone()) open a new handle on the already open file and share the index, so give each
    thread its own copy, e.g. with
    @code
    #pragma omp parallel for firstprivate(exp)
    @endcode

    @ingroup Kernel
  */
  class OPENMS_DLLAPI OnDiscThermoRawExperiment :
    public IMSExperimentRandomAccess
  {
  public:
    using Options = ThermoRawFile::Options;

    /// Default constructor; use openFile()
    OnDiscThermoRawExperiment();

    /// Shares the index of @p source and opens a new handle on the file
    OnDiscThermoRawExperiment(const OnDiscThermoRawExperiment& source);

    /// Move constructor
    OnDiscThermoRawExperiment(OnDiscThermoRawExperiment&& source) noexcept;

    /// Shares the index of @p source and opens a new handle on the file
    OnDiscThermoRawExperiment& operator=(const OnDiscThermoRawExperiment& source);

    /// Move assignment
    OnDiscThermoRawExperiment& operator=(OnDiscThermoRawExperiment&& source) noexcept;

    /// Destructor
    ~OnDiscThermoRawExperiment() override;

    /**
      @brief Open @p path and index its scans.

      @param[in] path Thermo RAW file
      @param[in] options Scan selection and conversion options (as for ThermoRawFile::load())

      @throws Exception::FileNotFound if @p path does not exist
      @throws Exception::ParseError if the file cannot be read
    */
    void openFile(const std::string& path, const Options& options = Options());

    /// True after a successful openFile()
    bool isOpen() const;

    /// The file given to openFile()
    const std::string& getFilename() const;

    /// The options given to openFile() (checksum is always false)
    const Options& getOptions() const;

    Size getNrSpectra() const override;
    Size getNrChromatograms() const override;

    /**
      @brief The spectrum with index @p index, as ThermoRawFile::load() gives it.

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSSpectrum getSpectrum(Size index) override;

    /// Alias for getSpectrum()
    MSSpectrum operator[](Size index)
    {
      return getSpectrum(index);
    }

    /**
      @brief The spectrum with index @p index without peaks.

      Same as ThermoRawFile::load() with Options::metadata_only, i.e. without the precursor
      intensity estimate and with the vendor's base peak values. Still reads the scan header.

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSSpectrum getSpectrumMetadata(Size index) override;

    /**
      @brief Peaks (and the charge array if requested) of spectrum @p index.

      Skips the scan header and the precursors; only RT, MS level and native ID are set.

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSSpectrum getSpectrumPeaks(Size index) override;

    /**
      @brief The spectrum with native ID @p native_id.

      @throws Exception::ElementNotFound if no indexed spectrum has this native ID
    */
    MSSpectrum getSpectrumByNativeId(const std::string& native_id);

    double getSpectrumRT(Size index) const override;
    UInt getSpectrumMSLevel(Size index) const override;

    /// Native ID of spectrum @p index
    const std::string& getSpectrumNativeID(Size index) const;

    /// Scan number of spectrum @p index (0 for PDA spectra)
    int getScanNumber(Size index) const;

    /// Index of the spectrum of scan @p scan, or getNrSpectra() if the scan is not indexed
    Size getSpectrumIndexByScanNumber(int scan) const;

    bool isSortedByRT() const override;
    Size findFirstRT(double rt) const override;

    /**
      @brief Chromatogram @p index: the TIC first (unless Options::skip_chromatograms), then detector chromatograms.

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSChromatogram getChromatogram(Size index) override;

    std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const override;

    std::unique_ptr<IMSExperimentRandomAccess> clone() const override;

  private:
    struct Index;
    struct Handle;
    void checkIndex_(Size index) const;
    std::shared_ptr<const Index> index_;
    std::unique_ptr<Handle> handle_;
  };
} // namespace OpenMS

#endif // WITH_THERMO_RAW
