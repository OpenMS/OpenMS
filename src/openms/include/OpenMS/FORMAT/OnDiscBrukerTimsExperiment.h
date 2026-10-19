// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/KERNEL/IMSExperimentRandomAccess.h>

#include <memory>
#include <string>

namespace OpenMS
{
  namespace Internal
  {
    class BrukerTimsRun;
    class BrukerTimsSpectrumReader;
  }

  /**
    @brief Random access to the spectra of a Bruker timsTOF .d directory without loading it.

    openFile() lists the spectra that BrukerTimsFile::load() gives with the same configuration, in the same order:
    their retention time, MS level and native ID. It reads the frame table and analysis.tdf, not the peaks.
    getSpectrum() then builds one spectrum on request, in any order, from the frames it needs: one frame for an MS1
    spectrum, a frame spectrum (export mode FRAME) or a DIA window spectrum, the frames of the precursor for a DDA MS2
    spectrum, and the neighboring frames with aggregation over retention time.

    The configuration selects and shapes the spectra as for load(). The one exception is Config::calibrate with the
    linear m/z model, whose recalibration is fitted to all precursors of the run: openFile() rejects it.

    A few spectra turn out to have no peaks only when they are built, and load() leaves those out: a DDA precursor or
    a DIA window without peaks (or without centroids), and a DDA MS2 spectrum whose ion mobility lies outside
    Config::im_min/im_max. They are listed here, and getSpectrum() gives them without peaks. So load() gives exactly
    the spectra of this experiment, minus the MS2 spectra without peaks that getSpectrum() gives for those. All other
    spectra are those of load(), at the same position relative to each other.

    @code
    OnDiscBrukerTimsExperiment exp;
    BrukerTimsFile::Config config;
    config.load_ms1 = false;
    exp.openFile("run.d", config);
    Size first = exp.findFirstRT(600.0);
    for (Size i = first; i < exp.size() && exp.getSpectrumRT(i) < 660.0; ++i)
    {
      MSSpectrum spectrum = exp.getSpectrum(i);
    }
    @endcode

    @note An instance holds one handle on the data and is not thread-safe. Copies (and clone()) share the list of
    spectra and open their own handle when they read their first spectrum, so give each thread its own copy, e.g. with
    @code
    #pragma omp parallel for firstprivate(exp)
    @endcode
    With the converters of the Bruker SDK (Config::tims_calibration_strategy), copies read one at a time.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI OnDiscBrukerTimsExperiment :
    public IMSExperimentRandomAccess
  {
  public:
    using Config = BrukerTimsFile::Config;

    /// Default constructor; use openFile()
    OnDiscBrukerTimsExperiment();

    /// Shares the list of spectra of @p source; opens its own handle on the data
    OnDiscBrukerTimsExperiment(const OnDiscBrukerTimsExperiment& source);

    /// Move constructor
    OnDiscBrukerTimsExperiment(OnDiscBrukerTimsExperiment&& source) noexcept;

    /// Shares the list of spectra of @p source; opens its own handle on the data
    OnDiscBrukerTimsExperiment& operator=(const OnDiscBrukerTimsExperiment& source);

    /// Move assignment
    OnDiscBrukerTimsExperiment& operator=(OnDiscBrukerTimsExperiment&& source) noexcept;

    /// Destructor
    ~OnDiscBrukerTimsExperiment() override;

    /**
      @brief Open @p path and list its spectra

      @param[in] path A Bruker .d directory, or a '.d.zip' archive of one (unpacked into a temporary directory that
                 this experiment and its copies keep)
      @param[in] config Which spectra to give and how to build them, as for BrukerTimsFile::load()

      @throws Exception::IllegalArgument for Config::calibrate with the linear m/z model (see above)
      @throws the exceptions of BrukerTimsFile::load(), e.g. Exception::ParseError for data that cannot be read
    */
    void openFile(const std::string& path, const Config& config = Config());

    /// True after a successful openFile()
    bool isOpen() const;

    /// The path given to openFile()
    const std::string& getFilename() const;

    /// The configuration given to openFile()
    const Config& getConfig() const;

    Size getNrSpectra() const override;

    /// No chromatograms (as load())
    Size getNrChromatograms() const override;

    /**
      @brief The spectrum with index @p index, as BrukerTimsFile::load() gives it (see above for spectra that load()
      leaves out)

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSSpectrum getSpectrum(Size index) override;

    /// Alias for getSpectrum()
    MSSpectrum operator[](Size index)
    {
      return getSpectrum(index);
    }

    /**
      @brief The spectrum with index @p index without peaks, from the list of spectra (no data is read)

      Retention time, MS level, native ID, spectrum type, ion mobility unit, and for MS2 spectra the precursor (DDA: with
      its meta value 'bruker_precursor_id'; DIA: the isolation window, with the meta values 'ion mobility lower limit'
      and 'ion mobility upper limit'). Not set are the values that come from the peaks: the ion mobility of a DDA MS2
      spectrum and its precursor, and the ion mobility peak type.

      @throws Exception::IndexOverflow if @p index is out of range
    */
    MSSpectrum getSpectrumMetadata(Size index) override;

    /**
      @brief The spectrum with native ID @p native_id

      @throws Exception::ElementNotFound if no spectrum has this native ID
    */
    MSSpectrum getSpectrumByNativeId(const std::string& native_id);

    double getSpectrumRT(Size index) const override;
    UInt getSpectrumMSLevel(Size index) const override;

    /// Native ID of spectrum @p index
    const std::string& getSpectrumNativeID(Size index) const;

    /// Index of the spectrum with native ID @p native_id, or getNrSpectra() if there is none
    Size getSpectrumIndexByNativeId(const std::string& native_id) const;

    bool isSortedByRT() const override;
    Size findFirstRT(double rt) const override;

    /// @throws Exception::IndexOverflow (there are no chromatograms)
    MSChromatogram getChromatogram(Size index) override;

    std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const override;

    std::unique_ptr<IMSExperimentRandomAccess> clone() const override;

  private:
    struct Index;
    void checkIndex_(Size index) const;
    std::shared_ptr<const Index> index_;
    std::unique_ptr<Internal::BrukerTimsSpectrumReader> reader_;
  };
} // namespace OpenMS

#endif // WITH_OPENTIMS
