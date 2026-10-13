// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/IMSExperimentRandomAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @brief @ref OpenSwath::ISpectrumAccess on an on-disc experiment (e.g. OnDiscThermoRawExperiment or OnDiscMSExperiment).

    Spectra are read from the file on request (IMSExperimentRandomAccess::getSpectrumPeaks()), nothing
    is cached. The accessor can expose a subset of the spectra, e.g. those of one SWATH window, so that
    several views share one opened file and its index (as SpectrumAccessSqMass does). The visible
    spectra are ordered by retention time (ties keep the order of the experiment); id @c i refers to
    the @c i-th visible spectrum.

    Chromatograms are those of the experiment (the subset applies to spectra only).

    @note lightClone() and the copy constructor clone the experiment (IMSExperimentRandomAccess::clone()),
    i.e. open a new handle on the file, and share the index and the subset. Clones can be used
    concurrently from different threads; a single instance must not.

    @code
    auto experiment = std::make_shared<OnDiscThermoRawExperiment>();
    experiment->openFile("run.raw");
    std::vector<Size> window; // spectra of one SWATH window
    // ...
    OpenSwath::SpectrumAccessPtr sptr = std::make_shared<SpectrumAccessOnDisc>(experiment, window);
    @endcode

    @ingroup Analysis_OpenSWATH
  */
  class OPENMS_DLLAPI SpectrumAccessOnDisc :
    public OpenSwath::ISpectrumAccess
  {
  public:
    /**
      @brief Expose all spectra of @p experiment.

      @param[in] experiment The opened experiment; this accessor uses it, clones open new handles.
    */
    explicit SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment);

    /**
      @brief Expose only the spectra of @p experiment with the given indices.

      @param[in] experiment The opened experiment; this accessor uses it, clones open new handles.
      @param[in] indices Spectrum indices of @p experiment; they are ordered by retention time.

      @throws Exception::IndexOverflow if an index is out of range
    */
    SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment, std::vector<Size> indices);

    /// Copy constructor; clones the experiment (new file handle) and shares the index and subset
    SpectrumAccessOnDisc(const SpectrumAccessOnDisc& rhs);

    /// Destructor
    ~SpectrumAccessOnDisc() override;

    /// A copy that can be used from another thread (see the copy constructor)
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    /// Peaks and data arrays of visible spectrum @p id (float and integer arrays, named after the array)
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    /// Retention time and MS level of visible spectrum @p id
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    /**
      @brief Visible spectra with RT in [@p RT - @p deltaRT, @p RT + @p deltaRT].

      Returns the first visible spectrum with RT >= @p RT - @p deltaRT and all following ones up to
      @p RT + @p deltaRT, as SpectrumAccessOpenMS does (so @p deltaRT = 0 gives the first spectrum at or
      after @p RT).
    */
    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    /// Number of visible spectra
    size_t getNrSpectra() const override;

    /// Chromatogram @p id of the experiment
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    /// Number of chromatograms of the experiment
    size_t getNrChromatograms() const override;

    /// Native ID of chromatogram @p id
    std::string getChromatogramNativeID(int id) const override;

    /// Index in the experiment of visible spectrum @p id
    Size getExperimentIndex(int id) const;

    /// The experiment this accessor reads from
    const IMSExperimentRandomAccess& getExperiment() const;

  private:
    struct View
    {
      std::vector<Size> indices; ///< experiment index of each visible spectrum
      std::vector<double> rts;   ///< and its retention time (ascending)
    };
    std::shared_ptr<IMSExperimentRandomAccess> experiment_;
    std::shared_ptr<const View> view_;
  };
} // namespace OpenMS
