// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/OpenMSConfig.h>

#include <memory>

namespace OpenMS
{
  class ExperimentalSettings;

  /**
    @brief Random access to the spectra and chromatograms of an MS run that stays on disk.

    Backends keep a small index in memory (e.g. retention time and MS level of every spectrum)
    and decode spectra on request. Implemented by OnDiscMSExperiment (indexed mzML) and
    OnDiscThermoRawExperiment (Thermo RAW); further formats (e.g. Bruker TDF, imzML) can be
    added the same way. Algorithms that only need random access, such as the OpenSWATH
    adapter SpectrumAccessOnDisc, are written against this interface.

    @note An instance is not thread-safe, as it typically holds one file handle. Give each
    thread its own clone(); clones share the immutable index.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI IMSExperimentRandomAccess
  {
  public:
    virtual ~IMSExperimentRandomAccess();

    /// Number of spectra
    virtual Size getNrSpectra() const = 0;

    /// Number of chromatograms
    virtual Size getNrChromatograms() const = 0;

    /// Alias for getNrSpectra()
    Size size() const
    {
      return getNrSpectra();
    }

    /// True if there are no spectra
    bool empty() const
    {
      return getNrSpectra() == 0;
    }

    /**
      @brief The spectrum with index @p index, including peaks and metadata.

      @param[in] index Spectrum index (0 <= index < getNrSpectra())
      @throws Exception::IndexOverflow if @p index is out of range (backend dependent)
    */
    virtual MSSpectrum getSpectrum(Size index) = 0;

    /**
      @brief The spectrum with index @p index without peaks and data arrays.

      @param[in] index Spectrum index (0 <= index < getNrSpectra())
    */
    virtual MSSpectrum getSpectrumMetadata(Size index) = 0;

    /**
      @brief Peaks and data arrays of spectrum @p index; other metadata may be missing.

      Backends can skip decoding metadata that consumers of peak data alone (e.g. OpenSWATH)
      do not need. The default returns getSpectrum().

      @param[in] index Spectrum index (0 <= index < getNrSpectra())
    */
    virtual MSSpectrum getSpectrumPeaks(Size index);

    /// Retention time (in seconds) of spectrum @p index, from the index
    virtual double getSpectrumRT(Size index) const = 0;

    /// MS level of spectrum @p index, from the index
    virtual UInt getSpectrumMSLevel(Size index) const = 0;

    /// True if the spectra are sorted by ascending retention time
    virtual bool isSortedByRT() const = 0;

    /**
      @brief Index of the first spectrum with retention time >= @p rt.

      Binary search if the spectra are sorted by retention time, linear search otherwise.

      @param[in] rt Retention time in seconds
      @return the index, or getNrSpectra() if there is no such spectrum
    */
    virtual Size findFirstRT(double rt) const;

    /**
      @brief The chromatogram with index @p index.

      @param[in] index Chromatogram index (0 <= index < getNrChromatograms())
    */
    virtual MSChromatogram getChromatogram(Size index) = 0;

    /// Run-level metadata (instrument, source files, sample, ...)
    virtual std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const = 0;

    /**
      @brief A new instance on the same data with its own file handle (for another thread).

      Cheap: the index is shared, not copied.
    */
    virtual std::unique_ptr<IMSExperimentRandomAccess> clone() const = 0;

  protected:
    IMSExperimentRandomAccess() = default;
    IMSExperimentRandomAccess(const IMSExperimentRandomAccess&) = default;
    IMSExperimentRandomAccess& operator=(const IMSExperimentRandomAccess&) = default;
  };
} // namespace OpenMS
