// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/CachedMzML.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief @ref OpenSwath::ISpectrumAccess implementation backed by an
           on-disk @ref OpenMS::CachedmzML cache.

    Combines the in-memory index of @ref OpenMS::CachedmzML with the
    OpenSWATH spectrum-access interface: spectrum and chromatogram
    metadata are served from the in-memory @c MSExperiment, and the
    binary payloads are read on demand from the side-car
    @c .mzML.cached file.

    @note This implementation is @b not thread-safe: a single file
          stream is maintained and seeked per request, so concurrent
          calls would race on the stream position. Callers must
          serialise access externally.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SpectrumAccessOpenMSCached :
    public OpenSwath::ISpectrumAccess,
    public OpenMS::CachedmzML
  {

public:
    typedef OpenMS::PeakMap MSExperimentType;
    typedef OpenMS::MSSpectrum MSSpectrumType;

    /**
      @brief Open the cached @c .mzML pair at @p filename.

      Delegates to @ref OpenMS::CachedmzML's filename constructor: the
      metadata is loaded from @p filename and the side-car index for
      @p filename + @c ".cached" is built; subsequent accessor calls
      read the payloads from that side-car on demand.

      @param[in] filename Path to the @c .mzML metadata file; the
                          side-car is expected next to it at
                          @p filename + @c ".cached".

      @throws Exception::FileNotFound when either file cannot be opened.
      @throws Exception::ParseError   when the metadata cannot be
                                      parsed or the side-car index
                                      cannot be built.
    */
    explicit SpectrumAccessOpenMSCached(const String& filename);

    /// Destructor; the inherited @ref OpenMS::CachedmzML destructor closes the side-car stream.
    ~SpectrumAccessOpenMSCached() override;

    /**
      @brief Copy constructor.

      Each copy keeps its own file handle on the @c .mzML.cached
      side-car (a fresh stream is opened on the same path by the
      inherited @ref OpenMS::CachedmzML copy constructor). The
      in-memory metadata and side-car index are duplicated; the
      spectrum and chromatogram binary payloads are not.

      @param[in] rhs Source accessor to copy.
    */
    SpectrumAccessOpenMSCached(const SpectrumAccessOpenMSCached & rhs);

    /**
      @brief Return a clone of this accessor as a new @ref OpenSwath::ISpectrumAccess.

      Equivalent to copy-constructing through @ref SpectrumAccessOpenMSCached(const SpectrumAccessOpenMSCached&);
      the clone reads payloads lazily from its own stream on the same
      side-car file.

      @return Shared pointer to the cloned accessor.
    */
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    /**
      @brief Read one spectrum from the side-car file.

      @param[in] id Spectrum index in @c [0, getNrSpectra()). Out-of-range
                    access is a programming error (checked in debug builds).
      @return Decoded spectrum.

      @throws Exception::ParseError when @c seekg on the side-car stream
                                    fails (e.g. an overflow on a 32-bit
                                    system reading a > 2 GB cache file).
    */
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    /**
      @brief Return the RT and MS level of one spectrum.

      Only @c RT and @c ms_level are filled on the returned
      @c SpectrumMeta; for the full spectrum settings use
      @ref getSpectraMetaInfo.

      @param[in] id Spectrum index in @c [0, getNrSpectra()). Out-of-range
                    access is a programming error (checked in debug builds).
      @return @c SpectrumMeta carrying only @c RT and @c ms_level.
    */
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    /**
      @brief Indices of cached spectra whose RT lies in @c [RT - deltaRT, RT + deltaRT].

      Walks the metadata @c MSExperiment forward from the first spectrum
      with @c RT >= @c RT - @c deltaRT.

      @param[in] RT      Centre of the RT window.
      @param[in] deltaRT Half-width of the RT window. Must be non-negative;
                         a negative value is a programming error
                         (checked in debug builds).
      @return Indices into the in-memory metadata for spectra inside the window.
    */
    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    /// Number of spectra in the loaded metadata.
    size_t getNrSpectra() const override;

    /**
      @brief Return the full @ref SpectrumSettings of one spectrum.

      Unlike @ref getSpectrumMetaById (which returns only @c RT and
      @c ms_level), this returns the full per-spectrum metadata read
      at construction time.

      @param[in] id Spectrum index. No bounds check is performed.
      @return Copy of the stored @c SpectrumSettings for that spectrum.
    */
    SpectrumSettings getSpectraMetaInfo(int id) const;

    /**
      @brief Read one chromatogram from the side-car file.

      @param[in] id Chromatogram index in @c [0, getNrChromatograms()).
                    Out-of-range access is a programming error (checked
                    in debug builds).
      @return Decoded chromatogram.

      @throws Exception::ParseError when @c seekg on the side-car
                                    stream fails.
    */
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    /// Number of chromatograms in the loaded metadata.
    size_t getNrChromatograms() const override;

    /**
      @brief Return the @ref ChromatogramSettings of one chromatogram.

      @param[in] id Chromatogram index in @c [0, getNrChromatograms()).
                    Out-of-range access is a programming error (checked
                    in debug builds).
      @return Copy of the stored @c ChromatogramSettings.
    */
    ChromatogramSettings getChromatogramMetaInfo(int id) const;

    /**
      @brief Native id of one cached chromatogram.

      @param[in] id Chromatogram index in @c [0, getNrChromatograms()).
                    Out-of-range access is a programming error (checked
                    in debug builds).
      @return Native-id string from the stored metadata.
    */
    std::string getChromatogramNativeID(int id) const override;
  };

} //end namespace

