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

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <memory>

namespace OpenMS
{
  /**
    @brief @ref OpenSwath::ISpectrumAccess implementation that holds the complete
    spectrum / chromatogram payload in memory.

    Constructed from any other @ref OpenSwath::ISpectrumAccess (file-backed, on-disk
    cached, etc.) — every spectrum and chromatogram is pulled through the source's
    interface at construction time and cached locally, so subsequent accesses do not
    touch the disk again. The view exposes exactly the same data as the source object.

    There is a fast-path optimisation when the source happens to be a
    @ref SpectrumAccessSqMass — the constructor recognises it via @c dynamic_cast and
    uses the batched @ref SpectrumAccessSqMass::getAllSpectra call instead of issuing
    per-spectrum reads (avoiding the per-SQLite-query overhead). For every other source
    the constructor falls back to a per-spectrum / per-chromatogram loop.

    @code
      OpenSwath::ISpectrumAccess * data_access;
      fillData(data_access); // assume that data_access points to some data
      OpenSwath::ISpectrumAccess * in_memory_data_access = new SpectrumAccessOpenMSInMemory(*data_access);
    @endcode

    After executing this code, two variables exist: @c data_access which provides access
    to the original data (in one of multiple ways which is not transparent to the user),
    and @c in_memory_data_access which exposes the same data with the guarantee that it
    is held in memory and not re-read from disk.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SpectrumAccessOpenMSInMemory :
    public OpenSwath::ISpectrumAccess
  {
public:
    typedef OpenMS::PeakMap MSExperimentType;
    typedef OpenMS::MSSpectrum MSSpectrumType;
    typedef OpenMS::MSChromatogram MSChromatogramType;

    /**
      @brief Construct by reading every spectrum and chromatogram from @p origin into memory.

      For a @ref SpectrumAccessSqMass source, takes the batched
      @ref SpectrumAccessSqMass::getAllSpectra fast path. For every other source, iterates
      the spectra (@c 0..getNrSpectra()-1) and chromatograms
      (@c 0..getNrChromatograms()-1) and pulls each through the per-id getter. After
      construction, @p origin is no longer referenced.

      @param[in] origin Source accessor; iterated once during construction and not retained.
    */
    explicit SpectrumAccessOpenMSInMemory(OpenSwath::ISpectrumAccess & origin);

    /// Destructor
    ~SpectrumAccessOpenMSInMemory() override;

    /**
      @brief Copy constructor (light copy).

      Only the @c SpectrumPtr / @c ChromatogramPtr vectors are copied, not the underlying
      spectrum / chromatogram payloads (those stay shared via @c shared_ptr).

      @param[in] rhs Source accessor to copy.
    */
    SpectrumAccessOpenMSInMemory(const SpectrumAccessOpenMSInMemory & rhs);

    /**
      @brief Return a clone of this accessor as a new @ref OpenSwath::ISpectrumAccess pointer.

      Equivalent to copy-constructing via the copy constructor — the clone shares the
      underlying spectrum / chromatogram @c shared_ptr payloads, so it is safe to hand to
      another thread without re-reading the source.

      @return Shared pointer to the cloned accessor.
    */
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    /**
      @brief Direct in-memory lookup of one spectrum by index.

      @param[in] id Spectrum index (0-based; OPENMS_PRECONDITION asserts @c 0 <= id < @c getNrSpectra()).
      @return Cached @c SpectrumPtr for the requested spectrum.
    */
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    /**
      @brief Direct in-memory lookup of one spectrum's metadata by index.

      @param[in] id Spectrum index (0-based; OPENMS_PRECONDITION asserts @c 0 <= id < @c getNrSpectra()).
      @return Cached @c SpectrumMeta for the requested spectrum.
    */
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    /**
      @brief Indices of cached spectra whose RT lies in [@p RT - @p deltaRT, @p RT + @p deltaRT].

      Uses @c std::lower_bound on the cached metadata (assumed RT-sorted from the source)
      to seek to the first matching spectrum, then walks forward while the RT remains in
      range.

      @param[in] RT      Centre of the RT window.
      @param[in] deltaRT Half-width of the RT window (OPENMS_PRECONDITION asserts @c deltaRT >= 0).
      @return Indices into the in-memory view of spectra whose RT falls in the window.
    */
    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    /// Number of spectra cached in memory.
    size_t getNrSpectra() const override;

    /**
      @brief Direct in-memory lookup of one chromatogram by index.

      @param[in] id Chromatogram index (0-based; OPENMS_PRECONDITION asserts @c 0 <= id < @c getNrChromatograms()).
      @return Cached @c ChromatogramPtr for the requested chromatogram.
    */
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    /// Number of chromatograms cached in memory.
    size_t getNrChromatograms() const override;

    /**
      @brief Native id of one cached chromatogram.

      @param[in] id Chromatogram index (0-based; OPENMS_PRECONDITION asserts @c 0 <= id < @c getNrChromatograms()).
      @return Native-id string captured at construction time.
    */
    std::string getChromatogramNativeID(int id) const override;

private:

    std::vector< OpenSwath::SpectrumPtr > spectra_;       ///< Spectrum payloads pulled from the source at construction; shared with any @ref lightClone derivative via @c shared_ptr
    std::vector< OpenSwath::SpectrumMeta > spectra_meta_; ///< Parallel metadata cache (native id, RT, MS level) used by @ref getSpectrumMetaById and the @ref getSpectraByRT binary search

    std::vector< OpenSwath::ChromatogramPtr > chromatograms_; ///< Chromatogram payloads cached at construction; shared with clones
    std::vector< std::string > chromatogram_ids_;             ///< Parallel chromatogram native-id cache used by @ref getChromatogramNativeID

  };

} //end namespace OpenMS

