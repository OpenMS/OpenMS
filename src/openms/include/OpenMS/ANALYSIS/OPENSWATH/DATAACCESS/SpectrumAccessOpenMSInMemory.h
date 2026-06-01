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

    Constructed from any other @ref OpenSwath::ISpectrumAccess (file-backed,
    on-disk cached, etc.); every spectrum and chromatogram is read through the
    source's interface during construction and cached locally, so subsequent
    accesses do not touch the disk again. The view exposes exactly the same
    data as the source. After construction the source is no longer referenced.

    @code
      OpenSwath::ISpectrumAccess * data_access;
      fillData(data_access); // assume data_access now points at some data
      OpenSwath::ISpectrumAccess * in_memory_data_access = new SpectrumAccessOpenMSInMemory(*data_access);
    @endcode

    After this snippet, @c data_access still serves the original data (its
    backing store is opaque to the caller), while @c in_memory_data_access
    serves the same data with the guarantee that it is held in memory and
    not re-read from disk.

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
      @brief Construct by reading every spectrum and chromatogram from
             @p origin into memory.

      After construction the accessor exposes the same spectra and
      chromatograms as @p origin, served entirely from memory; @p origin
      itself is no longer referenced.

      @param[in] origin Source accessor; read once during construction and
                        not retained.
    */
    explicit SpectrumAccessOpenMSInMemory(OpenSwath::ISpectrumAccess & origin);

    /// Destructor
    ~SpectrumAccessOpenMSInMemory() override;

    /**
      @brief Copy constructor (light copy).

      The copy shares the underlying spectrum / chromatogram payloads with
      @p rhs; only the per-element handles are duplicated. No data is
      re-read from the original source.

      @param[in] rhs Source accessor to copy.
    */
    SpectrumAccessOpenMSInMemory(const SpectrumAccessOpenMSInMemory & rhs);

    /**
      @brief Return a clone of this accessor as a new @ref OpenSwath::ISpectrumAccess.

      The clone shares this accessor's spectrum / chromatogram payloads, so
      it is cheap to create and safe to hand to another thread without
      re-reading the source.

      @return Shared pointer to the cloned accessor.
    */
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    /**
      @brief Look up one spectrum by index.

      @param[in] id Spectrum index in @c [0, getNrSpectra()). Out-of-range
                    access is a programming error (checked in debug builds).
      @return Cached @c SpectrumPtr for the requested spectrum.
    */
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    /**
      @brief Look up one spectrum's metadata by index.

      @param[in] id Spectrum index in @c [0, getNrSpectra()). Out-of-range
                    access is a programming error (checked in debug builds).
      @return Cached @c SpectrumMeta for the requested spectrum.
    */
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    /**
      @brief Indices of cached spectra whose RT lies in [@p RT - @p deltaRT, @p RT + @p deltaRT].

      Relies on the cached metadata being sorted by RT (inherited from the
      source); the lookup is logarithmic in the number of spectra.

      @param[in] RT      Centre of the RT window.
      @param[in] deltaRT Half-width of the RT window. Must be non-negative;
                         a negative value is a programming error (checked
                         in debug builds).
      @return Indices into the in-memory view of spectra whose RT falls in
              the window.
    */
    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    /// Number of spectra cached in memory.
    size_t getNrSpectra() const override;

    /**
      @brief Look up one chromatogram by index.

      @param[in] id Chromatogram index in @c [0, getNrChromatograms()).
                    Out-of-range access is a programming error (checked in
                    debug builds).
      @return Cached @c ChromatogramPtr for the requested chromatogram.
    */
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    /// Number of chromatograms cached in memory.
    size_t getNrChromatograms() const override;

    /**
      @brief Native id of one cached chromatogram.

      @param[in] id Chromatogram index in @c [0, getNrChromatograms()).
                    Out-of-range access is a programming error (checked in
                    debug builds).
      @return Native-id string captured at construction time.
    */
    std::string getChromatogramNativeID(int id) const override;

private:

    std::vector< OpenSwath::SpectrumPtr > spectra_;       ///< Spectrum payloads captured at construction; shared with light clones.
    std::vector< OpenSwath::SpectrumMeta > spectra_meta_; ///< Parallel metadata cache (native id, RT, MS level) keyed on the same index.

    std::vector< OpenSwath::ChromatogramPtr > chromatograms_; ///< Chromatogram payloads captured at construction; shared with light clones.
    std::vector< std::string > chromatogram_ids_;             ///< Parallel chromatogram native-id cache keyed on the same index.

  };

} //end namespace OpenMS

