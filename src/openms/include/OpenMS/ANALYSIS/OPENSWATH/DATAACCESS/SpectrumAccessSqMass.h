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

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <memory>

namespace OpenMS
{
  /**
    @brief @ref OpenSwath::ISpectrumAccess implementation backed by an sqMass (SQLite) spectrum store.

    Wraps an @ref OpenMS::Internal::MzMLSqliteHandler to expose its spectra through the
    OpenSWATH spectrum-access interface. Optionally restricts the visible set to a subset of
    spectra via a vector of indices, so the same on-disk store can back several views
    (e.g. one per SWATH window, or MS1-only / MS2-only).

    @section SpectrumAccessSqMass_subsetting Subsetting

    Three constructors are provided:
      - @c (handler) — expose every spectrum in @p handler.
      - @c (handler, indices) — expose only the spectra at the given indices.
      - @c (parent, indices) — narrow an existing @ref SpectrumAccessSqMass further; @p indices
        are interpreted against the parent's already-filtered view. Out-of-range entries
        throw @ref Exception::IllegalArgument.

    All @c getXxxById methods honour the subset: id @c i refers to the @c i-th element of the
    visible view, not the @c i-th row of the underlying SQLite table.

    @section SpectrumAccessSqMass_perf Performance and concurrency

    Per-spectrum access (@ref getSpectrumById / @ref getSpectrumMetaById) is implemented but
    each call performs one @ref OpenMS::Internal::MzMLSqliteHandler::readSpectra invocation
    and incurs the per-query SQLite overhead. Callers that need many spectra should prefer
    @ref getAllSpectra, which issues a single batched read.

    The class is read-only and SQLite supports parallel reads from independent connections,
    so cloning via @ref lightClone gives a clone safe to use from another thread.

    @section SpectrumAccessSqMass_chromos Chromatogram support

    Spectrum-side operations are fully supported. On the chromatogram side, only the
    count @ref getNrChromatograms is implemented; both @ref getChromatogramById and
    @ref getChromatogramNativeID throw @ref Exception::NotImplemented. Use a different
    accessor (e.g. read via @ref OpenMS::Internal::MzMLSqliteHandler directly) for
    chromatogram payloads.

    @section SpectrumAccessSqMass_example Sample usage

    @code
      // Obtain swath_map with boundaries first
      std::vector<int> indices = sql_mass_reader.readSpectraForWindow(swath_map);
      OpenMS::Internal::MzMLSqliteHandler handler(file);
      OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
      swath_maps[k].sptr = sptr;
    @endcode

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SpectrumAccessSqMass :
    public OpenSwath::ISpectrumAccess
  {

public:
    typedef OpenMS::MSSpectrum MSSpectrumType;
    typedef OpenMS::MSChromatogram MSChromatogramType;

    /**
      @brief Construct from an sqMass handler exposing every spectrum it contains.

      @param[in] handler Read-only handler to the underlying sqMass file.
    */
    SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler);

    /**
      @brief Construct from an sqMass handler exposing only the spectra at the given indices.

      @param[in] handler Read-only handler to the underlying sqMass file.
      @param[in] indices Absolute spectrum indices to expose; an empty vector falls back to "all spectra" semantics.
    */
    SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int> & indices);

    /**
      @brief Construct a narrower view of an existing accessor.

      If @p sp already has a subset configured, @p indices are interpreted against that subset
      (i.e. @p indices[k] is a position within @p sp's already-visible view, not within the
      raw sqMass file).

      @param[in] sp      Parent accessor.
      @param[in] indices Indices into @p sp's already-filtered view; empty inherits @p sp's subset unchanged.
      @throws Exception::IllegalArgument If any entry of @p indices is out of range of @p sp's subset.
    */
    SpectrumAccessSqMass(const SpectrumAccessSqMass& sp, const std::vector<int>& indices);

    /// Destructor
    ~SpectrumAccessSqMass() override;

    /**
      @brief Copy constructor.

      Duplicates the @ref OpenMS::Internal::MzMLSqliteHandler value and the optional
      subset-index vector; no shared mutable state remains between the copy and the
      original.

      @param[in] rhs Source accessor to copy.
    */
    SpectrumAccessSqMass(const SpectrumAccessSqMass & rhs);

    /**
      @brief Return a copy of this accessor as a new @ref OpenSwath::ISpectrumAccess pointer.

      The clone uses an independent @ref OpenMS::Internal::MzMLSqliteHandler value, so it is
      safe to use from another thread (SQLite supports parallel reads from independent
      connections).

      @return Shared pointer to a new accessor with the same visible view as this one.
    */
    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    /**
      @brief Load one spectrum by index in the visible view.

      Issues a single @ref OpenMS::Internal::MzMLSqliteHandler::readSpectra call per
      invocation, so the SQLite per-query overhead applies to each call.

      @param[in] id Spectrum index in the visible view (0-based; honours the optional subset).
      @return Pointer to the loaded spectrum.
      @warning Prefer @ref getAllSpectra when many spectra are needed; per-call cost is dominated by the SQLite query overhead.
    */
    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    /**
      @brief Load metadata (native id, RT, MS level) for one spectrum.

      Same per-call cost trade-off as @ref getSpectrumById.

      @param[in] id Spectrum index in the visible view (0-based).
      @return Metadata record for the requested spectrum.
    */
    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    /**
      @brief Batch-load every spectrum in the visible view into memory.

      Preferred over repeated @ref getSpectrumById calls — issues a single batched read.

      @param[out] spectra      Receives one @c SpectrumPtr per visible spectrum.
      @param[out] spectra_meta Receives the parallel metadata records.
    */
    void getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta) const;

    /**
      @brief Indices into the visible view of spectra whose RT lies within an absolute window.

      @param[in] RT      Centre of the RT window (same units as the underlying spectra).
      @param[in] deltaRT Half-width of the RT window; matches spectra in @f$[RT - \mathrm{deltaRT}, RT + \mathrm{deltaRT}]@f$.
      @return Indices into the visible view (already remapped through the subset, if any).
      @note Aborts via @c OPENMS_PRECONDITION if @p deltaRT is negative.
    */
    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    /**
      @brief Number of spectra in the visible view.

      @return Subset size when a subset is configured, otherwise the total spectrum count of the underlying sqMass file.
    */
    size_t getNrSpectra() const override;

    /**
      @brief Not supported on sqMass.

      @param[in] id Ignored.
      @return Never returns normally.
      @throws Exception::NotImplemented Always — sqMass-backed access does not implement per-chromatogram payload retrieval. Use a different accessor for chromatogram payloads.
    */
    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    /**
      @brief Number of chromatograms in the underlying sqMass file.

      The count is supported even though per-chromatogram payload / native-id access is not.

      @return Total chromatogram count.
    */
    size_t getNrChromatograms() const override;

    /**
      @brief Not supported on sqMass.

      @param[in] id Ignored.
      @return Never returns normally.
      @throws Exception::NotImplemented Always — sqMass-backed access does not expose chromatogram native ids.
    */
    std::string getChromatogramNativeID(int id) const override;

private:

    /// Handler for the underlying sqMass file (copied by value into each accessor / clone)
    OpenMS::Internal::MzMLSqliteHandler handler_;
    /// Optional subset of absolute spectrum indices; empty means "all spectra in the file are visible"
    std::vector<int> sidx_;
  };
} //end namespace OpenMS



