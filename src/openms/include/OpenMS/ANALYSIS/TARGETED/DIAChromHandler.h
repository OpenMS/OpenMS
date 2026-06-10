// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

namespace OpenMS
{
  /**
    @brief IChromatogramHandler implementation for DIA / SWATH / PRM /
    PASEF-DIA acquisitions, where chromatograms are @em extracted from
    full-scan MS2 data rather than being directly acquired as SRM/MRM
    traces.

    DIAChromHandler is the DIA-shaped sibling of @ref MRMChromHandler.
    The two implement the same @ref IChromatogramHandler interface but
    differ in where their chromatograms come from:

      - MRMChromHandler @em selects pre-acquired SRM/MRM chromatograms
        that were stored alongside the spectra and maps them to the
        target transitions.
      - DIAChromHandler @em extracts new chromatograms from full-scan
        MS2 (SWATH) data using @ref ChromatogramExtractor. Mapping
        back to transitions is handled implicitly by
        @ref ChromatogramExtractor::return_chromatogram, which assigns
        the transition native IDs to the extracted traces; the
        @p mrm_mapping_param argument of either method is therefore
        silently ignored.

    @ref DefaultChromHandler is a thin dispatcher that owns one of each
    and forwards to whichever variant matches the run's acquisition mode.

    @section DIAChromHandler_pasef PASEF support

    Both override methods accept a @p pasef flag. When set, the handler
    performs an extra precursor-matching pass that requires every
    transition to carry a valid precursor IM value (via
    @c LightTransition::getPrecursorIM); transitions are then routed to
    the SWATH window whose IM range covers their precursor (see
    @c OpenSwathHelper::selectSwathTransitionsPasef). Without the flag
    the standard m/z-only window selection is used.

    @section DIAChromHandler_memory In-memory acceleration

    When @p load_into_memory is true the swath maps are wrapped in
    @ref SpectrumAccessOpenMSInMemory before extraction, eliminating
    repeated disk I/O at the cost of holding the spectra in RAM for the
    duration of the call. Recommended for iRT calibration on small
    iRT-transition lists.

    @ingroup TargetedQuantitation
    @see IChromatogramHandler
    @see MRMChromHandler
    @see DefaultChromHandler
    @see ChromatogramExtractor
  */
  class OPENMS_DLLAPI DIAChromHandler : public IChromatogramHandler
  {
  public:
    DIAChromHandler() = default;
    ~DIAChromHandler() override = default;

    /**
      @brief Extract chromatograms for the iRT calibration transitions
      from each non-MS1 SWATH map.

      Per-map workflow:
        - MS1 maps are skipped.
        - In PASEF mode (@p pasef true), every transition is required to
          carry a valid precursor IM and @c OpenSwathHelper::selectSwathTransitionsPasef
          picks one map per transition by matching precursor IM/m/z.
        - Otherwise, @c OpenSwathHelper::selectSwathTransitions picks
          the transitions whose precursor m/z falls into the map's
          [lower, upper] window.
        - If no transitions match but the map is chromatogram-only
          (zero spectra, non-zero chromatograms), the full
          @p irt_transitions list is used as fallback so that native-ID
          matching against pre-acquired chromatograms can still happen.
        - The RT trafo is inverted and applied to the per-coordinate
          start/end times when @p cp.rt_extraction_window >= 0.
        - Chromatograms are extracted via
          @ref ChromatogramExtractor::extractChromatograms and
          materialised by @ref ChromatogramExtractor::return_chromatogram,
          which also assigns the transition native IDs (so a separate
          @ref MRMMapping pass is not needed).
        - Per-chromatogram TIC is checked; chromatograms with TIC == 0
          (empty extraction window) are dropped and logged.

      The outer loop over swath maps is OpenMP-parallel; the result
      vector is updated under a critical section.

      @param[in] swath_maps         List of SWATH maps to extract from
                                    (MS1 maps are skipped).
      @param[in] irt_transitions    iRT calibration transitions
                                    (precursor and product info).
      @param[in] mrm_mapping_param  Currently unused (no MRMMapping pass
                                    is run); kept only to satisfy the
                                    @ref IChromatogramHandler interface.
      @param[in] cp                 Chromatogram extraction parameters
                                    (window widths, RT range, etc.).
      @param[in] trafo              RT transformation; default is
                                    identity. Inverted before being
                                    applied to extraction coordinates.
      @param[in] pasef              If true, route transitions to maps
                                    by precursor IM in addition to m/z;
                                    every transition must have a valid
                                    @c getPrecursorIM value.
      @param[in] load_into_memory   If true, wrap each map in
                                    @ref SpectrumAccessOpenMSInMemory
                                    before extracting.
      @return Non-empty (TIC > 0) iRT chromatograms in extraction order.
      @throws Exception::InvalidParameter when @p pasef is true and any
              iRT transition has @c getPrecursorIM == -1.
    */
    std::vector<MSChromatogram> collectIrtChromatogramsForIrt(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & irt_transitions,
      const Param & mrm_mapping_param,
      const ChromExtractParams & cp,
      const TransformationDescription& trafo = TransformationDescription(),
      bool pasef = false,
      bool load_into_memory = false) override;

    /**
      @brief Extract chromatograms for the analytical transitions from each
      non-MS1 SWATH map.

      For every non-MS1 swath map with a non-null spectrum-access pointer,
      the entire @p transition_exp is handed to
      @ref ChromatogramExtractor::prepare_coordinates (no per-window
      precursor-m/z filtering is performed here), the chromatograms are
      extracted via @ref ChromatogramExtractor::extractChromatograms, then
      @ref ChromatogramExtractor::return_chromatogram populates them and
      assigns native IDs from the transition list directly. Extraction
      failures (any @c std::exception thrown by the extractor) on a single
      map are caught and logged at DEBUG, then the loop continues with the
      next map.

      Mapping to transitions is therefore done implicitly via the native
      IDs assigned by @ref ChromatogramExtractor — no separate
      @ref MRMMapping pass is run, and the result is the concatenation of
      everything extracted (no per-map filtering down to a "mapped subset"
      either).

      @param[in] swath_maps        List of SWATH maps to extract from. MS1
                                   maps and maps with a null
                                   @c sptr are skipped.
      @param[in] transition_exp    Analytical transition list (passed in
                                   full to every non-skipped map).
      @param[in] cp                Chromatogram extraction parameters
                                   (window widths, RT range, etc.).
      @param[in] mrm_mapping_param Currently unused (no MRMMapping pass is
                                   run in the DIA handler); kept only to
                                   satisfy the @ref IChromatogramHandler
                                   interface.
      @return All extracted chromatograms with their native IDs set by
              @ref ChromatogramExtractor::return_chromatogram.
    */
    std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const ChromExtractParams & cp,
      const Param & mrm_mapping_param) override;
  };

} // namespace OpenMS
