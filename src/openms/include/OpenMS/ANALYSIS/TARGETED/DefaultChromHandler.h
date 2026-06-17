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
#include <OpenMS/ANALYSIS/TARGETED/MRMChromHandler.h>
#include <OpenMS/ANALYSIS/TARGETED/DIAChromHandler.h>
#include <memory>

namespace OpenMS
{
  /**
    @brief Dispatching chromatogram handler that picks SRM/MRM or DIA at runtime based on the input maps.

    Wraps an @ref MRMChromHandler and a @ref DIAChromHandler and forwards each call to whichever
    is appropriate for the supplied SWATH maps. The mode is decided on every call by inspecting
    @p swath_maps:

      - If every map is chromatogram-only (neither flagged as MS1 nor carrying any spectra),
        the call is routed to the SRM/MRM handler.
      - Otherwise (any map flagged @c ms1 or carrying at least one spectrum), the call is routed
        to the DIA handler.

    SRM/MRM error handling: when the SRM/MRM path raises @ref Exception::IllegalArgument
    (typically "no chromatograms mapped to transitions" — common for chromatogram-only inputs
    whose precursor/product m/z don't match the supplied transition list), the exception is
    logged as a warning and an empty vector is returned. All other exceptions, and any exception
    raised on the DIA path, propagate to the caller.

    Construction is cheap (default-constructs both inner handlers up front); no input data is
    inspected until the first call.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI DefaultChromHandler : public IChromatogramHandler
  {
  public:
    /// Default constructor; eagerly constructs the inner SRM/MRM and DIA handlers
    DefaultChromHandler();

    /// Destructor
    ~DefaultChromHandler() override;

    /**
      @brief Collect iRT chromatograms by delegating to the SRM/MRM or DIA handler.

      Routing rule: SRM/MRM if every @p swath_maps entry is chromatogram-only (neither @c ms1
      nor any spectra), DIA otherwise. On the SRM/MRM path, @ref Exception::IllegalArgument
      ("no chromatograms mapped") is logged and converted to an empty result; all other
      exceptions propagate unchanged.

      @param[in] swath_maps        Input maps (mixed-mode allowed; the per-map check decides routing).
      @param[in] irt_transitions   iRT transition list used for mapping.
      @param[in] mrm_mapping_param Parameters forwarded to the @ref MRMMapping step.
      @param[in] cp                Chromatogram-extraction parameters (windows, RT range, etc.).
      @param[in] trafo             Optional RT transformation applied to iRT transitions.
      @param[in] pasef             If true, enable PASEF/ion-mobility-aware mapping in the delegate.
      @param[in] load_into_memory  If true, materialise on-disk data before extraction.
      @return Mapped iRT chromatograms (empty on the recoverable SRM/MRM "no mapping" case).
      @throws Exception::BaseException Propagated for non-recoverable SRM/MRM delegate failures and all DIA-path exceptions.
      @note On the SRM/MRM path, @ref Exception::IllegalArgument ("no chromatograms mapped") is logged and converted to an empty result.
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
      @brief Extract and map analytical transition chromatograms via the SRM/MRM or DIA handler.

      Routing rule and error handling are identical to collectIrtChromatogramsForIrt(): SRM/MRM
      for chromatogram-only inputs (where @ref Exception::IllegalArgument is downgraded to an
      empty result), DIA otherwise (all exceptions propagate).

      @param[in] swath_maps        Input maps (mixed-mode allowed; the per-map check decides routing).
      @param[in] transition_exp    Analytical transition list to map against.
      @param[in] cp                Chromatogram-extraction parameters.
      @param[in] mrm_mapping_param Parameters forwarded to the @ref MRMMapping step.
      @return Mapped, filtered chromatograms ready for downstream scoring.
      @throws Exception::BaseException Propagated for non-recoverable SRM/MRM delegate failures and all DIA-path exceptions.
      @note On the SRM/MRM path, @ref Exception::IllegalArgument ("no chromatograms mapped") is logged and converted to an empty result.
    */
    std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const ChromExtractParams & cp,
      const Param & mrm_mapping_param) override;

  private:
    /// Inner SRM/MRM handler; selected when all input maps are chromatogram-only
    std::unique_ptr<MRMChromHandler> mrm_;

    /// Inner DIA handler; selected when at least one input map carries MS1 or spectra
    std::unique_ptr<DIAChromHandler> dia_;
  };

} // namespace OpenMS
