// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

namespace OpenMS
{
  /**
    @brief Scoring of an spectrum given library intensities of a transition group.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    See also class DIAScoring.

    Simulate theoretical spectrum from library intensities of transition group
    and compute manhattan distance and dotprod score between spectrum intensities
    and simulated spectrum.
  */

  class OPENMS_DLLAPI DiaPrescore :
    public DefaultParamHandler
  {
    double dia_extract_window_; //done
    int nr_isotopes_;
    int nr_charges_;
public:
    /**
      @brief Cached transition-group-specific theoretical spectrum used by the DIA prescore.

      Stores precomputed theoretical m/z and intensity arrays plus rescaling
      bounds for repeated scoring of the same transition group.
    */
    struct OPENMS_DLLAPI TransitionGroupTheoreticalSpectrumCache
    {
      /// Theoretical fragment/isotope m/z positions used for window integration.
      std::vector<double> mz_theor;
      /// Non-negative normalized theoretical intensities used for Manhattan-style comparison.
      std::vector<double> int_theor;
      /// Signed normalized theoretical intensities including negatively weighted pre-isotope peaks.
      std::vector<double> int_theor_neg;
      /// Lower bound used to rescale the signed dot-product score.
      double neg_val = 0.0;
      /// Upper bound used to rescale the signed dot-product score.
      double pos_val = 0.0;

      /**
        @brief Reset all cached spectrum data.

        Clears the cached arrays and resets scaling values to zero.
      */
      void clear();
    };

    DiaPrescore();

    DiaPrescore(double dia_extract_window, int nr_isotopes = 4, int nr_charges = 4);

    void defineDefaults();

    void updateMembers_() override;

    /**
      @brief Score one or more observed spectra for a transition group.

      Builds the transition-group-specific theoretical spectrum cache from the
      library intensities and compares it against the observed spectrum
      sequence. The input may contain multiple adjacent MS2 spectra / frames
      around the apex or a single pre-merged spectrum, depending on the
      upstream merge mode.
    */
    void score(const SpectrumSequence& spec,
               const std::vector<OpenSwath::LightTransition>& lt,
               const RangeMobility& im_range,
               double& dotprod,
               double& manhattan) const;

    /**
      @brief Build the theoretical DIA isotope spectrum cache for repeated scoring.

      The transition-group-specific theoretical spectrum cache depends only on
      the transition group and DiaPrescore parameters, not on the observed
      spectrum.

      @param[in] lt Library transitions for the transition group
      @param[out] theoretical_spectrum_cache Precomputed theoretical spectrum cache
    */
    void buildTheoreticalSpectrum(const std::vector<OpenSwath::LightTransition>& lt,
                                  TransitionGroupTheoreticalSpectrumCache& theoretical_spectrum_cache) const;

    /**
      @brief Score an observed spectrum sequence against a precomputed theoretical DIA spectrum.

      The observed input may be a single merged spectrum or multiple spectra /
      frames around the feature apex, depending on the configured spectrum
      addition mode.

      @param[in] spec Observed spectrum sequence around the feature apex
      @param[in] theoretical_spectrum_cache Precomputed theoretical spectrum cache
      @param[in] dia_extract_window DIA extraction window in Th
      @param[in] im_range Ion mobility extraction range
      @param[out] dotprod Dot product score
      @param[out] manhattan Manhattan distance score
    */
    static void scorePrepared(const SpectrumSequence& spec,
                              const TransitionGroupTheoreticalSpectrumCache& theoretical_spectrum_cache,
                              double dia_extract_window,
                              const RangeMobility& im_range,
                              double& dotprod,
                              double& manhattan);

    /**
      @brief Compute Manhattan and dot-product scores for all spectra accessible through the SpectrumAccessPtr.

      Scores all transition groups in the LightTargetedExperiment against the
      spectra provided by @p swath_ptr and writes the per-spectrum score vectors
      to @p ivw.
    */
    void operator()(const OpenSwath::SpectrumAccessPtr& swath_ptr,
                    OpenSwath::LightTargetedExperiment& transition_exp_used, const RangeMobility& range_im,
                    OpenSwath::IDataFrameWriter* ivw) const;
  };


}
