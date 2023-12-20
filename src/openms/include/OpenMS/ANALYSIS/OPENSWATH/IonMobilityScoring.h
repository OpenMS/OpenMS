// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// data access
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <vector>

namespace OpenMS
{

  /** @brief A class that calls the ion mobility scoring routines
   *
   * Use this class to invoke the individual OpenSWATH ion mobility scoring
   * routines. These scores use the ion mobilograms from individual peptides in
   * one (or more) frames to compute additional scores.
   *
   * - driftScoring() performs scoring on fragment ion mobilograms extracted from a DIA frame
   * - driftScoringMS1() performs scoring on precursor ion mobilograms extracted from a MS1 frame
   * - driftScoringMS1Contrast() performs cross correlation (contrast) scoring between precursor and fragment ion mobilograms
   *
  */
  class OPENMS_DLLAPI IonMobilityScoring
  {
    typedef OpenSwath::LightCompound CompoundType;
    typedef OpenSwath::LightTransition TransitionType;

  public:

    /// Constructor
    IonMobilityScoring();

    /// Destructor
    ~IonMobilityScoring();

    /**
      @brief Performs scoring of the ion mobility dimension in MS2

      @param spectrum The DIA MS2 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoring(const OpenSwath::SpectrumPtr& spectrum,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_lower,
                             const double drift_upper,
                             const double drift_target,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1

      @param spectrum The DIA MS1 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extra extraction to use for drift time (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoringMS1(const OpenSwath::SpectrumPtr& spectrum,
                                const std::vector<TransitionType> & transitions,
                                OpenSwath_Scores & scores,
                                const double drift_lower,
                                const double drift_upper,
                                const double drift_target,
                                const double dia_extract_window_,
                                const bool dia_extraction_ppm_,
                                const bool use_spline,
                                const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1 and MS2 (contrast)

      @param spectrum The DIA MS2 spectrum found at the peak apex
      @param ms1spectrum The DIA MS1 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extra extraction to use for drift time (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoringMS1Contrast(const OpenSwath::SpectrumPtr& spectrum, const OpenSwath::SpectrumPtr& ms1spectrum,
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_lower,
                                        const double drift_upper,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const double drift_extra);
  };
}

