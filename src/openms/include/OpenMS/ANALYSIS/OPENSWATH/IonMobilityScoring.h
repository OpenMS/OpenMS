// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  struct RangeMobility;
  struct RangeMZ;

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

    struct MobilityPeak
    {
      double im;
      double intensity;
      MobilityPeak ();
      MobilityPeak (double im_, double int_) :
        im(im_),
        intensity(int_)
      {}
    };
    typedef std::vector< MobilityPeak > IonMobilogram;

  public:

    /// Constructor
    IonMobilityScoring();

    /// Destructor
    ~IonMobilityScoring();

    /**
      @brief Performs scoring of the ion mobility dimension in MS2

      Populates additional scores in the @p scores object

      @param spectra Sequence of segments of the DIA MS2 spectrum found at (and around) the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_target Ion Mobility extraction target
      @param im_range Ion Mobility extraction range
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)
    */
    static void driftScoring(const SpectrumSequence& spectra,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_target,
                             RangeMobility im_range,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1

      Populates additional scores in the @p scores object

      @param spectra vector containing the DIA MS1 spectra found at (or around) the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param im_range Ion Mobility extraction range
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extra extraction to use for drift time (in percent)
    */
    static void driftScoringMS1(const SpectrumSequence& spectra,
                                const std::vector<TransitionType> & transitions,
                                OpenSwath_Scores & scores,
                                const double drift_target,
                                RangeMobility im_range,
                                const double dia_extraction_window_,
                                const bool dia_extraction_ppm_,
                                const bool use_spline,
                                const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1 and MS2 (contrast)

      Populates additional scores in the @p scores object

      @param spectra Vector of the DIA MS2 spectrum found in SpectrumSequence object (can contain 1 or multiple spectra centered around peak apex)
      @param ms1spectrum The DIA MS1 spectrum found in SpectrumSequence object (can contain 1 or multiple spectra centered around peak apex)
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param im_range the ion mobility range
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param drift_extra Extra extraction to use for drift time (in percent)
    */
    static void driftScoringMS1Contrast(const SpectrumSequence& spectra, const SpectrumSequence& ms1spectrum,
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        RangeMobility im_range,
                                        const double dia_extraction_window_,
                                        const bool dia_extraction_ppm_,
                                        const double drift_extra);

    /**
     * @brief computes ion mobilogram to be used in scoring based on mz_range and im_range.
     * Also integrates intensity in the resulting ion mobility mobilogram in mz_range and im_range across all the entire SpectrumSequence.
     * @note If there is no signal, mz will be set to -1 and intensity to 0
     * @param[in] spectra Raw data in a spectrumSequence object (can contain 1 or multiple spectra centered around peak apex)
     * @param[in] mz_range the range across mz to extract
     * @param[in] im_range the range across im to extract
     * @param[out] im computed weighted average ion mobility
     * @param[out] intensity intensity computed intensity
     * @param[out] res outputted ion mobilogram
     * @param[in] eps minimum distance to allow for two seperate points
     */
    static void computeIonMobilogram(const SpectrumSequence& spectra,
                              const RangeMZ & mz_range,
                              const RangeMobility & im_range,
                              double & im,
                              double & intensity,
                              IonMobilogram& res,
                              double eps);


  private:
    /**
     * @brief helper function to computeIonMobilogram. Discretizes ion mobility values into a grid.
    **/
    static std::vector<double> computeGrid_(const std::vector< IonMobilogram >& mobilograms, double eps);


    /*
     @brief Extracts ion mobility values projected onto a grid

     For a given ion mobility profile and a grid, compute an ion mobilogram
     across the grid for each ion mobility data point. Returns two data arrays
     for the ion mobilogram: intensity (y) and ion mobility (x). Zero values are
     inserted if no data point was found for a given grid value.

     @param profile The ion mobility data
     @param im_grid The grid to be used
     @param al_int_values The intensity vector (y)
     @param al_im_values The ion mobility vector (x)
     @param eps Epsilon used for computing the ion mobility grid
     @param max_peak_idx The grid position of the maximum
    */
    static void alignToGrid_(const IonMobilogram& profile,
                 const std::vector<double>& im_grid,
                 std::vector< double >& al_int_values,
                 std::vector< double >& al_im_values,
                 double eps,
                 Size & max_peak_idx);

  };
}

