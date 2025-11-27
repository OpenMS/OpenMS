// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/KERNEL/Mobilogram.h>

// Kernel classes
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h>

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
    typedef MRMTransitionGroup< MSChromatogram, TransitionType> MRMTransitionGroupType;

  public:

    /// Constructor
    IonMobilityScoring();

    /// Destructor
    ~IonMobilityScoring();

    /**
      @brief Performs scoring of the ion mobility dimension in MS2

      Populates additional scores in the @p scores object

      If @p apply_im_peak_picking is set to true, peak picking is performed on the Savitzky-Golay smoothed ion mobilogram. This is useful for minimizing interference from co-eluting analytes in the ion mobility dimension (IM) that fall within the current extraction window. This process improves the specificity of analyte detection by dynamically adjusting the IM extraction window to extract only over the IM elution of the highest intensity species. If multiple peaks are present in the IM dimension, lower intensity peaks get discarded.

      @param spectra Sequence of segments of the DIA MS2 spectrum found at (and around) the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_target Ion Mobility extraction target
      @param im_range Ion Mobility extraction range
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param drift_extra Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)
      @param apply_im_peak_picking Apply peak picking on the ion mobilogram
    */
    static void driftScoring(const SpectrumSequence& spectra,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_target,
                             RangeMobility im_range,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const double drift_extra,
                             const bool apply_im_peak_picking);

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
      @param drift_extra Extra extraction to use for drift time (in percent)
    */
    static void driftScoringMS1(const SpectrumSequence& spectra,
                                const std::vector<TransitionType> & transitions,
                                OpenSwath_Scores & scores,
                                const double drift_target,
                                RangeMobility im_range,
                                const double dia_extraction_window_,
                                const bool dia_extraction_ppm_,
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
      @brief Performs scoring of the ion mobility dimension for identification transitions against detection transitions

      @param spectra Vector of the DIA MS2 spectrum found in SpectrumSequence object (can contain 1 or multiple spectra centered around peak apex)
      @param transitions The transitions used for scoring
      @param transition_group_detection The detection transition group
      @param scores The output scores
      @param drift_target Ion Mobility extraction target
      @param im_range Ion Mobility extraction range
      @param dia_extract_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param drift_extra Extra extraction to use for drift time (in percent)
      @param apply_im_peak_picking Apply peak pickng on the ion mobilogram
    */
    static void driftIdScoring(const SpectrumSequence& spectra,
                                const std::vector<TransitionType> & transitions,
                                MRMTransitionGroupType& transition_group_detection,
                                OpenSwath_Scores & scores,
                                const double drift_target,
                                RangeMobility im_range,
                                const double dia_extract_window_,
                                const bool dia_extraction_ppm_,
                                const double drift_extra,
                                const bool apply_im_peak_picking);

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
                              Mobilogram & res,
                              double eps);


  private:
    /**
     * @brief helper function to computeIonMobilogram. Discretizes ion mobility values into a grid.
    **/
    static std::vector<double> computeGrid_(const std::vector< Mobilogram >& mobilograms, double eps);


    /*
     @brief Extracts ion mobility values projected onto a grid

     For a given ion mobility profile and a grid, compute an ion mobilogram
     across the grid for each ion mobility data point. Returns two data arrays
     for the ion mobilogram: intensity (y) and ion mobility (x). Zero values are
     inserted if no data point was found for a given grid value.

     @param profile The ion mobility data
     @param im_grid The grid to be used
     @param eps Epsilon used for computing the ion mobility grid
     @param max_peak_idx The grid position of the maximum
    */
    static void alignToGrid_(const Mobilogram& profile,
                 const std::vector<double>& im_grid,
                 Mobilogram & aligned_profile,
                 double eps,
                 Size & max_peak_idx);

    /*
    @brief Extracts intensity values from a vector of Mobilogram objects

    This function takes a vector of Mobilogram objects and extracts the intensity
    values from each Mobilogram, storing them in a 2D vector of doubles. The
    resulting vector of intensity values is stored in the provided output parameter.

    @param[in] mobilograms A const reference to a vector of Mobilogram objects
                            from which to extract intensity values.
    @param[out] int_values A reference to a vector of vector of doubles where
                            the extracted intensity values will be stored. This
                            vector will be cleared and resized as necessary.
    */
    static void extractIntensities(const std::vector< Mobilogram >& mobilograms,
                                   std::vector<std::vector<double>>& int_values);

    /**
     * @brief Extracts intensity values from a single Mobilogram object
     *
     * This function takes a single Mobilogram object and extracts the intensity
     * values, returning them as a vector of doubles.
     *
     * @param mobilogram [in] A const reference to a Mobilogram object from which to extract intensity values.
     * @return A vector of doubles containing the extracted intensity values.
     */
    static std::vector<double> extractIntensities(const Mobilogram& mobilogram);

  };

  /*
  @brief Helper function to sum up aligned mobilograms

  This function takes a vector of aligned mobilograms and sums them up to create a single Mobilogram object. The resulting Mobilogram object will contain the sum of the intensity values from all the input mobilograms.

  @note the input vector of Mobilograms all need to have the same size, otherwise the function will throw a pre-condition exception.

  @param[in] aligned_mobilograms A vector of aligned mobilograms
  @return  A Mobilogram object that is the sum of the input mobilograms
  */
  OpenMS::Mobilogram sumAlignedMobilograms(const std::vector<OpenMS::Mobilogram>& aligned_mobilograms);

  } // namespace OpenMS
