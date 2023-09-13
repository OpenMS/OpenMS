// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <vector>

namespace OpenMS
{

  /**
    @brief This class implements a fast peak-picking algorithm best suited for
    high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
    no overlapping and therefore allow for a clear separation. Furthermore, ion
    signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in raw data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

    @htmlinclude OpenMS_PeakPickerHiRes.parameters

    @note The peaks must be sorted according to ascending m/z!

    @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerMaxima 
  {
public:

    /// Constructor
    PeakPickerMaxima(double signal_to_noise, double spacing_difference = 1.5,
        double spacing_difference_gap = 4.0, double sn_window_length = 200,
        unsigned missing = 2);

    /// Destructor
    virtual ~PeakPickerMaxima() {}

    /**
      @brief The PeakCandidate describes the output of the peak picker

      It contains the m/z and intensity value of the peak candidate.

      It also contains the original index in the m/z axis where the peak was
      found as well as an estimate of its right and left boundary. 
    */
    struct OPENMS_DLLAPI PeakCandidate 
    {
      /// index of the peak apex (relative to the input data) 
      int pos;
      /// index of the left boundary (relative to the input data) 
      int left_boundary;
      /// index of the right boundary (relative to the input data) 
      int right_boundary;
      /// m/z value of the peak apex
      double mz_max;
      /// intensity value of the peak apex
      double int_max;
    };

    /**
      @brief Will find local maxima in raw data

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates
      @param check_spacings check spacing constraints? (recommended settings: yes for spectra, no for chromatograms)

      @note This function will directly report peak apices with right and left
      boundaries but will not use any fitting to estimate the true m/z and
      intensity of the peak. Note that the mz_max and int_max fields will be
      empty in the result (set to -1).

    */
    void findMaxima(const std::vector<double>& mz_array,
        const std::vector<double>& int_array,
        std::vector<PeakCandidate>& pc,
        bool check_spacings = true) const;

    /**
      @brief Will pick peaks in a spectrum

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates
      @param check_spacings check spacing constraints? (recommended settings: yes for spectra, no for chromatograms)

      @note This function will first find maxima in the intensity domain and
      then use a spline function to estimate the best m/z and intensity for
      each peak candidate.
    */
    void pick(std::vector<double>& mz_array, 
        std::vector<double>& int_array, 
        std::vector<PeakCandidate>& pc,
        bool check_spacings = true);

protected:

    // signal-to-noise parameter
    double signal_to_noise_;

    // signal-to-noise window length 
    double sn_window_length_;

    // maximal spacing difference defining a missing data point
    double spacing_difference_;

    // maximal spacing difference defining a large gap
    double spacing_difference_gap_;

    // maximum number of missing points
    unsigned missing_;

  }; // end PeakPickerMaxima

} // namespace OpenMS

