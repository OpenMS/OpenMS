// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <map>

namespace OpenMS
{
  /**
    @brief Implementation of the PScore PSM scoring algorithm
  */

struct OPENMS_DLLAPI PScore
{
  /** @brief calculate local (windowed) peak ranks.
   
    The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
    The result can be used to efficiently filter spectra for top 1..n peaks in mass windows 
   
    @note ranks are zero based (highest intensity peak in window has rank 0)
   
    @param mz m/z positions of the peaks
    @param intensities of the peaks
    @param mz_window window in Thomson centered at each peak
   */  
  static std::vector<Size> calculateIntensityRankInMZWindow(const std::vector<double>& mz, const std::vector<double>& intensities, double mz_window);

  /** @brief precalculated, windowed peak ranks for a whole experiment. 
   
    The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
   
    -# Each spectrum is subdivided into windows of size \p mz_window.
    -# For each window, peak ranks are assigned using calculateIntensityRankInMZWindow().
    -# A rank map is returned
   
    @note ranks are zero based (top element has rank 0)
   
    @param peak_map Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra.
    @param mz_window window in Thomson centered at each peak
   */
  static std::vector<std::vector<Size> > calculateRankMap(const PeakMap& peak_map, double mz_window = 100);

  /** @brief Calculates spectra for peak level between min_level to max_level and stores them in the map
    A spectrum of peak level n retains the (n+1) top intensity peaks in a sliding mz_window centered at each peak.
   
    @note levels are zero based (level 0 has only the top intensity peaks for each window, level 1 the top and second most intensive one)
    @note min and max level are taken from the Andromeda publication but are similar to the AScore publication
   
   */ 
  static std::map<Size, PeakSpectrum > calculatePeakLevelSpectra(const PeakSpectrum& spec, const std::vector<Size>& ranks, Size min_level = 1, Size max_level = 9);

  /** @brief Computes the PScore for a vector of theoretical spectra
   
    Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
    The best score obtained by scoring all those theoretical spectra against the experimental ones is returned.
   
    @param fragment_mass_tolerance mass tolerance for matching peaks
    @param fragment_mass_tolerance_unit_ppm whether Thomson or ppm is used
    @param peak_level_spectra spectra for different peak levels (=filtered by maximum rank).
    @param theo_spectra theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
    @param mz_window window in Thomson centered at each peak
   */ 
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, PeakSpectrum>& peak_level_spectra, const std::vector<PeakSpectrum>& theo_spectra, double mz_window = 100.0);

  /** @brief Computes the PScore for a single theoretical spectrum
   
    @param fragment_mass_tolerance mass tolerance for matching peaks
    @param fragment_mass_tolerance_unit_ppm whether Thomson or ppm is used
    @param peak_level_spectra spectra for different peak levels (=filtered by maximum rank).
    @param theo_spectrum Theoretical spectrum as obtained e.g. from TheoreticalSpectrumGenerator
    @param mz_window window in Thomson centered at each peak
   */ 
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, PeakSpectrum>& peak_level_spectra, const PeakSpectrum& theo_spectrum, double mz_window = 100.0);

  /// additive correction terms used by Andromeda (pscore + massC + cleaveC + modC - 100). For reference see the Andromeda source code.
  /// @note constants used in the correction term might be instrument dependent
  static double massCorrectionTerm(double mass);

  /// correction term for type of cleavage. For reference see the Andromeda source code.
  /// @note constants used in the correction term might be instrument dependent
  static double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage);

  /// correction term for modification. For reference see the Andromeda source code.
  /// @note constants used in the correction term might be instrument dependent
  static double modificationCorrectionTerm(Size modifications);
};

}

