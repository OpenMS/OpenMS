// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <vector>

namespace OpenMS
{

/**
 *  @brief An implementation of the X!Tandem HyperScore PSM scoring function
 */               

struct OPENMS_DLLAPI HyperScore
{
  typedef std::pair<Size, double> IndexScorePair; 

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param exp_spectrum measured spectrum
   * @param theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   */
//  static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const RichPeakSpectrum& theo_spectrum);

  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum);

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  overload that returns some additional information on the match
   */
  struct PSMDetail
  {
    size_t matched_b_ions = 0;
    size_t matched_y_ions = 0;
    double mean_error = 0.0;
  };

  static double computeWithDetail(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum,
                        PSMDetail& d
                       );

  private:
    /// helper to compute the log factorial
    static double logfactorial_(const int x, int base = 2);
};

}
