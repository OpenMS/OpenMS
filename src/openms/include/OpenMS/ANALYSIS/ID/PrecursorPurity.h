// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <vector>


namespace OpenMS
{
  /**
      @brief Precursor purity or noise estimation

      This class computes metrics for precursor isolation window purity (or noise).
      The function extracts the peaks from an isolation window targeted for fragmentation
      and determines which peaks are isotopes of the target and which come from other sources.
      The intensities of the assumed target peaks are summed up as the target intensity.
      Using this information it calculates an intensity ratio for the relative intensity of the target
      compared to other sources.
      These metrics are combined over the previous and the next MS1 spectrum.
      @note: If an MS1 spectrum does not contain the target peak within the given tolerance, all values are returned as 0.

      @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI PrecursorPurity
  {

  public:

    struct PurityScores
    {
      double total_intensity = 0.0;
      double target_intensity = 0.0;
      double signal_proportion = 0.0;
      Size target_peak_count = 0;
      Size interfering_peak_count = 0;
      PeakSpectrum interfering_peaks; // peaks left after precursor (isotopic) peaks have been removed
    };

    /** @brief compute precursor purity metrics for each MS2 spectrum in a PeakMap
       This is the main function of this class. See class description.
       Note: Spectra annotated with charge 0 will be treated as charge 1.       

     * @param spectra A PeakMap containing MS1 and MS2 spectra in order of acquisition or measurement. The first spectrum must be an MS1.
     * @param precursor_mass_tolerance The precursor tolerance. Is used for determining the targeted peak and deisotoping.
     * @param precursor_mass_tolerance_unit_ppm The unit of the precursor tolerance
     * @param ignore_missing_precursor_spectra Allow MS2 spectra without a MS1 precursor spectrum (PurityScores for these spectra will be 0).
    */
    static std::map<String, PurityScores> computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, bool ignore_missing_precursor_spectra = false);

    /** @brief compute precursor purity metrics for one MS2 precursor

       @note This function is implemented in a general way and can also be used for e.g. MS3 precursor isolation windows in MS2 spectra.
             Spectra annotated with charge 0 will be treated as charge 1.

      @param ms1 The Spectrum containing the isolation window
      @param pre The precursor containing the definition the isolation window
      @param precursor_mass_tolerance The precursor tolerance. Is used for determining the targeted peak and deisotoping.
      @param precursor_mass_tolerance_unit_ppm The unit of the precursor tolerance
    */
    static PurityScores computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm);


  private:
    // simple helper to combine the metrics contained in two PurityScores
    static PurityScores combinePrecursorPurities(const PrecursorPurity::PurityScores& score1, const PrecursorPurity::PurityScores& score2);

  };
}
