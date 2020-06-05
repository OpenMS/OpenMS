// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

      @ingroup ID
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
      Size residual_peak_count = 0;
    };

    /** @brief compute precursor purity metrics for each MS2 spectrum in a PeakMap
       This is the main function of this class. See class description.

     * @param spectra A PeakMap containing MS1 and MS2 spectra in order of acquisition or measurement. The first spectrum must be an MS1.
     * @param precursor_mass_tolerance The precursor tolerance. Is used for determining the targeted peak and deisotoping.
     * @param precursor_mass_tolerance_unit_ppm The unit of the precursor tolerance
    */
    static std::map<String, PurityScores> computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);

    /** @brief compute precursor purity metrics for one MS2 precursor

       @note This function is implemented in a general way and can also be used for e.g. MS3 precursor isolation windows in MS2 spectra

     * @param ms1 The Spectrum containing the isolation window
     * @param pre The precursor containing the definition the isolation window
     * @param precursor_mass_tolerance The precursor tolerance. Is used for determining the targeted peak and deisotoping.
     * @param precursor_mass_tolerance_unit_ppm The unit of the precursor tolerance
    */
    static PurityScores computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm);


  private:
    // simple helper to combine the metrics contained in two PurityScores
    static PurityScores combinePrecursorPurities(const PrecursorPurity::PurityScores& score1, const PrecursorPurity::PurityScores& score2);

  };
}
