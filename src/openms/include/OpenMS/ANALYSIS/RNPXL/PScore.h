// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_RNPXL_PSCORE
#define OPENMS_ANALYSIS_RNPXL_PSCORE

#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

struct OPENMS_DLLAPI PScore
{
  // calculates for each peak, how many neighboring peaks (in the given window) have higher intensity
  // the result can be used to efficiently filter spectra for top 1..n peaks in mass windows
  static std::vector<Size> calculateIntensityRankInMZWindow(const std::vector<double>& mz, const std::vector<double>& intensities, double mz_window);

  // used to precalculate peak ranks for a whole experiment using the calculateIntensityRankInMZWindow function
  static std::vector<std::vector<Size> > calculateRankMap(const PeakMap& peak_map, double mz_window = 100);

  // Calculates spectra for peak level between min_level to max_level and stores them in the map
  // A spectrum of peak level n retains the top n intensity peaks in a sliding mz_window centered at each peak
  // min and max level are taken from the Andromeda publication but are similar to the AScore publication
  static std::map<Size, PeakSpectrum > calculatePeakLevelSpectra(const PeakSpectrum& spec, const std::vector<Size>& ranks, Size min_level = 2, Size max_level = 10);

  // Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence
  // The best score obtained by scoring all those theoretical spectra against the experimental ones is returned.
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, PeakSpectrum>& peak_level_spectra, const std::vector<RichPeakSpectrum>& theo_spectra, double mz_window = 100.0);

  // Single spectrum version of above
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, PeakSpectrum>& peak_level_spectra, const RichPeakSpectrum& theo_spectrum, double mz_window = 100.0);

  // additive correction terms used by Andromeda (pscore + massC + cleaveC + modC - 100)
  static double massCorrectionTerm(double mass);

  static double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage);

  static double modificationCorrectionTerm(Size modifications);

};

}
#endif

