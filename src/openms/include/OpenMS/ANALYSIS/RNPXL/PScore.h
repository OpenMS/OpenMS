// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
  /**
    @brief Implementation of the PScore PSM scoring algorithm
  */

struct OPENMS_DLLAPI PScore
{
  /* @brief calculate local (windowed) peak ranks.
   * The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
   * The result can be used to efficiently filter spectra for top 1..n peaks in mass windows 
   * @note: ranks are zero based (highest intensity peak in window has rank 0)
   * @param mz m/z positions of the peaks
   * @param intensities of the peaks
   * @param mz_window window in Thomson centered at each peak
   */  
  static std::vector<Size> calculateIntensityRankInMZWindow(const std::vector<double>& mz, const std::vector<double>& intensities, double mz_window);

  /* @brief precalculated, windowed peak ranks for a whole experiment. 
   * The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
   * 1. Each spectrum is subdivided into windows of size @param mz_window.
   * 2. For each window, peak ranks are assigned using calculateIntensityRankInMZWindow.
   * 3. A rank map is returned
   * @note: ranks are zero based (top element has rank 0)
   * @param peak_map Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra.
   * @param mz_window window in Thomson centered at each peak
   */
  template <typename PeakType>
  static std::vector<std::vector<Size> > calculateRankMap(const MSExperiment<PeakType>& peak_map, double mz_window = 100)
  {
    std::vector<std::vector<Size> > rank_map; // note: ranks are zero based
    rank_map.reserve(peak_map.size());
    for (Size i = 0; i != peak_map.size(); ++i)
    {
      const MSSpectrum<PeakType>& spec = peak_map[i];
      std::vector<double> mz;
      std::vector<double> intensities;
      for (Size j = 0; j != spec.size(); ++j)
      {
        mz.push_back(spec[j].getMZ());
        intensities.push_back(spec[j].getIntensity());
      }
      rank_map.push_back(calculateIntensityRankInMZWindow(mz, intensities, mz_window));
    }
    return rank_map;
  }

  static double computeCumulativeScore_(Size N, Size n, double p);

  /* @brief Calculates spectra for peak level between min_level to max_level and stores them in the map
   * A spectrum of peak level n retains the (n+1) top intensity peaks in a sliding mz_window centered at each peak.
   * @note: levels are zero based (level 0 has only the top intensity peaks for each window, level 1 the top and second most intensive one)
   * @note: min and max level are taken from the Andromeda publication but are similar to the AScore publication
   */ 
  template <typename SpectrumType>
  static std::map<Size, SpectrumType > calculatePeakLevelSpectra(const SpectrumType& spec, const std::vector<Size>& ranks, Size min_level = 1, Size max_level = 9)
  {
    typename std::map<Size, SpectrumType > peak_level_spectra;

    if (spec.empty()) return peak_level_spectra;

    // loop over all peaks and associated (zero-based) ranks
    for (Size i = 0; i != ranks.size(); ++i)
    {
      // start at the highest (less restrictive) level
      for (int j = static_cast<int>(max_level); j >= static_cast<int>(min_level); --j)
      {
        // if the current peak is annotated to have lower or equal rank then allowed for this peak level add it
        if (static_cast<int>(ranks[i]) <= j)
        {
          peak_level_spectra[j].push_back(spec[i]);
        }
        else
        {
          // if the current peak has higher rank than the current level then all it is also to high for the lower levels
          break;
        }
      }
    }
    return peak_level_spectra;
  }

  /* @brief Computes the PScore for a vector of theoretical spectra
   * Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
   * The best score obtained by scoring all those theoretical spectra against the experimental ones is returned.
   * @param fragment_mass_tolerance mass tolerance for matching peaks
   * @param fragment_mass_tolerance_unit_ppm whether Thomson or ppm is used
   * @param peak_level_spectra spectra for different peak levels (=filtered by maximum rank).
   * @param theo_spectra theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
   * @param mz_window window in Thomson centered at each peak
   */ 
  template <typename SpectrumType>
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, SpectrumType>& peak_level_spectra, const std::vector<RichPeakSpectrum>& theo_spectra, double mz_window = 100.0)
  {
//    OpenMS::AScore a_score_algorithm; // TODO: make the cumulative score function static

    double best_pscore = 0.0;

    for (typename std::vector<SpectrumType>::const_iterator theo_spectra_it = theo_spectra.begin(); theo_spectra_it != theo_spectra.end(); ++theo_spectra_it)
    {
      const RichPeakSpectrum& theo_spectrum = *theo_spectra_it;

      // number of theoretical ions for current spectrum
      Size N = theo_spectrum.size();

      for (typename std::map<Size, SpectrumType>::const_iterator l_it = peak_level_spectra.begin(); l_it != peak_level_spectra.end(); ++l_it)
      {
        const double level = static_cast<double>(l_it->first);
        const PeakSpectrum& exp_spectrum = l_it->second;

        Size matched_peaks(0);
        for (RichPeakSpectrum::ConstIterator theo_peak_it = theo_spectrum.begin(); theo_peak_it != theo_spectrum.end(); ++theo_peak_it)
        {
          const double& theo_mz = theo_peak_it->getMZ();

          double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);
          double exp_mz = exp_spectrum[index].getMZ();

          // found peak match
          if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            ++matched_peaks;
          }
        }

        // compute p score as e.g. in the AScore implementation or Andromeda
        const double p = level / mz_window;
        const double pscore = -10.0 * log10(computeCumulativeScore_(N, matched_peaks, p));
        if (pscore > best_pscore)
        {
          best_pscore = pscore;
        }
      }
    }

    return best_pscore;
  }

  /* @brief Computes the PScore for a single theoretical spectrum
   * @param fragment_mass_tolerance mass tolerance for matching peaks
   * @param fragment_mass_tolerance_unit_ppm whether Thomson or ppm is used
   * @param peak_level_spectra spectra for different peak levels (=filtered by maximum rank).
   * @param theo_spectra theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
   * @param mz_window window in Thomson centered at each peak
   */ 
  template <typename SpectrumType>
  static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map<Size, SpectrumType>& peak_level_spectra, const RichPeakSpectrum& theo_spectrum, double mz_window = 100.0)
  {
//    AScore a_score_algorithm; // TODO: make the cumulative score function static

    double best_pscore = 0.0;

    // number of theoretical ions for current spectrum
    Size N = theo_spectrum.size();

    for (typename std::map<Size, SpectrumType>::const_iterator l_it = peak_level_spectra.begin(); l_it != peak_level_spectra.end(); ++l_it)
    {
      const double level = static_cast<double>(l_it->first);
      const SpectrumType& exp_spectrum = l_it->second;

      Size matched_peaks(0);
      for (RichPeakSpectrum::ConstIterator theo_peak_it = theo_spectrum.begin(); theo_peak_it != theo_spectrum.end(); ++theo_peak_it)
      {
        const double& theo_mz = theo_peak_it->getMZ();

        double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size index = exp_spectrum.findNearest(theo_mz);
        double exp_mz = exp_spectrum[index].getMZ();

        // found peak match
        if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
        {
          ++matched_peaks;
        }
      }
      // compute p score as e.g. in the AScore implementation or Andromeda
      const double p = (level + 1) / mz_window;
      const double pscore = -10.0 * log10(computeCumulativeScore_(N, matched_peaks, p));

      if (pscore > best_pscore)
      {
        best_pscore = pscore;
      }
    }

    return best_pscore;
  }

  /// additive correction terms used by Andromeda (pscore + massC + cleaveC + modC - 100). For reference see the Andromeda source code.
  /// @note: constants used in the correction term might be instrument dependent
  static double massCorrectionTerm(double mass);

  /// correction term for type of cleavage. For reference see the Andromeda source code.
  /// @note: constants used in the correction term might be instrument dependent
  static double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage);

  /// correction term for modification. For reference see the Andromeda source code.
  /// @note: constants used in the correction term might be instrument dependent
  static double modificationCorrectionTerm(Size modifications);
};

}
#endif

