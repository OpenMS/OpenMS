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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/PScore.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>

#include <vector>
#include <map>

using std::map;
using std::vector;

namespace OpenMS
{

  vector<Size> PScore::calculateIntensityRankInMZWindow(const vector<double>& mz, const vector<double>& intensities, double mz_window = 100)
  {
    vector<Size> ranks;
    if (mz.empty())
    {
      return ranks;
    }
    ranks.reserve(mz.size());

    const double half_window = mz_window / 2.0;
    for (Size p = 0; p < mz.size(); ++p)
    {
      const double m = mz[p];
      const double i = intensities[p];

      // determine rank
      Size rank(0);

      // count neighbors to the left that have higher intensity
      for (Int j = p - 1; j >= 0; --j)
      {
        if (mz[j] < m - half_window) break;
        if (intensities[j] > i) ++rank;
      }

      // count neighbors to the right that have higher intensity
      for (Size j = p + 1; j < mz.size(); j++)
      {
        if (mz[j] > m + half_window) break;
        if (intensities[j] > i) ++rank;
      }
      ranks.push_back(rank);
    }

    return ranks;
  }

  map<Size, MSSpectrum<Peak1D> > PScore::calculatePeakLevelSpectra(const MSSpectrum<Peak1D>& spec, double mz_window = 100, Size min_level = 2, Size max_level = 10)
  {
    map<Size, MSSpectrum<Peak1D> > peak_level_spectra;

    if (spec.empty()) return peak_level_spectra;

    vector<double> mz;
    vector<double> intensities;
    mz.reserve(spec.size());
    intensities.reserve(spec.size());

    for (Size i = 0; i != spec.size(); ++i)
    {
      mz.push_back(spec[i].getMZ());
      intensities.push_back(spec[i].getIntensity());
    }

    vector<Size> ranks = calculateIntensityRankInMZWindow(mz, intensities, mz_window);

    // loop over all peaks and associated ranks
    for (Size i = 0; i != ranks.size(); ++i)
    {
      // start at the highest (less restrictive) level
      for (Size j = max_level; j >= min_level; --j)
      {
        // if the current peak is annotated to have lower or equal rank then allowed for this peak level add it
        if (ranks[i] <= j)
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

  double PScore::computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const map<Size, MSSpectrum<Peak1D> >& peak_level_spectra, const vector< MSSpectrum<RichPeak1D> >& theo_spectra, double mz_window = 100.0)
  {
    AScore a_score_algorithm; // TODO: make the cumulative score function static

    double best_pscore = 0.0;

    for (vector< MSSpectrum<RichPeak1D> >::const_iterator theo_spectra_it = theo_spectra.begin(); theo_spectra_it != theo_spectra.end(); ++theo_spectra_it)
    {
      const MSSpectrum<RichPeak1D>& theo_spectrum = *theo_spectra_it;

      // number of theoretical ions for current spectrum
      Size N = theo_spectrum.size();

      for (map<Size, MSSpectrum<Peak1D> >::const_iterator l_it = peak_level_spectra.begin(); l_it != peak_level_spectra.end(); ++l_it)
      {
        const double level = static_cast<double>(l_it->first);
        const MSSpectrum<Peak1D>& exp_spectrum = l_it->second;

        Size matched_peaks(0);
        for (MSSpectrum<RichPeak1D>::ConstIterator theo_peak_it = theo_spectrum.begin(); theo_peak_it != theo_spectrum.end(); ++theo_peak_it)
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
        const double pscore = -10.0 * log10(a_score_algorithm.computeCumulativeScore(N, matched_peaks, p));
        if (pscore > best_pscore)
        {
          best_pscore = pscore;
        }
      }
    }

    return best_pscore;
  }

}

