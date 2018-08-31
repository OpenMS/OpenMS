// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>

#include <OpenMS/CONCEPT/Constants.h>

// #include <iostream>

namespace OpenMS
{

  PrecursorPurity::PurityScores PrecursorPurity::computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm)
  {
    PrecursorPurity::PurityScores score;
    double target_mz = pre.getMZ();
    double lower = target_mz - pre.getIsolationWindowLowerOffset();
    double upper = target_mz + pre.getIsolationWindowUpperOffset();
    int charge = pre.getCharge();

    double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance * 1e-6) : precursor_mass_tolerance;

    auto lower_it = ms1.MZBegin(lower);
    auto upper_it = ms1.MZEnd(upper);

    // std::cout << "MS1: " << ms1.getNativeID() << " | charge: " << charge << " | lower: " << lower << " | target: " << target_mz << " | upper: " << upper << std::endl;
    // std::cout << "lower_it: " << (*lower_it).getMZ() << " | upper_it: " << (*upper_it).getMZ() << std::endl;
    // std::cout << "lower offset: " << pre.getIsolationWindowLowerOffset() << " | upper offset: " << pre.getIsolationWindowUpperOffset() << std::endl;

    PeakSpectrum isolated_window;
    while (lower_it != upper_it)
    {
      isolated_window.push_back(*lower_it);
      lower_it++;
    }

    // std::cout << "Isolation window peaks: " << isolated_window.size();
    // for (auto peak : isolated_window)
    // {
    //   std::cout << " | " << peak.getMZ();
    // }
    // std::cout << std::endl;

    // total intensity in isolation window
    double total_intensity(0);
    for (auto peak : isolated_window)
    {
      total_intensity += peak.getIntensity();
    }

    // search for the target peak, return scores with 0-values if it is not found
    if (isolated_window.empty())
    {
      return score;
    }
    int target_index = isolated_window.findNearest(target_mz, precursor_tolerance_abs);
    if (target_index == -1)
    {
      return score;
    }
    double target_intensity(0);
    target_intensity = isolated_window[target_index].getIntensity();
    isolated_window.erase(isolated_window.begin()+target_index);
    if (target_intensity == 0.0)
    {
      return score;
    }

    // deisotoping
    // std::cout << "iso peaks: ";
    Size target_peak_count(1);
    bool next_peak_found = true;
    double iso = 1;
    while(next_peak_found)
    {
      double next_peak = target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge);
      int next_iso_index = isolated_window.findNearest(next_peak, precursor_tolerance_abs);
      if (next_iso_index == -1)
      {
        next_peak_found = false;
      }
      else
      {
        target_intensity += isolated_window[next_iso_index].getIntensity();

        // std::cout << isolated_window[next_iso_index].getMZ() << " | ";

        isolated_window.erase(isolated_window.begin()+next_iso_index);
        iso++;
        target_peak_count++;
      }
    }

    next_peak_found = true;
    iso = -1;
    while(next_peak_found)
    {
      double next_peak = target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge);
      int next_iso_index = isolated_window.findNearest(next_peak, precursor_tolerance_abs);
      if (next_iso_index == -1)
      {
        next_peak_found = false;
      }
      else
      {
        target_intensity += isolated_window[next_iso_index].getIntensity();

        // std::cout << isolated_window[next_iso_index].getMZ() << " | ";

        isolated_window.erase(isolated_window.begin()+next_iso_index);
        iso--;
        target_peak_count++;
      }
    }
    // std::cout << std::endl;

    // std::cout << "noise peaks: ";
    double noise_intensity(0);
    for (auto peak : isolated_window)
    {
      noise_intensity += peak.getIntensity();
      // std::cout << peak.getMZ() << " | ";
    }
    // std::cout << std::endl;

    double rel_sig = target_intensity / total_intensity;

    score.total_intensity = total_intensity;
    score.target_intensity = target_intensity;
    score.residual_intensity = noise_intensity;
    score.signal_proportion = rel_sig;
    score.target_peak_count = target_peak_count;
    score.residual_peak_count = isolated_window.size();

    return score;
  }

  PrecursorPurity::PurityScores PrecursorPurity::combinePrecursorPurities(const PrecursorPurity::PurityScores& score1, const PrecursorPurity::PurityScores& score2)
  {
    PrecursorPurity::PurityScores score;
    score.total_intensity = score1.total_intensity + score2.total_intensity;
    score.target_intensity = score1.target_intensity + score2.target_intensity;
    score.residual_intensity = score1.residual_intensity + score2.residual_intensity;
    if (score.target_intensity > 0.0) // otherwise default value of 0 is used
    {
      score.signal_proportion = score.target_intensity / score.total_intensity;
    }
    score.target_peak_count = score1.target_peak_count + score2.target_peak_count;
    score.residual_peak_count = score1.residual_peak_count + score2.residual_peak_count;

    return score;
  }

  std::vector<PrecursorPurity::PurityScores> PrecursorPurity::computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
  {
    std::vector<PrecursorPurity::PurityScores> purityscores;

    Size current_parent_index = 0;
    Size next_parent_index = 0;
    for (Size i = 0; i < spectra.size(); ++i)
    {
      if (spectra[i].getMSLevel() == 1)
      {
        current_parent_index = i;
      }
      else if (spectra[i].getMSLevel() == 2)
      {
        if (next_parent_index < i)
        {
          for (Size j = i+1; j < spectra.size(); ++j)
          {
            if (spectra[j].getMSLevel() == 1)
            {
              next_parent_index = j;
              break;
            }
          }
        }

        PrecursorPurity::PurityScores score1 = PrecursorPurity::computePrecursorPurity(spectra[current_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        PrecursorPurity::PurityScores score2 = PrecursorPurity::computePrecursorPurity(spectra[next_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        PrecursorPurity::PurityScores score = PrecursorPurity::combinePrecursorPurities(score1, score2);
        purityscores.push_back(score);

        // std::cout << "Score1 | Spectrum: " << i << " | total intensity: " << score1.total_intensity << " | target intensity: " << score1.target_intensity << " | noise intensity: " << score1.residual_intensity << " | rel_sig: " << score1.signal_proportion << std::endl;
        // std::cout << "Score2 | Spectrum: " << i << " | total intensity: " << score2.total_intensity << " | target intensity: " << score2.target_intensity << " | noise intensity: " << score2.residual_intensity << " | rel_sig: " << score2.signal_proportion << std::endl;
        // std::cout << "Combin | Spectrum: " << i << " | total intensity: " << score.total_intensity << " | target intensity: " << score.target_intensity << " | noise intensity: " << score.residual_intensity << " | rel_sig: " << score.signal_proportion << std::endl;

      } // end of MS2 spectrum
    } // spectra loop
    return purityscores;
  } // end of function def

}
