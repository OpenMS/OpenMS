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
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  PrecursorPurity::PurityScores PrecursorPurity::computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm)
  {
    PrecursorPurity::PurityScores score;
    double target_mz = pre.getMZ();
    double lower = target_mz - pre.getIsolationWindowLowerOffset();
    double upper = target_mz + pre.getIsolationWindowUpperOffset();
    int charge = pre.getCharge();

    double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance*2 * 1e-6) : precursor_mass_tolerance*2;

    auto lower_it = ms1.MZBegin(lower);
    auto upper_it = ms1.MZEnd(upper);

    // std::cout << "charge: " << charge << " | lower: " << lower << " | target: " << target_mz << " | upper: " << upper << std::endl;
    // std::cout << "lower offset: " << pre.getIsolationWindowLowerOffset() << " | upper offset: " << pre.getIsolationWindowUpperOffset() << std::endl;
    // std::cout << "lower peak: " << (*lower_it).getMZ() << " | upper peak: " << (*upper_it).getMZ() << std::endl;

    PeakSpectrum isolated_window;
    while (lower_it != upper_it)
    {
      isolated_window.push_back(*lower_it);
      lower_it++;
    }

    // std::cout << "Isolation window peaks: " << isolated_window.size();
    // for (const auto& peak : isolated_window)
    // {
    //   std::cout << " | " << peak.getMZ();
    // }
    // std::cout << std::endl;

    // total intensity in isolation window
    double total_intensity(0);
    double target_intensity(0);
    Size target_peak_count(0);
    for (const auto& peak : isolated_window)
    {
      total_intensity += peak.getIntensity();
    }

    // search for the target peak, return scores with 0-values if it is not found
    if (isolated_window.empty())
    {
      return score;
    }

    // estimate a lower boundary for isotopic peaks
    int negative_isotopes((pre.getIsolationWindowLowerOffset() * charge));
    double iso = -negative_isotopes;
    // depending on the isolation window, the first estimated peak might be outside the window
    if (target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge) < lower)
    {
      iso++;
    }

    // std::cout << "target peaks: ";
    // deisotoping (try to find isotopic peaks of the precursor mass, even if the actual precursor peak is missing)
    while (true) // runs as long as the next mz is within the isolation window
    {
      double next_peak = target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge);
      // std::cout << iso << " : iso | " << next_peak << " : next peak | ";

      // stop loop when new mz is outside the isolation window
      // changes through the isotope index iso
      if (next_peak > upper)
      {
        break;
      }
      int next_iso_index = isolated_window.findNearest(next_peak, precursor_tolerance_abs);
      if (next_iso_index != -1)
      {
        target_intensity += isolated_window[next_iso_index].getIntensity();

        // std::cout << isolated_window[next_iso_index].getMZ() << " : matched | ";

        isolated_window.erase(isolated_window.begin()+next_iso_index);
        target_peak_count++;
      }
      // always increment iso to progress the loop
      iso++;
    }
    // std::cout << std::endl;

    // // std::cout << "noise peaks: ";
    // double noise_intensity(0);
    // for (const auto& peak : isolated_window)
    // {
    //   noise_intensity += peak.getIntensity();
    //   // std::cout << peak.getMZ() << " | ";
    // }
    // // std::cout << std::endl;

    double rel_sig(0);
    if (target_intensity > 0.0)
    {
      rel_sig = target_intensity / total_intensity;
    }
    double noise_intensity = total_intensity - target_intensity;

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

    // TODO throw an exception or at least a warning?
    // if there is no MS1 before the first MS2, the spectra datastructure is not suitable for this function
    if (spectra[0].getMSLevel() == 2)
    {
      LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. Will be skipped!\n";
      return purityscores;
    }

    // keep the index of the two MS1 spectra flanking the current group of MS2 spectra
    Size current_parent_index = 0;
    Size next_parent_index = 0;
    bool lastMS1(false);
    for (Size i = 0; i < spectra.size(); ++i)
    {
      // change current parent index if a new MS1 spectrum is reached
      if (spectra[i].getMSLevel() == 1)
      {
        current_parent_index = i;
      }
      else if (spectra[i].getMSLevel() == 2)
      {
        // update next MS1 index, if it is lower than the current MS2 index
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
          // if the next MS1 index was not updated,
          // the end of the PeakMap was reached right after this current group of MS2 spectra
          if (next_parent_index < i)
          {
            lastMS1 = true;
          }
        }

        // std::cout << "MS1 Spectrum: " << spectra[current_parent_index].getNativeID() << " | MS2 : " << spectra[i].getNativeID() << std::endl;
        PrecursorPurity::PurityScores score1 = PrecursorPurity::computePrecursorPurity(spectra[current_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        // use default values of 0, if there is an MS2 spectrum without an MS1 after it (may happen for the last group of MS2 spectra)
        PrecursorPurity::PurityScores score2;
        if (!lastMS1) // there is an MS1 after this MS2
        {
          score2 = PrecursorPurity::computePrecursorPurity(spectra[next_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        }
        PrecursorPurity::PurityScores score = PrecursorPurity::combinePrecursorPurities(score1, score2);
        purityscores.push_back(score);

        // std::cout << "Score1 | Spectrum: " << i << " | total intensity: " << score1.total_intensity << " | target intensity: " << score1.target_intensity << " | noise intensity: " << score1.residual_intensity << " | rel_sig: " << score1.signal_proportion << std::endl;
        // std::cout << "Score2 | Spectrum: " << i << " | total intensity: " << score2.total_intensity << " | target intensity: " << score2.target_intensity << " | noise intensity: " << score2.residual_intensity << " | rel_sig: " << score2.signal_proportion << std::endl;
        // std::cout << "Combin | Spectrum: " << i << " | total intensity: " << score.total_intensity << " | target intensity: " << score.target_intensity << " | noise intensity: " << score.residual_intensity << " | rel_sig: " << score.signal_proportion << std::endl;
        // std::cout << "#################################################################################################################" << std::endl;

      } // end of MS2 spectrum
    } // spectra loop
    return purityscores;
  } // end of function def

}
