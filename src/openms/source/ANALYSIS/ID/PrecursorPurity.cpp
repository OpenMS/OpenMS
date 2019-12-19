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

// #ifdef _OPENMP
// #include <omp.h>
// #endif

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

    PeakSpectrum isolated_window;
    while (lower_it != upper_it)
    {
      isolated_window.push_back(*lower_it);
      lower_it++;
    }

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

    // deisotoping (try to find isotopic peaks of the precursor mass, even if the actual precursor peak is missing)
    while (true) // runs as long as the next mz is within the isolation window
    {
      double next_peak = target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge);

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

        isolated_window.erase(isolated_window.begin()+next_iso_index);
        target_peak_count++;
      }
      // always increment iso to progress the loop
      iso++;
    }

    double rel_sig(0);
    if (target_intensity > 0.0)
    {
      rel_sig = target_intensity / total_intensity;
    }

    score.total_intensity = total_intensity;
    score.target_intensity = target_intensity;
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
    if (score.target_intensity > 0.0) // otherwise default value of 0 is used
    {
      score.signal_proportion = score.target_intensity / score.total_intensity;
    }
    score.target_peak_count = score1.target_peak_count + score2.target_peak_count;
    score.residual_peak_count = score1.residual_peak_count + score2.residual_peak_count;

    return score;
  }

  std::map<String, PrecursorPurity::PurityScores> PrecursorPurity::computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
  {
    std::map<String, PrecursorPurity::PurityScores> purityscores;
    std::pair<std::map<String, PrecursorPurity::PurityScores>::iterator, bool> insert_return_value;

    if (spectra[0].getMSLevel() != 1)
    {
      OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. First Spectrum is not MS1. Precursor Purity info will not be calculated!\n";
      return purityscores;
    }
    if (spectra[0].getNativeID().empty())
    {
      OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. Spectrum without an ID. Precursor Purity info will not be calculated!\n";
      return purityscores;
    }

    // keep the index of the two MS1 spectra flanking the current group of MS2 spectra
    int current_parent_index = 0;
    int next_parent_index = 0;
    bool lastMS1(false);
    int spectra_size = static_cast<int>(spectra.size());
    bool break_loop(false);

// TODO make this multi-threaded, after making sure it is thread safe
// current_parent_index and next_parent_index variables are not thread safe, determine flanking MS1 spectra directly from the current MS2
//#pragma omp parallel for schedule(guided)
    for (int i = 0; i < spectra_size; ++i)
    {
      // cannot return or break from OpenMP for loop, so continue to the end and return empty vector after loop
      if (break_loop)
      {
        continue;
      }
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
          for (int j = i+1; j < spectra_size; ++j)
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

        PrecursorPurity::PurityScores score;
        PrecursorPurity::PurityScores score1 = PrecursorPurity::computePrecursorPurity(spectra[current_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        // use default values of 0, if there is an MS2 spectrum without an MS1 after it (may happen for the last group of MS2 spectra)
        PrecursorPurity::PurityScores score2;
        if (!lastMS1) // there is an MS1 after this MS2, combine values from previous and next MS1
        {
          score2 = PrecursorPurity::computePrecursorPurity(spectra[next_parent_index], spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
          score = PrecursorPurity::combinePrecursorPurities(score1, score2);
        }
        else // no later MS1 spectrum, use only numbers from previous MS1
        {
          score = score1;
        }

//#pragma omp critical (purityscores_access)
        {
          insert_return_value = purityscores.insert(std::pair<String, PrecursorPurity::PurityScores>(spectra[i].getNativeID(), score));
          // check if the key had already been used before, and only set break_loop if really necessary
          // (do not set it to false every time, might overwrite a true value in parallel processing)
          if (!insert_return_value.second)
          {
            break_loop = true;
          }
        }
        if (break_loop)
        {
          OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. Duplicate Spectrum IDs. Precursor Purity info will not be calculated!\n";
        }
      } // end of MS2 spectrum
    } // end of parallelized spectra loop
    if (break_loop)
    {
      return std::map<String, PrecursorPurity::PurityScores>();
    }
    else
    {
      return purityscores;
    }
  } // end of function def

}
