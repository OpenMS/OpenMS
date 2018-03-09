// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_ANALYSIS_XLMS_OPXLSPECTRUMPROCESSINGALGORITHMS_H
#define OPENMS_ANALYSIS_XLMS_OPXLSPECTRUMPROCESSINGALGORITHMS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <numeric>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI OPXLSpectrumProcessingAlgorithms
  {
    public:

    /**
       * @brief Merges two spectra into one while correctly considering metainfo in DataArrays
       * @param first_spectrum
       * @param second_spectrum
       * @return A PeakSpectrum containing all peaks from both input spectra
       */
      static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

      /**
       * @brief Preprocesses spectra
       *
       * Filters out spectra with too few peaks (based on peptide_min_size) and those that do not fit into the precursor charge range.
       * Removes zero intensity peaks and normalizes intensities.
       * If the given tolerance is low enough, deisotoping is performed. Otherwise only the 500 most intense peaks are kept, if the param labeled is false.
       * The number of returned spectra is equal to the number of input spectra for labeled data (otherwise not necessarily).
       *
       * @param exp
       * @param fragment_mass_tolerance_xlinks
       * @param fragment_mass_tolerance_unit_ppm
       * @param peptide_min_size
       * @param min_precursor_charge
       * @param max_precursor_charge
       * @param labeled
       * @return A PeakMap of preprocessed spectra
       */
      static PeakMap preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance_xlinks, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, bool labeled);

      /**
       * @brief Computes a spectrum alignment while considering fragment charges stored in a IntegerDataArray and an intensity difference ratio
       * @param alignment The empty alignment, that will be filled by the algorithm
       * @param s1 The first spectrum to be aligned
       * @param s2 the second spectrum to be aligned
       * @param tolerance The peak mass tolerance
       * @param relative_tolerance True if the given tolerance is a ppm tolerance, false if tolerance is in Da
       * @param intensity_cutoff Peaks will only be aligned if intensity1 / intensity2 > intensity_cutoff, with intensity1 being the lower of the two compared peaks and intensity2 the higher one. Set to 0 to ignore intensity differences.
       */
      static void getSpectrumAlignment(std::vector <std::pair <Size, Size> >& alignment, const PeakSpectrum & s1, const PeakSpectrum & s2, double tolerance, bool relative_tolerance, DataArrays::FloatDataArray & ppm_error_array, double intensity_cutoff = 0.0);

      template <typename SpectrumType1, typename SpectrumType2>
          static void getSpectrumAlignmentFastCharge(
            std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance,
            bool fragment_mass_tolerance_unit_ppm,
            const SpectrumType1& theo_spectrum,
            const SpectrumType2& exp_spectrum,
            const typename SpectrumType1::IntegerDataArray& theo_charges,
            const typename SpectrumType2::IntegerDataArray& exp_charges)
          {
            OPENMS_PRECONDITION(exp_spectrum.isSorted(), "Spectrum needs to be sorted.");
            OPENMS_PRECONDITION(theo_spectrum.isSorted(), "Spectrum needs to be sorted.");
            OPENMS_PRECONDITION((alignment.empty() == true), "Alignment result vector needs to be empty.");

            const Size n_t(theo_spectrum.size());
            const Size n_e(exp_spectrum.size());
            const bool has_charge = !(exp_charges.empty() || theo_charges.empty());

            if (n_t == 0 || n_e == 0) { return; }

            Size t(0), e(0);
            alignment.reserve(theo_spectrum.size());

            while (t < n_t && e < n_e)
            {
              const double theo_mz = theo_spectrum[t].getMZ();
              const double exp_mz = exp_spectrum[e].getMZ();

              double tz(0), ez(0);
              if (has_charge)
              {
                tz = theo_charges[t];
                ez = exp_charges[e];
              }
              const bool tz_matches_ez = (ez == tz || !ez || !tz);

              double d = exp_mz - theo_mz;
              const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

              if (fabs(d) <= max_dist_dalton) // match in tolerance window?
              {
                // get first peak with matching charge in tolerance window
                if (!tz_matches_ez)
                {
                  Size e_candidate(e);
                  while (true)
                  {
                    ++e_candidate;
                    double new_ez = has_charge ? exp_charges[e_candidate] : 0;
                    const bool charge_matches = (new_ez == tz || !new_ez || !tz);
                    double new_d = exp_spectrum[e].getMZ() - theo_mz;
                    if (charge_matches && new_d <= max_dist_dalton)
                    { // found a match
                      break;
                    }
                    else if (new_d > max_dist_dalton)
                    { // no match found
                      e_candidate = e;
                      break;
                    }
                  }
                  if (e == e_candidate)
                  { // no match found continue with next theo. peak
                    ++t;
                    continue;
                  }
                  else
                  { // match found
                    e = e_candidate;
                  }
                }

                // Invariant: e now points to the first peak in tolerance window, that matches in charge

                // last peak? there can't be a better one in this tolerance window
                if (e >= n_e - 1) { alignment.emplace_back(std::make_pair(t, e)); return; }

                Size closest_exp_peak(e);

                // Invariant: closest_exp_peak always point to best match

                double new_d, new_ez(0);
                double best_d = exp_spectrum[closest_exp_peak].getMZ() - theo_mz;

      //        std::cerr << "first peak in window:" <<  exp_spectrum[closest_exp_peak].getMZ() << "\t theo: " << theo_mz << "\t best d:" << best_d << "\n";

                do // check for better match in tolerance window
                {
                  // advance to next exp. peak
                  ++e;

                  // determine distance of next peak
                  new_d = exp_spectrum[e].getMZ() - theo_mz;
                  const bool in_tolerance_window = (fabs(new_d) < max_dist_dalton);
      //          std::cerr << "  new peak:" << exp_spectrum[e].getMZ() << "\t theo: " << theo_mz << "\t new d: " << new_d << "\t in window:" << in_tolerance_window << "\n";

                  if (!in_tolerance_window) { break; }

                  // Invariant: e is in tolerance window

                  // check if charge of next peak matches
                  if (has_charge) { new_ez = exp_charges[e]; }
                  const bool charge_matches = (new_ez == tz || !new_ez || !tz);
                  if (!charge_matches) { continue; }

                  // Invariant: charge matches

                  const bool better_distance = (fabs(new_d) <= fabs(best_d));

                  // better distance (and matching charge)? better match found
                  if (better_distance)
                  { // found a better match
                    closest_exp_peak = e;
                    best_d = new_d;
                  }
                  else
                  { // distance got worse -> no additional matches!
                    break;
                  }
                }
                while (e < n_e - 1);

      //        std::cerr << "added peak:" << exp_spectrum[closest_exp_peak].getMZ() << "\t theo: " << theo_mz << "\n";

                // search in tolerance window for an experimental peak closer to theoretical one
                alignment.emplace_back(std::make_pair(t, closest_exp_peak));

                e = closest_exp_peak + 1;  // advance experimental peak to 1-after the best match
                ++t; // advance theoretical peak
              }
              else if (d < 0) // exp. peak is left of theo. peak (outside of tolerance window)
              {
                ++e;
              }
              else if (d > 0) // theo. peak is left of exp. peak (outside of tolerance window)
              {
                ++t;
              }
            }
      }

      /**
       * @brief Deisotopes a spectrum and stores the determined charges in an IntegerDataArray

          If keep_only_deisotoped is false, the peaks that could not be deisotoped are assigned the charge 0.
          If an isotopic pattern contains more peaks than max_isopeaks, the rest are ignored for the current pattern.

       * @param old_spectrum The spectrum to be deisotoped
       * @param min_charge Minimal charge to consider for the isotope patterns
       * @param max_charge Maximal charge to consider for the isotope patterns
       * @param fragment_tolerance The mass tolerance for matching peaks of an isotope pattern
       * @param fragment_tolerance_unit_ppm True, if the given tolerance is in ppm, false if it is in Da
       * @param keep_only_deisotoped True if the peaks that could not be deisotoped should be discarded
       * @param min_isopeaks The minimal number of consecutive peaks in an isotopic pattern, before it gets acknowledged as an isotopic pattern
       * @param max_isopeaks The maximal number of consecutive peaks in an isotopic pattern.
       * @param make_single_charged If true, all peaks with charges larger than 1 are replaced with peaks with their corresponding single charged MZ
       * @return A PeakSpectrum annotated with charges
       */
      static PeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false);

  };

}

#endif // OPXLSPECTRUMPROCESSINGALGORITHMS_H
