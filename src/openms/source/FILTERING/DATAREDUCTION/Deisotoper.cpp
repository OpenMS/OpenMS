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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{

// static
void Deisotoper::deisotopeAndSingleCharge(MSSpectrum& spec,
                      double fragment_tolerance,
                      bool fragment_unit_ppm,
                      int min_charge,
                      int max_charge,
                      bool keep_only_deisotoped,
                      unsigned int min_isopeaks,
                      unsigned int max_isopeaks,
                      bool make_single_charged,
                      bool annotate_charge,
                      bool annotate_iso_peak_count,
                      bool use_decreasing_model,
                      bool add_up_intensity)
{
  OPENMS_PRECONDITION(spec.isSorted(), "Spectrum must be sorted.");

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__,
		    __LINE__,
		    OPENMS_PRETTY_FUNCTION,
		    "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spec.empty()) { return; }

  Size charge_index(0);
  Size iso_peak_count_index(0);

  // reserve integer data array to store charge of peaks
  if (annotate_charge)
  {
    // expand to hold one additional integer data array to hold the charge
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("charge");
    charge_index = spec.getIntegerDataArrays().size()-1;
  }
  // reserve integer data array to store number of isotopic peaks for each isotopic pattern
  if (annotate_iso_peak_count)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("iso_peak_count");
    iso_peak_count_index = spec.getIntegerDataArrays().size()-1;
  }

  // during discovery phase, work on a constant reference (just to make sure we do not modify spec)
  const MSSpectrum& old_spectrum = spec;

  // determine charge seeds and extend them
  std::vector<size_t> mono_isotopic_peak(old_spectrum.size(), 0);
  std::vector<int> features(old_spectrum.size(), -1);
  std::vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
  std::vector<Size> iso_peak_count(old_spectrum.size(), 1);
  int feature_number = 0;

  std::vector<size_t> extensions;

  for (size_t current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    const double current_mz = old_spectrum[current_peak].getMZ();
    if (add_up_intensity)
    {
      mono_iso_peak_intensity[current_peak] = old_spectrum[current_peak].getIntensity();
    }

    for (int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
    {
      // TODO stop looking for charges q with m/q < precursor_mass/q

      // try to extend isotopes from mono-isotopic peak
      // if extension larger then min_isopeaks possible:
      //   - save charge q in mono_isotopic_peak[]
      //   - annotate_charge all isotopic peaks with feature number
      if (features[current_peak] == -1) // only process peaks which have no assigned feature number
      {
        bool has_min_isopeaks = true;
        const double tolerance_dalton = fragment_unit_ppm ? Math::ppmToMass(fragment_tolerance, current_mz) : fragment_tolerance;
        extensions.clear();
        extensions.push_back(current_peak);
        for (unsigned int i = 1; i < max_isopeaks; ++i)
        {
          const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
          const int p = old_spectrum.findNearest(expected_mz, tolerance_dalton);
          if (p == -1) // test for missing peak
          {
            has_min_isopeaks = (i >= min_isopeaks);
            break;
          }
          else
          {
            // Possible improvement: include proper averagine model filtering. for now start at the second peak to test hypothesis
            // Note: this is a common approach used in several other search engines
            if (use_decreasing_model && (old_spectrum[p].getIntensity() > old_spectrum[extensions.back()].getIntensity()))
            {
              has_min_isopeaks = (i >= min_isopeaks);
              break;
            }
            // averagine check passed or skipped
            extensions.push_back(p);
            if (annotate_iso_peak_count)
            {
              iso_peak_count[current_peak] = i + 1; // with "+ 1" the monoisotopic peak is counted as well
            }
          }
        }

        if (has_min_isopeaks)
        {
          // std::cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << std::endl;
          mono_isotopic_peak[current_peak] = q;
          for (unsigned int i = 0; i != extensions.size(); ++i)
          {
            features[extensions[i]] = feature_number;
            // monoisotopic peak intensity is already set above, add up the other intensities here
            if (add_up_intensity && (i != 0))
            {
              mono_iso_peak_intensity[current_peak] += old_spectrum[extensions[i]].getIntensity();
            }
          }
          ++feature_number;
        }
      }
    }
  }

  // apply changes, i.e. select the indices which should survive
  std::vector<Size> select_idx;

  for (size_t i = 0; i != spec.size(); ++i)
  {
    Size z = mono_isotopic_peak[i];
    if (annotate_charge)
    {
      spec.getIntegerDataArrays()[charge_index].push_back((int)z);
    }
    if (annotate_iso_peak_count)
    {
      spec.getIntegerDataArrays()[iso_peak_count_index].push_back((int)iso_peak_count[i]);
    }
    if (add_up_intensity)
    {
      spec[i].setIntensity(mono_iso_peak_intensity[i]);
    }

    if (!keep_only_deisotoped)
    { // keep all unassigned peaks
      if (features[i] < 0)
      {
        select_idx.push_back(i);
        continue;
      }
    }

    if (z == 0) continue;

    // convert mono-isotopic peak with charge assigned by deisotoping
    if (make_single_charged)
    {
      spec[i].setMZ(spec[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
    }
    select_idx.push_back(i);

  }
  // properly subsets all datapoints (incl. dataArrays)
  spec.select(select_idx);
  spec.sortByPosition();
  return;
}
} // namespace
