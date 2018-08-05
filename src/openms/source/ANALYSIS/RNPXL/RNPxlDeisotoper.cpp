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
#include <OpenMS/ANALYSIS/RNPXL/RNPxlDeisotoper.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{

// static
void Deisotoper::deisotopeAndSingleCharge(MSSpectrum& spectra, 
                      double fragment_tolerance, 
      				        bool fragment_unit_ppm, 
                      int min_charge, 
					            int max_charge,
                      bool keep_only_deisotoped, 
                      unsigned int min_isopeaks, 
      				        unsigned int max_isopeaks, 
                      bool make_single_charged,
                      bool annotate_charge)
{
  OPENMS_PRECONDITION(spectra.isSorted(), "Spectrum must be sorted.");

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__, 
		    __LINE__, 
		    OPENMS_PRETTY_FUNCTION, 
		    "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spectra.empty()) { return; }

  MSSpectrum old_spectrum = spectra;

  // reserve integer data array to store charge of peaks
  if (annotate_charge) 
  {
    // expand to hold one additional integer data array to hold the charge
    spectra.getIntegerDataArrays().resize(spectra.getIntegerDataArrays().size() + 1);
    spectra.getIntegerDataArrays().back().setName("charge");
  }

  // determine charge seeds and extend them
  std::vector<unsigned int> mono_isotopic_peak(old_spectrum.size(), 0);
  std::vector<int> features(old_spectrum.size(), -1);
  int feature_number = 0;

  for (unsigned int current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    const double current_mz = old_spectrum[current_peak].getMZ();

    for (int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
    {
      // try to extend isotopes from mono-isotopic peak
      // if extension larger then min_isopeaks possible:
      //   - save charge q in mono_isotopic_peak[]
      //   - annotate_charge all isotopic peaks with feature number
      if (features[current_peak] == -1) // only process peaks which have no assigned feature number
      {
        bool has_min_isopeaks = true;
        std::vector<unsigned int> extensions;
        for (unsigned int i = 0; i < max_isopeaks; ++i)
        {
          const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
          const unsigned int p = old_spectrum.findNearest(expected_mz);
          const double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * expected_mz * 1e-6 : fragment_tolerance;
          const double distance_to_closest = fabs(old_spectrum[p].getMZ() - expected_mz);
	        if (distance_to_closest > tolerance_dalton) // test for missing peak
          {
            if (i < min_isopeaks) { has_min_isopeaks = false;}
	          break;
          }
          else
          {
            // Possible improvement: include proper averagine model filtering. for now start at the second peak to test hypothesis
	          // Note: this is a common approach used in several other search engines
            unsigned int n_extensions = extensions.size();
            if (n_extensions != 0)
            {
              if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
              {
                if (i < min_isopeaks) { has_min_isopeaks = false; }
               break;
	            }
            }

            // averagine check passed
            extensions.push_back(p);
          }
        }

        if (has_min_isopeaks)
        {
	        // std::cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << std::endl;
          mono_isotopic_peak[current_peak] = q;
          for (unsigned int i = 0; i != extensions.size(); ++i)
          {
            features[extensions[i]] = feature_number;
          }
          ++feature_number;
        }
      }
    }
  }

  spectra.clear(false);
  for (unsigned int i = 0; i != old_spectrum.size(); ++i)
  {
    int z = mono_isotopic_peak[i];
    if (keep_only_deisotoped)
    {
      if (z == 0) { continue; }

      // if already single charged or no decharging selected keep peak as it is
      if (!make_single_charged)
      {
        spectra.push_back(old_spectrum[i]);

        // add peak charge to annotation array
        if (annotate_charge)
        {
          spectra.getIntegerDataArrays().back().push_back(z);
        }
      }
      else
      {
        Peak1D p = old_spectrum[i];
        p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
        spectra.push_back(p);

        // add peak charge to annotation array
        if (annotate_charge) { spectra.getIntegerDataArrays().back().push_back(z); }
      }
    }
    else
    {
      // keep all unassigned peaks
      if (features[i] < 0)
      {
        spectra.push_back(old_spectrum[i]);

        // add peak charge to annotation array
        if (annotate_charge) { spectra.getIntegerDataArrays().back().push_back(z); }
        continue;
      }

      // convert mono-isotopic peak with charge assigned by deisotoping
      if (z != 0)
      {
        if (!make_single_charged)
        {
          spectra.push_back(old_spectrum[i]);

          if (annotate_charge) { spectra.getIntegerDataArrays().back().push_back(z); }
        }
        else // make single charged
        {
          Peak1D p = old_spectrum[i];
          p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          spectra.push_back(p);

          // annotate the original charge
          if (annotate_charge)
          {
            spectra.getIntegerDataArrays().back().push_back(z);
          }
        }
      }
    }
  }
  spectra.sortByPosition();
}
} // namespace

