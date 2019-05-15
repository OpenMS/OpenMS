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
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>

using std::vector;

namespace OpenMS
{
  inline double HyperScore::logfactorial_(const int x, int base)
  {
    double z(0);
    base = std::max(base, 2);
    for (int i = base; i <= x; ++i)
    {
      z += log(i);
    }
    return z;
  }


  double HyperScore::compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const PeakSpectrum& theo_spectrum)
  {
    if (exp_spectrum.size() < 1 || theo_spectrum.size() < 1)
    {
      std::cout << "Warning: HyperScore: One of the given spectra is empty." << std::endl;
      return 0.0;
    }

    // TODO this assumes only one StringDataArray is present and it is the right one
    const PeakSpectrum::StringDataArray* ion_names;
    if (theo_spectrum.getStringDataArrays().size() > 0)
    {
      ion_names = &theo_spectrum.getStringDataArrays()[0];
    }
    else
    {
      std::cout << "Error: HyperScore: Theoretical spectrum without StringDataArray (\"IonNames\" annotation) provided." << std::endl;
      return 0.0;
    }

    int y_ion_count = 0;
    int b_ion_count = 0;
    double dot_product = 0.0;
    if (fragment_mass_tolerance_unit_ppm) 
    {
      MatchedIterator<PeakSpectrum, PpmTrait, true> it(theo_spectrum, exp_spectrum, fragment_mass_tolerance);
      for (; it != it.end(); ++it)
      {
        dot_product += (*it).getIntensity() * it.ref().getIntensity(); /* * mass_error */;
        // fragment annotations in XL-MS data are more complex and do not start with the ion type, but the ion type always follows after a $
        auto i = it.refIdx();
        if ((*ion_names)[i][0] == 'y' || (*ion_names)[i].hasSubstring("$y"))
        {
          ++y_ion_count;
        }
        else if ((*ion_names)[i][0] == 'b' || (*ion_names)[i].hasSubstring("$b"))
        {
          ++b_ion_count;
        }
      }
    }
    else
    {
      MatchedIterator<PeakSpectrum, DaTrait, true> it(theo_spectrum, exp_spectrum, fragment_mass_tolerance);
      for (; it != it.end(); ++it)
      {
        dot_product += (*it).getIntensity() * it.ref().getIntensity(); /* * mass_error */;
        // fragment annotations in XL-MS data are more complex and do not start with the ion type, but the ion type always follows after a $
        auto i = it.refIdx();
        if ((*ion_names)[i][0] == 'y' || (*ion_names)[i].hasSubstring("$y"))
        {
          ++y_ion_count;
        }
        else if ((*ion_names)[i][0] == 'b' || (*ion_names)[i].hasSubstring("$b"))
        {
          ++b_ion_count;
        }
      }

    }

    // inefficient: calculates logs repeatedly
    //const double yFact = logfactorial_(y_ion_count);
    //const double bFact = logfactorial_(b_ion_count);
    //const double hyperScore = log1p(dot_product) + yFact + bFact;

    const int i_min = std::min(y_ion_count, b_ion_count);
    const int i_max = std::max(y_ion_count, b_ion_count);
    const double hyperScore = log1p(dot_product) + 2*logfactorial_(i_min) + logfactorial_(i_max, i_min + 1);
    return hyperScore;
  }

}

