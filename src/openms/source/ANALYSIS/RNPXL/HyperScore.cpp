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
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>

#include <vector>
#include <map>

using std::vector;

namespace OpenMS
{
  double HyperScore::logfactorial(UInt x)
  {
    UInt y;

    if (x < 2)
      return 1;
    else
    {
      double z = 0;
      for (y = 2; y <= x; y++)
      {
        z = log((double)y) + z;
      }

      return z;
    }
  }

  double HyperScore::compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const RichPeakSpectrum& theo_spectrum)
  {
    double dot_product = 0.0;
    UInt y_ion_count = 0;
    UInt b_ion_count = 0;

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
        dot_product += exp_spectrum[index].getIntensity();
        if (theo_peak_it->getMetaValue("IonName").toString()[0] == 'y')
        {
          std::cout << theo_peak_it->getMetaValue("IonName").toString() << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          ++y_ion_count;
        }
        else if (theo_peak_it->getMetaValue("IonName").toString()[0] == 'b')
        {
          std::cout << theo_peak_it->getMetaValue("IonName").toString() << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          ++b_ion_count;
        }
      }
    }

    // discard very low scoring hits (basically no matching peaks)
    if (dot_product > 1e-1)
    {
      double yFact = logfactorial(y_ion_count);
      double bFact = logfactorial(b_ion_count);
      double hyperScore = log(dot_product) + yFact + bFact;
      return hyperScore;
    }
    else
    {
      return 0;
    }
  }

  HyperScore::IndexScorePair HyperScore::compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const vector<RichPeakSpectrum>& theo_spectra)
  {
    double best_score(0.0);
    Size best_index(0);

    for (Size i = 0; i != theo_spectra.size(); ++i)
    {
      const RichPeakSpectrum& theo_spectrum = theo_spectra[i];
      const double score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);
      if (score > best_score)
      {
        best_score = score;
        best_index = i;
      }
    }
    return std::make_pair(best_index, best_score);
  }
}

