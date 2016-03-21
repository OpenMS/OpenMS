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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>

#include <vector>
#include <map>

using std::vector;

namespace OpenMS
{
  double HyperScore::logfactorial_(UInt x)
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
      const double& theo_intensity = theo_peak_it->getIntensity();

      double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
      Size index = exp_spectrum.findNearest(theo_mz);
      double exp_mz = exp_spectrum[index].getMZ();

      // found peak match
      if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
      {
        dot_product += exp_spectrum[index].getIntensity() * theo_intensity;
        if (theo_peak_it->getMetaValue("IonName").toString()[0] == 'y')
        {
          #ifdef DEBUG_HYPERSCORE
            std::cout << theo_peak_it->getMetaValue("IonName").toString() << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          #endif
          ++y_ion_count;
        }
        else if (theo_peak_it->getMetaValue("IonName").toString()[0] == 'b')
        {
          #ifdef DEBUG_HYPERSCORE
            std::cout << theo_peak_it->getMetaValue("IonName").toString() << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          #endif
          ++b_ion_count;
        }
      }
    }

    // discard very low scoring hits (basically no matching peaks)
    if (dot_product > 1e-1)
    {
      double yFact = logfactorial_(y_ion_count);
      double bFact = logfactorial_(b_ion_count);
      double hyperScore = log(dot_product) + yFact + bFact;
      return hyperScore;
    }
    else
    {
      return 0;
    }
  }

}

