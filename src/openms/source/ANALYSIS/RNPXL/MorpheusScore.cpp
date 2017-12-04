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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/MorpheusScore.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <cmath>

namespace OpenMS
{
  MorpheusScore::Result MorpheusScore::compute(double fragment_mass_tolerance, 
                                bool fragment_mass_tolerance_unit_ppm, 
                                const PeakSpectrum& exp_spectrum, 
                                const PeakSpectrum& theo_spectrum)
  {
    const Size n_t(theo_spectrum.size());
    const Size n_e(exp_spectrum.size());

    MorpheusScore::Result psm = MorpheusScore::Result{0, 0, 0, 0, 0};

    if (n_t == 0 || n_e == 0) { return psm; }

    Size t(0), e(0), matches(0);
    double total_intensity(0);

    // count matching peaks and make sure that every theoretical peak is matched at most once
    while (t < n_t && e < n_e)
    {
      const double theo_mz = theo_spectrum[t].getMZ();
      const double exp_mz = exp_spectrum[e].getMZ();
      const double d = exp_mz - theo_mz;
      const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
      if (fabs(d) <= max_dist_dalton) // match in tolerance window? 
      {
        ++matches;
        ++t;  // count theoretical peak only once
      }
      else if (d < 0) // exp. peak is left of theo. peak (outside of tolerance window)
      {
        total_intensity += exp_spectrum[e].getIntensity();
        ++e; 
      }
      else if (d > 0) // theo. peak is left of exp. peak (outside of tolerance window)
      {
        ++t;
      }
    }

    for (; e < n_e; ++e) { total_intensity += exp_spectrum[e].getIntensity(); }

    // similar to above but we now make sure that the intensity of every matched experimental peak is summed up to form match_intensity
    t = 0; 
    e = 0;
    double match_intensity(0);

    while (t < n_t && e < n_e)
    {
      const double theo_mz = theo_spectrum[t].getMZ();
      const double exp_mz = exp_spectrum[e].getMZ();
      const double d = exp_mz - theo_mz;
      const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
      if (fabs(d) <= max_dist_dalton) // match in tolerance window? 
      {
        match_intensity += exp_spectrum[e].getIntensity();
        ++e; // sum up experimental peak intensity only once
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

    const double intensity_fraction = match_intensity / total_intensity; 

    psm.score = static_cast<double>(matches) + intensity_fraction;
    psm.n_peaks = theo_spectrum.size();
    psm.matches = matches;
    psm.MIC = match_intensity;
    psm.TIC = total_intensity;
    return psm;
  }
}

