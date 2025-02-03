// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/MorpheusScore.h>
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

    MorpheusScore::Result psm = {};

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
    double match_intensity(0.0);
    double sum_error(0.0);

    while (t < n_t && e < n_e)
    {
      const double theo_mz = theo_spectrum[t].getMZ();
      const double exp_mz = exp_spectrum[e].getMZ();
      const double d = exp_mz - theo_mz;
      const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
      if (fabs(d) <= max_dist_dalton) // match in tolerance window? 
      {
        match_intensity += exp_spectrum[e].getIntensity();
        sum_error += fabs(d);
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
    psm.err = matches > 0 ? sum_error / static_cast<double>(matches) : 1e10;
    return psm;
  }
}

