// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/HyperScore.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>


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
    if (exp_spectrum.empty() || theo_spectrum.empty())
    {
      std::cout << "Warning: HyperScore: One of the given spectra is empty." << std::endl;
      return 0.0;
    }

    // TODO this assumes only one StringDataArray is present and it is the right one
    const PeakSpectrum::StringDataArray* ion_names;
    if (!theo_spectrum.getStringDataArrays().empty())
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

  double HyperScore::computeWithDetail(double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const PeakSpectrum& exp_spectrum, 
    const PeakSpectrum& theo_spectrum,
    PSMDetail& d)
  {
    if (exp_spectrum.empty() || theo_spectrum.empty())
    {
      std::cout << "Warning: HyperScore: One of the given spectra is empty." << std::endl;
      return 0.0;
    }

    // TODO this assumes only one StringDataArray is present and it is the right one
    const PeakSpectrum::StringDataArray* ion_names;
    if (!theo_spectrum.getStringDataArrays().empty())
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
    double abs_error = 0.0;
    if (fragment_mass_tolerance_unit_ppm) 
    {
      MatchedIterator<PeakSpectrum, PpmTrait, true> it(theo_spectrum, exp_spectrum, fragment_mass_tolerance);
      for (; it != it.end(); ++it)
      {
        const double exp_mz{it.ref().getMZ()};
        const double theo_mz{(*it).getMZ()};
        const double exp_int{it.ref().getIntensity()};
        const double theo_int{(*it).getIntensity()};
        abs_error += Math::getPPMAbs(theo_mz, exp_mz);
        dot_product += theo_int * exp_int; /* * mass_error */;
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
        abs_error += abs((*it).getMZ() - it.ref().getMZ());
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
    d.matched_b_ions = b_ion_count;
    d.matched_y_ions = y_ion_count;
    d.mean_error = (b_ion_count + y_ion_count) > 0 ? abs_error / (double)(b_ion_count + y_ion_count) : 0.0;
    return hyperScore;
  }

}

