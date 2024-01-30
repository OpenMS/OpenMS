// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  double HyperScore::compute(double fragment_mass_tolerance, 
                             bool fragment_mass_tolerance_unit_ppm, 
                             const PeakSpectrum& exp_spectrum, 
                             const DataArrays::IntegerDataArray& exp_charges,
                             const PeakSpectrum& theo_spectrum,
                             const DataArrays::IntegerDataArray& theo_charges)
  {
    double dot_product = 0.0;
    UInt y_ion_count = 0;
    UInt b_ion_count = 0;

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

    if (theo_charges.size() != theo_spectrum.size())
    {
      std::cout << "Error: HyperScore: #charges != #peaks in theoretical spectrum." << std::endl;
      return 0.0;
    }

    if (exp_charges.size() != exp_spectrum.size())
    {
      std::cout << "Error: HyperScore: #charges != #peaks in experimental spectrum." << std::endl;
      return 0.0;
    }

    for (Size i = 0; i < theo_spectrum.size(); ++i)
    {
      const double theo_mz = theo_spectrum[i].getMZ();

      double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
      Size index = exp_spectrum.findNearest(theo_mz);

      const double exp_mz = exp_spectrum[index].getMZ();
      const double theo_intensity = theo_spectrum[i].getIntensity();
      const int exp_z = exp_charges[index];
      const int theo_z = theo_charges[i];

      #ifdef DEBUG_HYPERSCORE
        if (exp_z != theo_z) std::cout << "exp_z != theo_z " << exp_z << "\t" << theo_z << std::endl;
      #endif

      // found peak match
      if (std::abs(theo_mz - exp_mz) < max_dist_dalton 
          && exp_z == theo_z)
      {
        // fragment annotations in XL-MS data are more complex and do not start with the ion type, but the ion type always follows after a $
        if ((*ion_names)[i][0] == 'y' || (*ion_names)[i].hasSubstring("$y"))
        {
          dot_product += exp_spectrum[index].getIntensity() * theo_intensity;
          #ifdef DEBUG_HYPERSCORE
            std::cout << (*ion_names)[i] << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          #endif
          ++y_ion_count;
        }
        else if ((*ion_names)[i][0] == 'b' || (*ion_names)[i].hasSubstring("$b"))
        {
          dot_product += exp_spectrum[index].getIntensity() * theo_intensity;

          #ifdef DEBUG_HYPERSCORE
            std::cout << (*ion_names)[i] << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
          #endif
          ++b_ion_count;
        }
      }
    }

    if (y_ion_count == 0 && b_ion_count == 0) return 0;

    const double yFact = logfactorial_(y_ion_count);
    const double bFact = logfactorial_(b_ion_count);
    const double hyperScore = log1p(dot_product) + yFact + bFact;
    #ifdef DEBUG_HYPERSCORE
      std::cout << "HyperScore/#y/#b: " << hyperScore << "/" << y_ion_count << "/" << b_ion_count << std::endl;
    #endif
    return hyperScore;    
  } 


  double HyperScore::compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        const PeakSpectrum& theo_spectrum,
                        const DataArrays::IntegerDataArray& theo_charges,
                        std::vector<double>& intensity_sum)
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

    if (theo_charges.size() != theo_spectrum.size())
    {
      std::cout << "Error: HyperScore: #charges != #peaks in theoretical spectrum." << std::endl;
      return 0.0;
    }

    if (exp_charges.size() != exp_spectrum.size())
    {
      std::cout << "Error: HyperScore: #charges != #peaks in experimental spectrum." << std::endl;
      return 0.0;
    }

    double dot_product = 0.0;
    const Size N = intensity_sum.size(); // length of peptide
    std::vector<double> b_ions(N, 0.0);
    std::vector<double> y_ions(N, 0.0);

    for (Size i = 0; i < theo_spectrum.size(); ++i)
    {
      const double theo_mz = theo_spectrum[i].getMZ();

      double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
      Size index = exp_spectrum.findNearest(theo_mz);

      const double exp_mz = exp_spectrum[index].getMZ();
      const double theo_intensity = theo_spectrum[i].getIntensity();
      const int exp_z = exp_charges[index];
      const int theo_z = theo_charges[i];

      #ifdef DEBUG_HYPERSCORE
        if (exp_z != theo_z) std::cout << "exp_z != theo_z " << exp_z << "\t" << theo_z << std::endl;
      #endif

      // found peak match
      if (std::abs(theo_mz - exp_mz) < max_dist_dalton 
          && exp_z == theo_z)
      {
        // fragment annotations in XL-MS data are more complex and do not start with the ion type, but the ion type always follows after a $
        if ((*ion_names)[i][0] == 'y')
        {
          auto start = (*ion_names)[i].begin() + 1;
          auto end = (*ion_names)[i].end();
          int ii(0);
          if (!StringUtils::extractInt(start, end, ii)) continue;
          dot_product += exp_spectrum[index].getIntensity() * theo_intensity;


          #ifdef DEBUG_HYPERSCORE
            std::cout << (*ion_names)[i] << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
            std::cout << "N:"  << N << "\t" << ii << "\t" << N-1 - (ii-1) << std::endl;            
          #endif
          // we observed the suffix (N-1-ii, N-1] in 0-based AA coordinates
          y_ions[N-1 - (ii-1)] += exp_spectrum[index].getIntensity();
        }
        else if ((*ion_names)[i][0] == 'b')
        {
          auto start = (*ion_names)[i].begin() + 1;
          auto end = (*ion_names)[i].end();
          int ii(0);
          if (!StringUtils::extractInt(start, end, ii)) continue;
          dot_product += exp_spectrum[index].getIntensity() * theo_intensity;

          #ifdef DEBUG_HYPERSCORE
            std::cout << (*ion_names)[i] << " intensity: " << exp_spectrum[index].getIntensity() << std::endl;
            std::cout << "N:"  << N << "\t" << ii << "\t" << (ii-1) << std::endl;
          #endif
          // we observed the prefix [0, ii) in 0-based AA coordinates
          b_ions[ii - 1] += exp_spectrum[index].getIntensity();
        }
      }
    }

    UInt y_ion_count(0), 
      b_ion_count(0);
    for (Size i = 0; i != b_ions.size(); ++i) 
    {
      if (b_ions[i] > 0) 
      {
        ++b_ion_count;
        intensity_sum[i] += b_ions[i];
      }       
    } 

    for (Size i = 0; i != y_ions.size(); ++i) 
    {
      if (y_ions[i] > 0) 
      {
        ++y_ion_count;
        intensity_sum[i] += y_ions[i];
      }       
    }

    if (y_ion_count == 0 && b_ion_count == 0) return 0;

    const double yFact = logfactorial_(y_ion_count);
    const double bFact = logfactorial_(b_ion_count);
    const double hyperScore = log1p(dot_product) + yFact + bFact;
    #ifdef DEBUG_HYPERSCORE
      std::cout << "HyperScore/#y/#b: " << hyperScore << "/" << y_ion_count << "/" << b_ion_count << std::endl;
    #endif
    return hyperScore;    
  } 

}

