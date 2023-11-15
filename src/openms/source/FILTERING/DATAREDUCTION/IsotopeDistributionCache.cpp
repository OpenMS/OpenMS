// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  IsotopeDistributionCache::IsotopeDistributionCache(double max_mass, double mass_window_width, double intensity_percentage, double intensity_percentage_optional) :
    mass_window_width_(mass_window_width)
  {
    Size num_isotopes = std::ceil(max_mass / mass_window_width) + 1;

    //reserve enough space
    isotope_distributions_.resize(num_isotopes);

    //calculate distribution if necessary
    for (Size index = 0; index < num_isotopes; ++index)
    {
      //log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
      CoarseIsotopePatternGenerator solver(20);
      auto d = solver.estimateFromPeptideWeight(0.5 * mass_window_width + index * mass_window_width);

      //trim left and right. And store the number of isotopes on the left, to reconstruct the monoisotopic peak
      Size size_before = d.size();
      d.trimLeft(intensity_percentage_optional);
      isotope_distributions_[index].trimmed_left = size_before - d.size();
      d.trimRight(intensity_percentage_optional);

      for (Peak1D& it : d)
      {
        isotope_distributions_[index].intensity.push_back(it.getIntensity());
        //log_ << " - " << it.second << std::endl;
      }

      //determine the number of optional peaks at the beginning/end
      Size begin = 0;
      Size end = 0;
      bool is_begin = true;
      bool is_end = false;
      for (Size i = 0; i < isotope_distributions_[index].intensity.size(); ++i)
      {
        if (isotope_distributions_[index].intensity[i] < intensity_percentage)
        {
          if (!is_end && !is_begin)
          {
            is_end = true;
          }
          if (is_begin)
          {
            ++begin;
          }
          else if (is_end)
          {
            ++end;
          }
        }
        else if (is_begin)
        {
          is_begin = false;
        }
      }
      isotope_distributions_[index].optional_begin = begin;
      isotope_distributions_[index].optional_end = end;

      //scale the distribution to a maximum of 1
      double max = 0.0;
      for (Size i = 0; i < isotope_distributions_[index].intensity.size(); ++i)
      {
        if (isotope_distributions_[index].intensity[i] > max)
        {
          max = isotope_distributions_[index].intensity[i];
        }
      }

      isotope_distributions_[index].max = max;

      for (Size i = 0; i < isotope_distributions_[index].intensity.size(); ++i)
      {
        isotope_distributions_[index].intensity[i] /= max;
      }
    }
  }

  // Returns the isotope distribution for a certain mass window
  const IsotopeDistributionCache::TheoreticalIsotopePattern& IsotopeDistributionCache::getIsotopeDistribution(double mass) const
  {
    //calculate index in the vector
    Size index = static_cast<Size>(std::floor(mass / mass_window_width_));

    if (index >= isotope_distributions_.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsotopeDistribution not precalculated. Maximum allowed index is " + String(isotope_distributions_.size()), String(index));
    }

    //Return distribution
    return isotope_distributions_[index];
  }

}
