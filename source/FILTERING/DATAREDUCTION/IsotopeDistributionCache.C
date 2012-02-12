// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>

using namespace OpenMS;

IsotopeDistributionCache::IsotopeDistributionCache(DoubleReal max_mass, DoubleReal mass_window_width, DoubleReal intensity_percentage, DoubleReal intensity_percentage_optional)
  : mass_window_width_(mass_window_width)
{
  Size num_isotopes = std::ceil(max_mass/mass_window_width) + 1;

  //reserve enough space
  isotope_distributions_.resize(num_isotopes);

  //calculate distribution if necessary
  for (Size index = 0; index < num_isotopes; ++index)
  {
    //log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
    IsotopeDistribution d;
    d.setMaxIsotope(20);
    d.estimateFromPeptideWeight(0.5*mass_window_width + index * mass_window_width);

    //trim left and right. And store the number of isotopes on the left, to reconstruct the monoisotopic peak
    Size size_before = d.size();
    d.trimLeft(intensity_percentage_optional);
    isotope_distributions_[index].trimmed_left = size_before - d.size();
    d.trimRight(intensity_percentage_optional);

    for (IsotopeDistribution::Iterator it=d.begin(); it!=d.end(); ++it)
    {
      isotope_distributions_[index].intensity.push_back(it->second);
      //log_ << " - " << it->second << std::endl;
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
        if (!is_end && !is_begin) is_end = true;
        if (is_begin) ++begin;
        else if (is_end) ++end;
      }
      else if (is_begin)
      {
        is_begin = false;
      }
    }
    isotope_distributions_[index].optional_begin = begin;
    isotope_distributions_[index].optional_end = end;

    //scale the distibution to a maximum of 1
    DoubleReal max = 0.0;
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

/// Returns the isotope distribution for a certain mass window
const IsotopeDistributionCache::TheoreticalIsotopePattern& IsotopeDistributionCache::getIsotopeDistribution(DoubleReal mass) const
{
  //calculate index in the vector
  Size index = (Size)std::floor(mass / mass_window_width_);

  if (index >= isotope_distributions_.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsotopeDistribution not precalculated. Maximum allowed index is " + String(isotope_distributions_.size()), String(index));
  }

  //Return distribution
  return isotope_distributions_[index];
}
