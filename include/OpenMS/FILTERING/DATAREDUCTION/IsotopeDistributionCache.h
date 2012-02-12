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

#ifndef OPENMS_FILTERING_DATAREDUCTION_ISOTOPEDISTRIBUTIONCACHE_H
#define OPENMS_FILTERING_DATAREDUCTION_ISOTOPEDISTRIBUTIONCACHE_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

namespace OpenMS
{
  /**
   * @brief Prealculate isotope distributions for interesting mass ranges
   */
  class OPENMS_DLLAPI IsotopeDistributionCache
  {
    public:
      typedef FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;

      IsotopeDistributionCache(DoubleReal max_mass, DoubleReal mass_window_width, DoubleReal intensity_percentage = 0, DoubleReal intensity_percentage_optional = 0);

      /// Returns the isotope distribution for a certain mass window
      const TheoreticalIsotopePattern& getIsotopeDistribution(DoubleReal mass) const;

    private:
      /// Vector of precalculated isotope distributions for several mass winows
      std::vector< TheoreticalIsotopePattern > isotope_distributions_;

      DoubleReal mass_window_width_;
  };
}

#endif
