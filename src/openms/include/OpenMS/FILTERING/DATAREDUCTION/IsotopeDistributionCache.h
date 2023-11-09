// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

namespace OpenMS
{
  /**
   * @brief Pre-calculate isotope distributions for interesting mass ranges
   */
  class OPENMS_DLLAPI IsotopeDistributionCache
  {
public:
    typedef FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;

    IsotopeDistributionCache(double max_mass, double mass_window_width, double intensity_percentage = 0, double intensity_percentage_optional = 0);

    /// Returns the isotope distribution for a certain mass window
    const TheoreticalIsotopePattern & getIsotopeDistribution(double mass) const;

private:
    /// Vector of pre-calculated isotope distributions for several mass windows
    std::vector<TheoreticalIsotopePattern> isotope_distributions_;

    double mass_window_width_;
  };
}

