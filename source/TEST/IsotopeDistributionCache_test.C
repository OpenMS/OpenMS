// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>

using namespace OpenMS;

START_TEST(IsotopeDistributionCache, "$Id$")

START_SECTION(IsotopeDistributionCache(DoubleReal max_mass, DoubleReal mass_window_width, DoubleReal intensity_percentage=0, DoubleReal intensity_percentage_optional=0))
  IsotopeDistributionCache c(100, 1);
END_SECTION 

START_SECTION(const TheoreticalIsotopePattern& getIsotopeDistribution(DoubleReal mass) const)
  IsotopeDistributionCache c(1000, 10);
  const IsotopeDistributionCache::TheoreticalIsotopePattern &p(c.getIsotopeDistribution(500));
  TEST_REAL_SIMILAR(p.intensity[0], 1);
  TEST_REAL_SIMILAR(p.intensity[1], 0.2668);
  TEST_REAL_SIMILAR(p.intensity[2], 0.04864);
  TEST_REAL_SIMILAR(p.intensity[3], 0.006652);
  TEST_EQUAL(&p == &c.getIsotopeDistribution(509.9), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(510.0), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(499.9), true);
END_SECTION

END_TEST

