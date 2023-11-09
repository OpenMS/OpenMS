// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>

using namespace OpenMS;

START_TEST(IsotopeDistributionCache, "$Id$")

START_SECTION(IsotopeDistributionCache(double max_mass, double mass_window_width, double intensity_percentage=0, double intensity_percentage_optional=0))
  IsotopeDistributionCache c(100, 1);
END_SECTION 

START_SECTION(const TheoreticalIsotopePattern& getIsotopeDistribution(double mass) const)
  IsotopeDistributionCache c(1000, 10);
  const IsotopeDistributionCache::TheoreticalIsotopePattern &p(c.getIsotopeDistribution(500));
  TEST_REAL_SIMILAR(p.intensity[0], 1);
  TEST_REAL_SIMILAR(p.intensity[1], 0.267834);
  TEST_REAL_SIMILAR(p.intensity[2], 0.048924);
  TEST_REAL_SIMILAR(p.intensity[3], 0.006703);
  TEST_EQUAL(&p == &c.getIsotopeDistribution(509.9), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(510.0), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(499.9), true);
END_SECTION

END_TEST

