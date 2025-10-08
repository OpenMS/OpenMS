// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/RANSAC/RANSAC.h>
///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(RANSAC, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION((static std::vector<std::pair<double, double> > ransac(const std::vector<std::pair<double, double> >& pairs, 
                                                                     size_t n, 
                                                                     size_t k, 
                                                                     double t, 
                                                                     size_t d, 
                                                                     int (*rng)(int) = NULL)))
{
  NOT_TESTABLE // tested in the models, e.g. RANSACModelLinear
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

