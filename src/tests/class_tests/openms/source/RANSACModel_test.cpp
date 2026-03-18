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
#include <OpenMS/ML/RANSAC/RANSACModel.h>
///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(RANSACModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((ModelParameters rm_fit(const DVecIt& begin, const DVecIt& end) const))
  NOT_TESTABLE // since base class; test derived classes
END_SECTION

START_SECTION((double rm_rsq(const DVecIt& begin, const DVecIt& end) const))
  NOT_TESTABLE // since base class; test derived classes
END_SECTION

START_SECTION((double rm_rss(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients) const))
  NOT_TESTABLE // since base class; test derived classes
END_SECTION

START_SECTION((DVec rm_inliers(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold) const))
  NOT_TESTABLE // since base class; test derived classes
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

