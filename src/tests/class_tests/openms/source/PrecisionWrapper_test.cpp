// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <sstream>

///////////////////////////

START_TEST(PrecisionWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((template <typename FloatingPointType> explicit PrecisionWrapper(const FloatingPointType rhs)))
{
  PrecisionWrapper<double> w(3.5);
  TEST_REAL_SIMILAR(w.ref_, 3.5)
}
END_SECTION

START_SECTION((PrecisionWrapper(const PrecisionWrapper& rhs)))
{
  PrecisionWrapper<double> w(2.25);
  PrecisionWrapper<double> copy(w);
  TEST_REAL_SIMILAR(copy.ref_, 2.25)
}
END_SECTION

START_SECTION((template <typename FloatingPointType> const PrecisionWrapper<FloatingPointType> precisionWrapper(const FloatingPointType rhs)))
{
  // the make-function deduces the type and stores the value
  auto wd = precisionWrapper(1.5);
  TEST_REAL_SIMILAR(wd.ref_, 1.5)
  auto wf = precisionWrapper(0.5f);
  TEST_REAL_SIMILAR(wf.ref_, 0.5f)
}
END_SECTION

START_SECTION((template <typename FloatingPointType> std::ostream& operator<<(std::ostream& os, const PrecisionWrapper<FloatingPointType>& rhs)))
{
  // full-precision output for double
  {
    std::ostringstream os;
    os << precisionWrapper(0.1234567890123456789);
    TEST_EQUAL(os.str(), "0.123456789012346")
  }
  // full-precision output for float
  {
    std::ostringstream os;
    os << precisionWrapper(0.1234567890123456789f);
    TEST_EQUAL(os.str(), "0.123457")
  }
  // round numbers keep a decimal point (StringUtils::toStr full-precision form)
  {
    std::ostringstream os;
    os << precisionWrapper(1.0);
    TEST_EQUAL(os.str(), "1.0")
  }
  {
    std::ostringstream os;
    os << precisionWrapper(100.5);
    TEST_EQUAL(os.str(), "100.5")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
