// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

double targetFunction(double x)
{
  return 10 + 20*x + 40*x*x;
}

START_TEST(LowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LowessSmoothing* ptr = nullptr;
LowessSmoothing* null_ptr = nullptr;
START_SECTION(LowessSmoothing())
{
	ptr = new LowessSmoothing();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~LowessSmoothing())
{
	delete ptr;
}
END_SECTION

/////
std::vector<double> x, y, y_noisy, out;

//exact data
for (double i = 1.0; i <= 20.0; i += 1.0)
{
    x.push_back(i);
    y.push_back(targetFunction(i));
}

//noisy data
// make some noise

boost::random::mt19937 rnd_gen_;
for (Size i = 0; i < y.size(); ++i)
{
  boost::normal_distribution<float> udist (y.at(i), 0.05);
  y_noisy.push_back( udist(rnd_gen_) );
}

LowessSmoothing lowsmooth;
Param lowpar;
lowpar.setValue("window_size", 15);

TOLERANCE_RELATIVE(1.0004);
TOLERANCE_ABSOLUTE(0.07);
START_SECTION(void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
    y.push_back(-1.0);
    TEST_EXCEPTION(Exception::InvalidValue, lowsmooth.smoothData(x, y, out));
    y.pop_back();
    out.clear();

    lowsmooth.smoothData(x, y, out);
    Size idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

    out.clear();
    lowsmooth.setParameters(lowpar);
    lowsmooth.smoothData(x, y, out);
    idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

    out.clear();
    lowsmooth.smoothData(x, y_noisy, out);
    idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



