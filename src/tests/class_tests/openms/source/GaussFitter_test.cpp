// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GaussFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter* ptr = nullptr;
GaussFitter* nullPointer = nullptr;
START_SECTION(GaussFitter())
{
	ptr = new GaussFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~GaussFitter()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

double mz[] = {
  240.1000470172,
  240.1002675493,
  240.1004880817,
  240.1007086145,
  240.1009291475,
  240.1011496808,
  240.1013702145};

double ints[] = {
  61134.39453125,
  111288.5390625,
  163761.46875,
  165861.4375,
  162133.46875,
  120060.5234375,
  71102.1328125,
};

// initial guesses
double max_peak_int = 168324;
double max_peak_mz =  240.10051;
double sigma = 0.000375375;

Math::GaussFitter::GaussFitResult gfi(max_peak_int, max_peak_mz, sigma);


START_SECTION((GaussFitResult fit(std::vector< DPosition< 2 > >& points) const))
{
  DPosition<2> pos;
	pos.setX(0.0);
	pos.setY(0.01);
	vector<DPosition<2> > points;
	points.push_back(pos);
	pos.setX(0.05);
	pos.setY(0.2);
	points.push_back(pos);
	pos.setX(0.16);
	pos.setY(0.63);
	points.push_back(pos);
	pos.setX(0.28);
	pos.setY(0.99);
	points.push_back(pos);
	pos.setX(0.66);
	pos.setY(0.03);
	points.push_back(pos);
	pos.setX(0.50);
	pos.setY(0.36);
	points.push_back(pos);
	
	ptr = new GaussFitter;
	GaussFitter::GaussFitResult result = ptr->fit(points);

	//TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result.A, 1.01898275662372)
	TEST_REAL_SIMILAR(result.x0, 0.300612870901173)
	TEST_REAL_SIMILAR(result.sigma, 0.136316330927453)

  ////////////////////////////////////////
  // second case which results in a negative sigma internally (requires using fabs())
  ////////////////////////////////////////
  std::vector< DPosition< 2 > > gp;
  for (Size iii = 0; iii < 7; ++iii)
  {
    DPosition< 2 > d(mz[iii], ints[iii]);
    gp.push_back(d);
  }

  Math::GaussFitter gf;
  gf.setInitialParameters(gfi);
  Math::GaussFitter::GaussFitResult gfr = gf.fit(gp);
  /*
  x0:      240.10051 --> 240.1007246725147
  sigma: 0.000375375 --> 0.00046642320683761701
  A:          168324 --> 175011.8930067491
  */
  TEST_REAL_SIMILAR(gfr.A, 175011.893006749)
  TEST_REAL_SIMILAR(gfr.x0, 240.1007246725147)
  TEST_REAL_SIMILAR(gfr.sigma, 0.00046642320683761701)
  delete ptr;
}
END_SECTION

START_SECTION((void setInitialParameters(const GaussFitResult& result)))
{
  GaussFitter f1;
  GaussFitter::GaussFitResult result (-1,-1,-1);
  f1.setInitialParameters(result);

	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((static std::vector<double> eval(const std::vector<double>& evaluation_points, const GaussFitResult& model)))
  GaussFitter f1;
  std::vector<double> rnd = Math::GaussFitter::eval(std::vector<double>(&mz[0], &mz[0] + 7), gfi);

  double int_fitted[] = {
    78670.515322697669,
    136633.77791868619,
    168037.29915800504,
    146337.00743127937,
    90240.802825824489,
    39405.008909696895,
    12184.248044493703
  };
  for (Size i=0; i < rnd.size(); ++i)
  {
    TEST_REAL_SIMILAR(int_fitted[i], rnd[i])
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



