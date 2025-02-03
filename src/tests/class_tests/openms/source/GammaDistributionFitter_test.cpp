// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GammaDistributionFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GammaDistributionFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GammaDistributionFitter* ptr = nullptr;
GammaDistributionFitter* nullPointer = nullptr;
START_SECTION(GammaDistributionFitter())
{
	ptr = new GammaDistributionFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~GammaDistributionFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((GammaDistributionFitResult fit(std::vector< DPosition< 2 > > & points)))
{
	DPosition<2> pos;
  vector<DPosition<2> > points;

	pos.setX(0.0001); pos.setY(0.1); points.push_back(pos);
	pos.setX(0.0251); pos.setY(0.3); points.push_back(pos);
	pos.setX(0.0501); pos.setY(0); points.push_back(pos);
	pos.setX(0.0751); pos.setY(0.7); points.push_back(pos);
	pos.setX(0.1001); pos.setY(0); points.push_back(pos);
	pos.setX(0.1251); pos.setY(1.6); points.push_back(pos);
	pos.setX(0.1501); pos.setY(0); points.push_back(pos);
	pos.setX(0.1751); pos.setY(2.1); points.push_back(pos);
	pos.setX(0.2001); pos.setY(0); points.push_back(pos);
	pos.setX(0.2251); pos.setY(3.7); points.push_back(pos);
	pos.setX(0.2501); pos.setY(0); points.push_back(pos);
	pos.setX(0.2751); pos.setY(4); points.push_back(pos);
	pos.setX(0.3001); pos.setY(0); points.push_back(pos);
	pos.setX(0.3251); pos.setY(3); points.push_back(pos);
	pos.setX(0.3501); pos.setY(0); points.push_back(pos);
	pos.setX(0.3751); pos.setY(2.6); points.push_back(pos);
	pos.setX(0.4001); pos.setY(0); points.push_back(pos);
	pos.setX(0.4251); pos.setY(3); points.push_back(pos);
	pos.setX(0.4501); pos.setY(0); points.push_back(pos);
	pos.setX(0.4751); pos.setY(3); points.push_back(pos);
	pos.setX(0.5001); pos.setY(0); points.push_back(pos);
	pos.setX(0.5251); pos.setY(2.5); points.push_back(pos);
	pos.setX(0.5501); pos.setY(0); points.push_back(pos);
	pos.setX(0.5751); pos.setY(1.7); points.push_back(pos);
	pos.setX(0.6001); pos.setY(0); points.push_back(pos);
	pos.setX(0.6251); pos.setY(1); points.push_back(pos);
	pos.setX(0.6501); pos.setY(0); points.push_back(pos);
	pos.setX(0.6751); pos.setY(0.5); points.push_back(pos);
	pos.setX(0.7001); pos.setY(0); points.push_back(pos);
	pos.setX(0.7251); pos.setY(0.3); points.push_back(pos);
	pos.setX(0.7501); pos.setY(0); points.push_back(pos);
	pos.setX(0.7751); pos.setY(0.4); points.push_back(pos);
	pos.setX(0.8001); pos.setY(0); points.push_back(pos);
	pos.setX(0.8251); pos.setY(0); points.push_back(pos);
	pos.setX(0.8501); pos.setY(0); points.push_back(pos);
	pos.setX(0.8751); pos.setY(0.1); points.push_back(pos);
	pos.setX(0.9001); pos.setY(0); points.push_back(pos);
	pos.setX(0.9251); pos.setY(0.1); points.push_back(pos);
	pos.setX(0.9501); pos.setY(0); points.push_back(pos);
	pos.setX(0.9751); pos.setY(0.2); points.push_back(pos);

	GammaDistributionFitter::GammaDistributionFitResult init_param(1.0, 3.0);

  ptr = new GammaDistributionFitter;
	ptr->setInitialParameters(init_param);
  GammaDistributionFitter::GammaDistributionFitResult result = ptr->fit(points);

  TOLERANCE_ABSOLUTE(0.01)
  TEST_REAL_SIMILAR(result.b, 7.25)
  TEST_REAL_SIMILAR(result.p, 3.11)
}
END_SECTION

START_SECTION((void setInitialParameters(const GammaDistributionFitResult & result)))
{
  GammaDistributionFitter f1;
  GammaDistributionFitter::GammaDistributionFitResult result (1.0, 5.0);
  f1.setInitialParameters(result);

	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
