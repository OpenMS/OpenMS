// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GammaDistributionFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GammaDistributionFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GammaDistributionFitter* ptr = 0;
CHECK(GammaDistributionFitter())
{
	ptr = new GammaDistributionFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(virtual ~GammaDistributionFitter())
{
	delete ptr;
}
RESULT

CHECK((GammaDistributionFitter(const GammaDistributionFitter & rhs)))
{
  GammaDistributionFitter::GammaDistributionFitResult result;
	result.b = 0.3;
	result.p = 0.7;
	GammaDistributionFitter f1;
	f1.setInitialParameters(result);
	GammaDistributionFitter f2(f1);
	TEST_REAL_EQUAL(f1.getInitialParameters().b, result.b)
	TEST_REAL_EQUAL(f1.getInitialParameters().p, result.p)
	TEST_REAL_EQUAL(f1.getInitialParameters().b, f2.getInitialParameters().b)
	TEST_REAL_EQUAL(f1.getInitialParameters().p, f2.getInitialParameters().p)
}
RESULT

CHECK((GammaDistributionFitter& operator=(const GammaDistributionFitter & rhs)))
{
  GammaDistributionFitter::GammaDistributionFitResult result;
	result.b = 0.3;
	result.p = 0.7;
	GammaDistributionFitter f1;
	f1.setInitialParameters(result);
	GammaDistributionFitter f2;
	f2 = f1;
	TEST_EQUAL(f1.getInitialParameters().b, result.b)
	TEST_EQUAL(f1.getInitialParameters().p, result.p)
	TEST_EQUAL(f2.getInitialParameters().b, result.b)
	TEST_EQUAL(f2.getInitialParameters().p, result.p)
	TEST_EQUAL(f1.getInitialParameters().b, f2.getInitialParameters().b)
	TEST_EQUAL(f1.getInitialParameters().p, f2.getInitialParameters().p)
}
RESULT

CHECK((GammaDistributionFitResult fit(std::vector< DPosition< 2 > > & points)))
{
	/*
	DPosition<2> pos;
  pos.setX(0.0);
  pos.setY(0.01);
  vector<DPosition<2> > points;
  points.push_back(pos);
  pos.setX(0.05); pos.setY(0.2*4); points.push_back(pos);
	pos.setX(0.08); pos.setY(0.3*4); points.push_back(pos);
	pos.setX(0.12); pos.setY(0.5*4); points.push_back(pos);
  pos.setX(0.16); pos.setY(0.63*4); points.push_back(pos);
  pos.setX(0.28); pos.setY(0.99*4); points.push_back(pos);
	pos.setX(0.43); pos.setY(0.9*4); points.push_back(pos);
  pos.setX(0.66); pos.setY(0.83*4); points.push_back(pos);
  pos.setX(0.78); pos.setY(0.10*4); points.push_back(pos);
	pos.setX(0.99); pos.setY(0.02*4); points.push_back(pos);
	

  ptr = new GammaDistributionFitter;
  GammaDistributionFitter::GammaDistributionFitResult result = ptr->fit(points);

  PRECISION(0.1)
  TEST_REAL_EQUAL(result.b, 1.0)
  TEST_REAL_EQUAL(result.p, 0.3)
	*/
}
RESULT

CHECK((const GammaDistributionFitResult& getInitialParameters() const ))
{
  NOT_TESTABLE // tested above
}
RESULT

CHECK((void setInitialParameters(const GammaDistributionFitResult & result)))
{
  GammaDistributionFitter f1;
	GammaDistributionFitter::GammaDistributionFitResult result = f1.getInitialParameters();
	TEST_REAL_EQUAL(result.b, 1.0)
	TEST_REAL_EQUAL(result.p, 5.0)
	result.b = 0.15;
	result.p = 0.24;
	f1.setInitialParameters(result);

	TEST_REAL_EQUAL(f1.getInitialParameters().b, 0.15)
	TEST_REAL_EQUAL(f1.getInitialParameters().p, 0.24)
}
RESULT

CHECK((const String& getGnuplotFormula() const ))
	// TODO
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



