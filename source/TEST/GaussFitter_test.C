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
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GaussFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter* ptr = 0;
CHECK(GaussFitter())
{
	ptr = new GaussFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK((GaussFitter(const GaussFitter & rhs)))
{
  GaussFitter::GaussFitResult result;
	result.A = 0.3;
	result.x0 = 0.4;
	result.sigma = 0.7;
	GaussFitter g1;
	g1.setInitialParameters(result);
	GaussFitter g2(g1);
	TEST_REAL_EQUAL(g1.getInitialParameters().A, result.A)
	TEST_REAL_EQUAL(g1.getInitialParameters().x0, result.x0)
	TEST_REAL_EQUAL(g1.getInitialParameters().sigma, result.sigma)
	TEST_REAL_EQUAL(g1.getInitialParameters().A, g2.getInitialParameters().A)
	TEST_REAL_EQUAL(g1.getInitialParameters().x0, g2.getInitialParameters().x0)
	TEST_REAL_EQUAL(g1.getInitialParameters().sigma, g2.getInitialParameters().sigma)
}
RESULT

CHECK((virtual ~GaussFitter()))
{
  delete ptr;
	NOT_TESTABLE
}
RESULT

CHECK((GaussFitter& operator=(const GaussFitter & rhs)))
{
  GaussFitter::GaussFitResult result;
  result.A = 0.3;
  result.x0 = 0.7;
	result.sigma = 0.9;
  GaussFitter f1;
  f1.setInitialParameters(result);
  GaussFitter f2;
  f2 = f1;

	TEST_REAL_EQUAL(f1.getInitialParameters().A, result.A)
  TEST_REAL_EQUAL(f1.getInitialParameters().x0, result.x0)
	TEST_REAL_EQUAL(f1.getInitialParameters().sigma, result.sigma)
  
	TEST_REAL_EQUAL(f2.getInitialParameters().A, result.A)
  TEST_REAL_EQUAL(f2.getInitialParameters().x0, result.x0)
	TEST_REAL_EQUAL(f2.getInitialParameters().sigma, result.sigma)
	
  TEST_REAL_EQUAL(f1.getInitialParameters().A, f2.getInitialParameters().A)
  TEST_REAL_EQUAL(f1.getInitialParameters().x0, f2.getInitialParameters().x0)
	TEST_REAL_EQUAL(f1.getInitialParameters().sigma, f2.getInitialParameters().sigma)
}
RESULT

CHECK((GaussFitResult fit(std::vector< DPosition< 2 > >& points)))
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

	PRECISION(0.1)
	TEST_REAL_EQUAL(result.A, 1.0)
	TEST_REAL_EQUAL(result.x0, 0.3)
	TEST_REAL_EQUAL(result.sigma, 0.2)
}
RESULT

CHECK((const GaussFitResult& getInitialParameters() const ))
{
  NOT_TESTABLE // tested above
}
RESULT

CHECK((void setInitialParameters(const GaussFitResult& result)))
{
  GaussFitter f1;
  GaussFitter::GaussFitResult result = f1.getInitialParameters();
  TEST_REAL_EQUAL(result.A, 0.06)
  TEST_REAL_EQUAL(result.x0, 3.0)
	TEST_REAL_EQUAL(result.sigma, 0.5)
  result.A = 0.15;
  result.x0 = 0.24;
	result.sigma = 0.35;
  f1.setInitialParameters(result);

  TEST_REAL_EQUAL(f1.getInitialParameters().A, 0.15)
  TEST_REAL_EQUAL(f1.getInitialParameters().x0, 0.24)
	TEST_REAL_EQUAL(f1.getInitialParameters().sigma, 0.35)
}
RESULT

CHECK((const String& getGnuplotFormula() const ))
{
  String formula = ptr->getGnuplotFormula();
	// f(x)=1.01775 * exp(-(x - 0.300549) ** 2 / 2 / (0.136341) ** 2
	TEST_EQUAL(formula.hasSubstring("f(x)="), true)
	TEST_EQUAL(formula.hasSubstring(" * exp(-(x - 0.3"), true)
	TEST_EQUAL(formula.hasSubstring(") ** 2 / 2 / (0.1"), true)
	TEST_EQUAL(formula.hasSubstring(") ** 2"), true)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



