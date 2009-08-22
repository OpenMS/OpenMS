// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GaussFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter* ptr = 0;
START_SECTION(GaussFitter())
{
	ptr = new GaussFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~GaussFitter()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((GaussFitResult fit(std::vector< DPosition< 2 > >& points)))
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

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result.A, 1.0)
	TEST_REAL_SIMILAR(result.x0, 0.3)
	TEST_REAL_SIMILAR(result.sigma, 0.2)
}
END_SECTION

START_SECTION((void setInitialParameters(const GaussFitResult& result)))
{
  GaussFitter f1;
  GaussFitter::GaussFitResult result;
  f1.setInitialParameters(result);
	
	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((const String& getGnuplotFormula() const ))
{
  String formula = ptr->getGnuplotFormula();
	// f(x)=1.01775 * exp(-(x - 0.300549) ** 2 / 2 / (0.136341) ** 2
	TEST_EQUAL(formula.hasSubstring("f(x)="), true)
	TEST_EQUAL(formula.hasSubstring(" * exp(-(x - 0.3"), true)
	TEST_EQUAL(formula.hasSubstring(") ** 2 / 2 / (0.1"), true)
	TEST_EQUAL(formula.hasSubstring(") ** 2"), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



