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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GumbelDistributionFitter, "$Id: GaussFitter_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GumbelDistributionFitter* ptr = 0;
START_SECTION(GumbelDistributionFitter())
{
	ptr = new GumbelDistributionFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~GumbelDistributionFitter()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((GumbelDistributionFitResult GumbelDistributionFitter::fit(std::vector< DPosition< 2 > >& points)))
{

		DPosition<2> pos;
  vector<DPosition<2> > points;

	pos.setX(-2.7); pos.setY(0.017); points.push_back(pos);
	pos.setX(-2.5); pos.setY(0.025); points.push_back(pos);
	pos.setX(-2); pos.setY(0.052); points.push_back(pos);
	pos.setX(-1); pos.setY(0.127); points.push_back(pos);
	pos.setX(-0.7); pos.setY(0.147); points.push_back(pos);
	pos.setX(-0.01); pos.setY(0.178); points.push_back(pos);
	pos.setX(0); pos.setY(0.178); points.push_back(pos);
	pos.setX(0.2); pos.setY(0.182); points.push_back(pos);
	pos.setX(0.5); pos.setY(0.184); points.push_back(pos);
	pos.setX(1); pos.setY(0.179); points.push_back(pos);
	pos.setX(1.3); pos.setY(0.171); points.push_back(pos);
	pos.setX(1.9); pos.setY(0.151); points.push_back(pos);
	pos.setX(2.5); pos.setY(0.127); points.push_back(pos);
	pos.setX(2.6); pos.setY(0.123); points.push_back(pos);
	pos.setX(2.7); pos.setY(0.119); points.push_back(pos);
	pos.setX(2.8); pos.setY(0.115); points.push_back(pos);
	pos.setX(2.9); pos.setY(0.111); points.push_back(pos);
	pos.setX(3); pos.setY(0.108); points.push_back(pos);
	pos.setX(3.5); pos.setY(0.089); points.push_back(pos);
	pos.setX(3.9); pos.setY(0.076); points.push_back(pos);
	pos.setX(4.01); pos.setY(0.073); points.push_back(pos);
	pos.setX(4.22); pos.setY(0.067); points.push_back(pos);
	pos.setX(4.7); pos.setY(0.054); points.push_back(pos);
	pos.setX(4.9); pos.setY(0.05); points.push_back(pos);
	pos.setX(5); pos.setY(0.047); points.push_back(pos);
	pos.setX(6); pos.setY(0.03); points.push_back(pos);
	pos.setX(7); pos.setY(0.017); points.push_back(pos);
	pos.setX(7.5); pos.setY(0.015); points.push_back(pos);
	pos.setX(7.9); pos.setY(0.012); points.push_back(pos);
	pos.setX(8.03); pos.setY(0.011); points.push_back(pos);
	//a= 0.5, b = 2
	

	ptr = new GumbelDistributionFitter;
		GumbelDistributionFitter::GumbelDistributionFitResult init_param;
	init_param.a = 1.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result = ptr->fit(points);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result.a, 0.5)
	TEST_REAL_SIMILAR(result.b, 2.0)	
	
	vector<DPosition<2> > points2;
	pos.setX(0); pos.setY(0.18); points2.push_back(pos);
	pos.setX(0.2); pos.setY(0.24); points2.push_back(pos);
	pos.setX(0.5); pos.setY(0.32); points2.push_back(pos);
	pos.setX(1); pos.setY(0.37); points2.push_back(pos);
	pos.setX(1.3); pos.setY(0.35); points2.push_back(pos);
	pos.setX(1.9); pos.setY(0.27); points2.push_back(pos);
	pos.setX(2.5); pos.setY(0.18); points2.push_back(pos);
	pos.setX(2.6); pos.setY(0.16); points2.push_back(pos);
	pos.setX(3); pos.setY(0.12); points2.push_back(pos);
	pos.setX(5); pos.setY(0.02); points2.push_back(pos);
	//a = 1, b = 1
	
	init_param.a = 3.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result2 = ptr->fit(points2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result2.a, 1.0)
	TEST_REAL_SIMILAR(result2.b, 1.0)	
	
}
END_SECTION

START_SECTION((void GumbelDistributionFitter::setInitialParameters(const GumbelDistributionFitResult& result)))
{
  GumbelDistributionFitter f1;
  GumbelDistributionFitter::GumbelDistributionFitResult result;
  f1.setInitialParameters(result);
	
	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((const String& GumbelDistributionFitter::getGnuplotFormula() const ))
{
  String formula = ptr->getGnuplotFormula();
	// f(x)=(1/1) * exp(-(x - 1)/1) * exp(-1 * exp(-(x-1)/1))
	TEST_EQUAL(formula.hasSubstring("f(x)="), true)
	TEST_EQUAL(formula.hasSubstring("* exp(-(x -"), true)
	TEST_EQUAL(formula.hasSubstring("* exp(-1*exp(-(x-"), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



