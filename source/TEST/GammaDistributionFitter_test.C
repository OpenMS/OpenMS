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
using namespace Math;
using namespace std;

START_TEST(GammaDistributionFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GammaDistributionFitter* ptr = 0;
START_SECTION(GammaDistributionFitter())
{
	ptr = new GammaDistributionFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~GammaDistributionFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((GammaDistributionFitResult fit(std::vector< DPosition< 2 > > & points)))
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

  TOLERANCE_ABSOLUTE(0.1)
  TEST_REAL_SIMILAR(result.b, 1.0)
  TEST_REAL_SIMILAR(result.p, 0.3)
	*/
}
END_SECTION

START_SECTION((void setInitialParameters(const GammaDistributionFitResult & result)))
{
  GammaDistributionFitter f1;
  GammaDistributionFitter::GammaDistributionFitResult result;
  f1.setInitialParameters(result);
	
	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((const String& getGnuplotFormula() const ))
	// TODO
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



