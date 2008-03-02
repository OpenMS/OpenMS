// -*- Mode: C++; tab-width: 2; -*-
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

CHECK(~GammaDistributionFitter())
{
	delete ptr;
}
RESULT

CHECK((GammaDistributionFitter(const GammaDistributionFitter &)))
{
  // TODO
}
RESULT

CHECK((virtual ~GammaDistributionFitter()))
{
  // TODO
}
RESULT

CHECK((GammaDistributionFitter& operator=(const GammaDistributionFitter &)))
{
  // TODO
}
RESULT

CHECK((GammaDistributionFitResult fit(std::vector< DPosition< 2 > > &)))
{
  // TODO
}
RESULT

CHECK((const GammaDistributionFitResult& getInitialParameters() const ))
{
  // TODO
}
RESULT

CHECK((void setInitialParameters(const GammaDistributionFitResult &)))
{
  // TODO
}
RESULT

CHECK((const String& getGnuplotFormula() const ))
{
  // TODO
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



