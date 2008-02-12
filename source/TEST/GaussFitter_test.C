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

CHECK(~GaussFitter())
{
	delete ptr;
}
RESULT

CHECK((GaussFitter(const GaussFitter &)))
{
  // TODO
}
RESULT

CHECK((virtual ~GaussFitter()))
{
  // TODO
}
RESULT

CHECK((GaussFitter& operator=(const GaussFitter &)))
{
  // TODO
}
RESULT

CHECK((GaussFitResult fit(std::vector< DPosition< 2 > > &)))
{
  // TODO
}
RESULT

CHECK((const GaussFitResult& getInitialParameters() const ))
{
  // TODO
}
RESULT

CHECK((void setInitialParameters(const GaussFitResult &)))
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



