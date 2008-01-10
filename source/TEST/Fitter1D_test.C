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
// $Maintainer: $ Marcel Grunert
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

using namespace OpenMS;
using namespace std;

START_TEST(Fitter1D, "$Id: Fitter1D_test.C 2642 2008-01-09 11:53:14Z grunert $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Fitter1D* ptr = 0;
CHECK(Fitter1D())
{
	ptr = new Fitter1D();
        TEST_EQUAL(ptr->getName(), "Fitter1D")
        TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~Fitter1D())
{
	delete ptr;
}
RESULT

CHECK((Fitter1D(const  Fitter1D &source)))
{
  // TODO
}
RESULT

CHECK((virtual ~Fitter1D()))
{
  // TODO
}
RESULT

CHECK((virtual Fitter1D& operator=(const  Fitter1D &source)))
{
  // TODO
}
RESULT

CHECK((virtual QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)=0))
{
  // TODO
}
RESULT

CHECK((void registerChildren()))
{
  // TODO
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



