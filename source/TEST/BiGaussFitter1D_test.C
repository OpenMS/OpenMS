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

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BiGaussFitter1D, "$Id: BiGaussFitter1D_test.C 2642 2008-01-09 15:22:14Z grunert $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BiGaussFitter1D* ptr = 0;
CHECK(BiGaussFitter1D())
{
  ptr = new BiGaussFitter1D();
  TEST_EQUAL(ptr->getName(), "BiGaussFitter1D")
      TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~BiGaussFitter1D())
{
  delete ptr;
}
RESULT

CHECK((BiGaussFitter1D(const  BiGaussFitter1D &source)))
{
  // TODO
}
RESULT

CHECK((virtual ~BiGaussFitter1D()))
{
  // TODO
}
RESULT

CHECK((virtual BiGaussFitter1D& operator=(const  BiGaussFitter1D &source)))
{
  // TODO
}
RESULT

CHECK((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
{
  // TODO
}
RESULT

CHECK((Fitter1D* create()))
{
  Fitter1D* ptr = BiGaussFitter1D::create();
  TEST_EQUAL(ptr->getName(), "BiGaussFitter1D")
  TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK((const String getProductName()))
{
  TEST_EQUAL(BiGaussFitter1D::getProductName(),"BiGaussFitter1D")
  TEST_EQUAL(BiGaussFitter1D().getName(),"BiGaussFitter1D")
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


