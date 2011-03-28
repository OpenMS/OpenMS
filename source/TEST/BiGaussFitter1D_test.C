// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BiGaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BiGaussFitter1D* ptr = 0;
BiGaussFitter1D* nullPointer = 0;
START_SECTION(BiGaussFitter1D())
{
  ptr = new BiGaussFitter1D();
  TEST_EQUAL(ptr->getName(), "BiGaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((BiGaussFitter1D(const  BiGaussFitter1D &source)))
{
	BiGaussFitter1D bgf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "statistics:variance1", 2.0 );
  param.setValue( "statistics:variance2", 5.0 );
  param.setValue( "interpolation_step", 1.0 );
	bgf1.setParameters(param);

	BiGaussFitter1D bgf2(bgf1);
  BiGaussFitter1D bgf3;
	bgf3.setParameters(param);
  bgf1 = BiGaussFitter1D();
	TEST_EQUAL(bgf3.getParameters(), bgf2.getParameters())
}
END_SECTION

START_SECTION((virtual ~BiGaussFitter1D()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual BiGaussFitter1D& operator=(const  BiGaussFitter1D &source)))
{
  BiGaussFitter1D bgf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "statistics:variance1", 2.0 );
  param.setValue( "statistics:variance2", 5.0 );
  param.setValue( "interpolation_step", 1.0 );
	bgf1.setParameters(param);

  BiGaussFitter1D bgf2;
  bgf2 = bgf1;

  BiGaussFitter1D bgf3;
	bgf3.setParameters(param);

  bgf1 = BiGaussFitter1D();
	TEST_EQUAL(bgf3.getParameters(), bgf2.getParameters())
}
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
{
  Fitter1D* ptr = BiGaussFitter1D::create();
  TEST_EQUAL(ptr->getName(), "BiGaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((const String getProductName()))
{
  TEST_EQUAL(BiGaussFitter1D::getProductName(),"BiGaussFitter1D")
  TEST_EQUAL(BiGaussFitter1D().getName(),"BiGaussFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


