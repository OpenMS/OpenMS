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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter1D* ptr = 0;
GaussFitter1D* nullPointer = 0;
START_SECTION(GaussFitter1D())
{
	ptr = new GaussFitter1D();
  TEST_EQUAL(ptr->getName(), "GaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((GaussFitter1D(const  GaussFitter1D &source)))
	GaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
	gf1.setParameters(param);

	GaussFitter1D gf2(gf1);
  GaussFitter1D gf3;
	gf3.setParameters(param);
  gf1 = GaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((virtual ~GaussFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual GaussFitter1D& operator=(const  GaussFitter1D &source)))
 	GaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
	gf1.setParameters(param);

  GaussFitter1D gf2;
  gf2 = gf1;

  GaussFitter1D gf3;
	gf3.setParameters(param);

  gf1 = GaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
{
  Fitter1D* ptr = GaussFitter1D::create();
  TEST_EQUAL(ptr->getName(), "GaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((const String getProductName()))
{
  TEST_EQUAL(GaussFitter1D::getProductName(),"GaussFitter1D")
  TEST_EQUAL(GaussFitter1D().getName(),"GaussFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



