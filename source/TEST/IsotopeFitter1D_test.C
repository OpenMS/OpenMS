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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsotopeFitter1D* ptr = 0;
IsotopeFitter1D* nullPointer = 0;
START_SECTION(IsotopeFitter1D())
{
	ptr = new IsotopeFitter1D();
  TEST_EQUAL(ptr->getName(), "IsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((IsotopeFitter1D(const  IsotopeFitter1D &source)))
	IsotopeFitter1D isof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	isof1.setParameters(param);

	IsotopeFitter1D isof2(isof1);
  IsotopeFitter1D isof3;
	isof3.setParameters(param);
  isof1 = IsotopeFitter1D();
	TEST_EQUAL(isof3.getParameters(), isof2.getParameters())
END_SECTION

START_SECTION((virtual ~IsotopeFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual IsotopeFitter1D& operator=(const  IsotopeFitter1D &source)))
 	IsotopeFitter1D isof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	isof1.setParameters(param);

  IsotopeFitter1D isof2;
  isof2 = isof1;

  IsotopeFitter1D isof3;
	isof3.setParameters(param);

  isof1 = IsotopeFitter1D();
	TEST_EQUAL(isof3.getParameters(), isof3.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest 
	IsotopeFitter1D if1;
	if1 = IsotopeFitter1D();
	TEST_EQUAL(if1.getParameters(), if1.getParameters())
END_SECTION

START_SECTION((Fitter1D* create()))
  Fitter1D* ptr = IsotopeFitter1D::create();
  TEST_EQUAL(ptr->getName(), "IsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((const String getProductName()))
  TEST_EQUAL(IsotopeFitter1D::getProductName(),"IsotopeFitter1D")
  TEST_EQUAL(IsotopeFitter1D().getName(),"IsotopeFitter1D")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



