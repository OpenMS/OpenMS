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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExtendedIsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExtendedIsotopeFitter1D* ptr = 0;
ExtendedIsotopeFitter1D* nullPointer = 0;
START_SECTION(ExtendedIsotopeFitter1D())
{
	ptr = new ExtendedIsotopeFitter1D();
        TEST_EQUAL(ptr->getName(), "ExtendedIsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((ExtendedIsotopeFitter1D(const  ExtendedIsotopeFitter1D &source)))
	ExtendedIsotopeFitter1D eisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	eisof1.setParameters(param);

	ExtendedIsotopeFitter1D eisof2(eisof1);
  ExtendedIsotopeFitter1D eisof3;
	eisof3.setParameters(param);
  eisof1 = ExtendedIsotopeFitter1D();
	TEST_EQUAL(eisof3.getParameters(), eisof2.getParameters())

END_SECTION

START_SECTION((virtual ~ExtendedIsotopeFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual ExtendedIsotopeFitter1D& operator=(const  ExtendedIsotopeFitter1D &source)))
	ExtendedIsotopeFitter1D eisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	eisof1.setParameters(param);

  ExtendedIsotopeFitter1D eisof2;
  eisof2 = eisof1;

  ExtendedIsotopeFitter1D eisof3;
	eisof3.setParameters(param);

  eisof1 = ExtendedIsotopeFitter1D();
	TEST_EQUAL(eisof3.getParameters(), eisof3.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest TODO
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
  Fitter1D* ptr = ExtendedIsotopeFitter1D::create();
  TEST_EQUAL(ptr->getName(), "ExtendedIsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((const String getProductName()))
  TEST_EQUAL(ExtendedIsotopeFitter1D::getProductName(),"ExtendedIsotopeFitter1D")
  TEST_EQUAL(ExtendedIsotopeFitter1D().getName(),"ExtendedIsotopeFitter1D")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



