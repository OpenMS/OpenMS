// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LmaIsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


LmaIsotopeFitter1D* ptr = 0;
START_SECTION(LmaIsotopeFitter1D())
{
  ptr = new LmaIsotopeFitter1D();
  TEST_EQUAL(ptr->getName(), "LmaIsotopeFitter1D")
  TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((LmaIsotopeFitter1D(const  LmaIsotopeFitter1D &source)))
	LmaIsotopeFitter1D lisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 2 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	lisof1.setParameters(param);

	LmaIsotopeFitter1D lisof2(lisof1);
  LmaIsotopeFitter1D lisof3;
	lisof3.setParameters(param);
  lisof1 = LmaIsotopeFitter1D();
	TEST_EQUAL(lisof3.getParameters(), lisof2.getParameters())
END_SECTION

START_SECTION((virtual ~LmaIsotopeFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual LmaIsotopeFitter1D& operator=(const  LmaIsotopeFitter1D &source)))
	LmaIsotopeFitter1D lisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 2 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	lisof1.setParameters(param);
	
  LmaIsotopeFitter1D lisof2;
  lisof2 = lisof1;

  LmaIsotopeFitter1D lisof3;
	lisof3.setParameters(param);

  lisof1 = LmaIsotopeFitter1D();
	TEST_EQUAL(lisof3.getParameters(), lisof3.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
{
  Fitter1D* ptr = LmaIsotopeFitter1D::create();
  TEST_EQUAL(ptr->getName(), "LmaIsotopeFitter1D")
  TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((const String getProductName()))
{
  TEST_EQUAL(LmaIsotopeFitter1D::getProductName(),"LmaIsotopeFitter1D")
  TEST_EQUAL(LmaIsotopeFitter1D().getName(),"LmaIsotopeFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



