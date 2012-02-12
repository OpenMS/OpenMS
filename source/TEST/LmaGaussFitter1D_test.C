// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LmaGaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LmaGaussFitter1D* ptr = 0;
LmaGaussFitter1D* nullPointer = 0;
START_SECTION(LmaGaussFitter1D())
	ptr = new LmaGaussFitter1D();
  TEST_EQUAL(ptr->getName(), "LmaGaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((LmaGaussFitter1D(const  LmaGaussFitter1D &source)))
	LmaGaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	gf1.setParameters(param);

	LmaGaussFitter1D gf2(gf1);
  LmaGaussFitter1D gf3;
	gf3.setParameters(param);
  gf1 = LmaGaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((virtual ~LmaGaussFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual LmaGaussFitter1D& operator=(const  LmaGaussFitter1D &source)))
	LmaGaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	gf1.setParameters(param);

  LmaGaussFitter1D gf2;
  gf2 = gf1;

  LmaGaussFitter1D gf3;
	gf3.setParameters(param);

  gf1 = LmaGaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *& model)))
	// dummy subtest 
	LmaGaussFitter1D gf1;
	gf1 = LmaGaussFitter1D();
	TEST_EQUAL(gf1.getParameters(), gf1.getParameters())
END_SECTION

START_SECTION((Fitter1D* create()))
  Fitter1D* ptr = LmaGaussFitter1D::create();
  TEST_EQUAL(ptr->getName(), "LmaGaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((const String getProductName()))
  TEST_EQUAL(LmaGaussFitter1D::getProductName(),"LmaGaussFitter1D")
  TEST_EQUAL(LmaGaussFitter1D().getName(),"LmaGaussFitter1D")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



