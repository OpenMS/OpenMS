// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter1D* ptr = nullptr;
GaussFitter1D* nullPointer = nullptr;
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
	delete ptr;
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



