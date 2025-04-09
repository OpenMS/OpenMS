// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/BiGaussModel.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BiGaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BiGaussFitter1D* ptr = nullptr;
BiGaussFitter1D* nullPointer = nullptr;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


