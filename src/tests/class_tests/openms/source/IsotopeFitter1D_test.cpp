// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/IsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsotopeFitter1D* ptr = nullptr;
IsotopeFitter1D* nullPointer = nullptr;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



