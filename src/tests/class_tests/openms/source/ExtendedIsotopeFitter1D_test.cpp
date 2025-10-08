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

#include <OpenMS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExtendedIsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExtendedIsotopeFitter1D* ptr = nullptr;
ExtendedIsotopeFitter1D* nullPointer = nullptr;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



