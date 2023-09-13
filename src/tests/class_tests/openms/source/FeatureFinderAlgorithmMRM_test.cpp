// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderAlgorithmMRM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderAlgorithmMRM* ptr = nullptr;
FeatureFinderAlgorithmMRM* nullPointer = nullptr;
FeatureFinderAlgorithm* ffA_nullPointer = nullptr;

START_SECTION(FeatureFinderAlgorithmMRM())
{
  ptr = new FeatureFinderAlgorithmMRM();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FeatureFinderAlgorithmMRM())
{
	delete ptr;
}
END_SECTION

ptr = new FeatureFinderAlgorithmMRM();

START_SECTION((virtual void run()))
{
	FeatureFinder ff;
  ff.setLogType(ProgressLogger::NONE);

  PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmMRM_input.mzML"), exp);

	FeatureMap features, seeds;
	Param ff_param(ptr->getParameters());
	ff.run("mrm", exp, features, ff_param, seeds);

	TEST_EQUAL(exp.getChromatograms().size(), 3)

	FeatureMap new_features;
	for (Size i = 0; i != features.size(); ++i)
	{
		if (features[i].getQuality(0) > 0.99)
		{
			new_features.push_back(features[i]);
		}
	}

	TEST_EQUAL(new_features.size(), 3)

	for (Size i = 0; i != new_features.size(); ++i)
	{
		TEST_EQUAL(new_features[i].getIntensity() > 100000, true)
	}
}
END_SECTION

START_SECTION((static FeatureFinderAlgorithm<PeakType>* create()))
{
  FeatureFinderAlgorithm* ptr2 = nullptr;
  ptr2 = FeatureFinderAlgorithmMRM::create();
  TEST_NOT_EQUAL(ptr2, ffA_nullPointer)
  delete ptr2;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_STRING_EQUAL(ptr->getProductName(), "mrm")
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



