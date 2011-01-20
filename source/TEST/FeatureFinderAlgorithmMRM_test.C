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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderAlgorithmMRM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderAlgorithmMRM<Peak1D, Feature>* ptr = 0;
START_SECTION(FeatureFinderAlgorithmMRM())
{
	ptr = new FeatureFinderAlgorithmMRM<Peak1D, Feature>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~FeatureFinderAlgorithmMRM())
{
	delete ptr;
}
END_SECTION

ptr = new FeatureFinderAlgorithmMRM<Peak1D, Feature>();

START_SECTION((virtual void run()))
{
	FeatureFinder ff;
  ff.setLogType(ProgressLogger::NONE);

	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmMRM_input.mzML"), exp);

	FeatureMap<> features, seeds;
	Param ff_param(ptr->getParameters());
	ff.run("mrm", exp, features, ff_param, seeds);

	TEST_EQUAL(exp.getChromatograms().size(), 3)

	FeatureMap<> new_features;
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

START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
{
  FeatureFinderAlgorithm<Peak1D, Feature>* ptr2 = 0;
	ptr2 = FeatureFinderAlgorithmMRM<Peak1D, Feature>::create();
	TEST_NOT_EQUAL(ptr2, 0)
	delete ptr;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_STRING_EQUAL(ptr->getProductName(), "mrm")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



