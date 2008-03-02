// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $ Marcel Grunert
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModelFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef RawDataPoint1D PeakType;
typedef Feature FeatureType;
typedef ModelFitter<PeakType,FeatureType> ModelFitterType;
typedef FeatureFinderDefs::ChargedIndexSet ChargedIndexSet;

ModelFitterType* ptr = 0;
CHECK((ModelFitter(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)))
	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;
	ptr = new ModelFitterType(&input,&features,&ff);
  TEST_EQUAL(ptr->getName(), "ModelFitter")
	TEST_NOT_EQUAL(ptr, 0)
RESULT


CHECK(~ModelFitter())
{
	delete ptr;
}
RESULT

CHECK((virtual ~ModelFitter()))
{
}
RESULT

CHECK((static const String getName()))
{
	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;
	ModelFitterType model_fitter(&input,&features,&ff);
	TEST_EQUAL(model_fitter.getName(),"ModelFitter");
}
RESULT

CHECK((Feature fit(const  ChargedIndexSet &index_set)))
{
  // TODO
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



