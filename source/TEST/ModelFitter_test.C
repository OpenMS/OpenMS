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
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModelFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModelFitter* ptr = 0;
CHECK(ModelFitter())
{
	ptr = new ModelFitter();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~ModelFitter())
{
	delete ptr;
}
RESULT

CHECK((ModelFitter(const  MSExperiment< PeakType > *map, FeatureMap< FeatureType > *features, FeatureFinder *ff)))
{
  // TODO
}
RESULT

CHECK((virtual ~ModelFitter()))
{
  // TODO
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



