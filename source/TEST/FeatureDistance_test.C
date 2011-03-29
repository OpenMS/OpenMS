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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

///////////////////////////

START_TEST(FeatureDistance, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeatureDistance* d_ptr = 0;
FeatureDistance* d_nullPointer = 0;
START_SECTION((FeatureDistance(DoubleReal max_intensity=1.0, bool force_constraints=false)))
{
	d_ptr = new FeatureDistance();
  TEST_NOT_EQUAL(d_ptr, d_nullPointer);
}
END_SECTION

START_SECTION((~FeatureDistance()))
{
	delete d_ptr;
}
END_SECTION

START_SECTION((std::pair<bool, DoubleReal> operator()(const BaseFeature& left, const BaseFeature& right)))
{
	FeatureDistance dist(1000.0, false);
	Param param = dist.getDefaults();
	param.setValue("distance_RT:max_difference", 100.0);
	param.setValue("distance_MZ:max_difference", 1.0);
	param.setValue("distance_MZ:exponent", 1.0);
	param.setValue("distance_intensity:weight", 1.0);
	dist.setParameters(param);
	BaseFeature left, right;
	left.setRT(100.0);
	left.setMZ(100.0);
	left.setIntensity(100.0);
	// all distance components vary by 10% of the maximum:
	right.setRT(110.0);
	right.setMZ(100.1);
	right.setIntensity(200.0);
	pair<bool, DoubleReal> result = dist(left, right);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 0.1);
	// no differences:
	result = dist(left, left);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 0.0);
	// differences at maximum:
	right.setRT(200.0);
	right.setMZ(101.0);
	right.setIntensity(1000.0);
	left.setIntensity(0.0);
	result = dist(left, right);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 1.0);
	// differences beyond maximum:
	right.setRT(300.0);
	result = dist(left, right);
	TEST_EQUAL(result.first, false);
	TEST_REAL_SIMILAR(result.second, 1.33333333);
	FeatureDistance dist2(1000.0, true);
	result = dist2(left, right);
	TEST_EQUAL(result.first, false);
	TEST_EQUAL(result.second, FeatureDistance::infinity);
}
END_SECTION

START_SECTION((FeatureDistance& operator=(const FeatureDistance& other)))
{
	FeatureDistance dist(1000.0, true);
	Param param = dist.getDefaults();
	param.setValue("distance_RT:max_difference", 100.0);
	param.setValue("distance_MZ:max_difference", 1.0);
	param.setValue("distance_MZ:exponent", 1.0);
	param.setValue("distance_intensity:weight", 1.0);
	dist.setParameters(param);
	FeatureDistance dist2;
	dist2 = dist;
	TEST_EQUAL(dist.getParameters(), dist2.getParameters());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
