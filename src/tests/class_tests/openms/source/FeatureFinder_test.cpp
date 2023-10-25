// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
///////////////////////////


#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinder* ptr = nullptr;
FeatureFinder* nullPointer = nullptr;
START_SECTION((FeatureFinder()))
	ptr = new FeatureFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureFinder()))
	delete ptr;
END_SECTION

START_SECTION((template <class PeakType, class FeatureType> void run(const String &algorithm_name, MSExperiment< PeakType > const &input_map, FeatureMap< FeatureType > &features, const Param &param, const FeatureMap<FeatureType>& seeds)))
	FeatureFinder ff;
	FeatureMap features;

	//empty map works -> nothing to do
	PeakMap map;
	ff.run("none", map, features, Param(), FeatureMap());

	//no updateRanges -> exception
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap()))

	//updateRanges -> it works again
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());

	//MS2 scans -> exception
	map[0].setMSLevel(1);
	map[0].setMSLevel(2);
	map.updateRanges();
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap()))
END_SECTION

START_SECTION((const Flag& getPeakFlag(const IndexPair& index) const))
	FeatureFinder ff;
	FeatureMap features;
	PeakMap map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::UNUSED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Flag& getPeakFlag(const IndexPair& index)))
	FeatureFinder ff;
	FeatureMap features;
	PeakMap map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap());
	ff.getPeakFlag(make_pair(0,0)) = FeatureFinderDefs::USED;
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::USED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Param getParameters(const String& algorithm_name) const))
	FeatureFinder ff;
	TEST_EQUAL(ff.getParameters("none").empty(),true)
	TEST_EQUAL(ff.getParameters("centroided").empty(),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



