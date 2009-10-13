// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinder, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinder* ptr = 0;
START_SECTION((FeatureFinder()))
	ptr = new FeatureFinder();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureFinder()))
	delete ptr;
END_SECTION

START_SECTION((template <class PeakType, class FeatureType> void run(const String &algorithm_name, MSExperiment< PeakType > const &input_map, FeatureMap< FeatureType > &features, const Param &param, const FeatureMap<FeatureType>& seeds)))
	FeatureFinder ff;
	FeatureMap<Feature> features;
	
	//empty map works -> nothing to do
	MSExperiment<Peak1D> map;
	ff.run("none", map, features, Param(), FeatureMap<Feature>());
	
	//no updateRanges -> exception
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap<Feature>()))
	
	//updateRanges -> it works again
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap<Feature>());
	
	//MS2 scans -> exception
	map[0].setMSLevel(1);
	map[0].setMSLevel(2);
	map.updateRanges();
	TEST_EXCEPTION(Exception::IllegalArgument, ff.run("none", map, features, Param(), FeatureMap<Feature>()))
END_SECTION

START_SECTION((const Flag& getPeakFlag(const IndexPair& index) const))
	FeatureFinder ff;
	FeatureMap<Feature> features;
	MSExperiment<Peak1D> map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap<Feature>());
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::UNUSED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Flag& getPeakFlag(const IndexPair& index)))
	FeatureFinder ff;
	FeatureMap<Feature> features;
	MSExperiment<Peak1D> map;
	map.resize(2);
	map[0].resize(1);
	map[1].resize(1);
	map.updateRanges();
	ff.run("none", map, features, Param(), FeatureMap<Feature>());
	ff.getPeakFlag(make_pair(0,0)) = FeatureFinderDefs::USED;
	TEST_EQUAL(ff.getPeakFlag(make_pair(0,0)),FeatureFinderDefs::USED)
	TEST_EQUAL(ff.getPeakFlag(make_pair(1,0)),FeatureFinderDefs::UNUSED)
END_SECTION

START_SECTION((Param getParameters(const String& algorithm_name) const))
	FeatureFinder ff;
	TEST_EQUAL(ff.getParameters("none")==Param(),true)
	TEST_EQUAL(ff.getParameters("simple")==Param(),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



