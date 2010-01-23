// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
///////////////////////////

START_TEST(SimpleSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
SimpleSeeder<Peak1D,Feature>* ptr = 0;
START_SECTION(SimpleSeeder(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff))
	MSExperiment<Peak1D> exp;
	ptr = new SimpleSeeder<Peak1D,Feature>(&exp,0,0);
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(virtual ~SimpleSeeder())
	delete ptr;
END_SECTION


START_SECTION(IndexPair nextSeed())

	//create map
	MSExperiment<Peak1D> exp;
	MSExperiment<Peak1D>::PeakType p;
	
	exp.resize(1);
	exp.back().setRT(1.0);
	p.setMZ(500.0);
	p.setIntensity(10.0);
	exp.back().push_back(p);
	p.setMZ(600.0);
	p.setIntensity(20.0);
	exp.back().push_back(p);
	p.setMZ(800.0);
	p.setIntensity(30.0);
	exp.back().push_back(p);
	p.setMZ(1000.0);
	p.setIntensity(40.0);
	exp.back().push_back(p);
	p.setMZ(1200.0);
	p.setIntensity(110.0);
	exp.back().push_back(p);
	
	exp.resize(2);
	exp.back().setRT(2.0);
	p.setMZ(500.0);
	p.setIntensity(101.0);
	exp.back().push_back(p);
	p.setMZ(600.0);
	p.setIntensity(81.0);
	exp.back().push_back(p);
	p.setMZ(800.0);
	p.setIntensity(31.0);
	exp.back().push_back(p);
	p.setMZ(1000.0);
	p.setIntensity(11.0);
	exp.back().push_back(p);
	p.setMZ(1200.0);
	p.setIntensity(111.0);
	exp.back().push_back(p);
	
	exp.updateRanges();
	
	//create dummy FF instance (for flags)
  FeatureFinder ff;
  FeatureMap<Feature> features;
  ff.run("none", exp,features, Param(), FeatureMap<>());
  
  //make two peaks used
	ff.getPeakFlag(make_pair(0,4)) = FeatureFinderDefs::USED;
	ff.getPeakFlag(make_pair(1,4)) = FeatureFinderDefs::USED;
	
	//First test (3 unused peaks with intensity greater than 35)
	{
		SimpleSeeder<Peak1D,Feature> seeder(&exp, &features, &ff);
		Param param;
		param.setValue("min_intensity",35.0);
		param.setValue("signal_to_noise",0.0);
		
		seeder.setParameters(param);
		
		FeatureFinderDefs::IndexPair peak;
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,0);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,1);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,0);
		TEST_EQUAL(peak.second,3);
	
		TEST_EXCEPTION( FeatureFinderDefs::NoSuccessor , seeder.nextSeed() )
	}
	
	
	// Second test (2 unused peaks with intensity greater than 27,75)
	{
		SimpleSeeder<Peak1D,Feature> seeder(&exp, &features, &ff);
		Param param;
		param.setValue("min_intensity",0.0);
		param.setValue("signal_to_noise",0.1);
    param.setValue("SignalToNoiseEstimationParameter:win_len",1.0);
    param.setValue("SignalToNoiseEstimationParameter:bin_count",4);
    param.setValue("SignalToNoiseEstimationParameter:min_required_elements",1);
		
		seeder.setParameters(param);
		
		FeatureFinderDefs::IndexPair peak;
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,1);
		TEST_EQUAL(peak.second,2);
		
		peak = seeder.nextSeed();
		TEST_EQUAL(peak.first,0);
		TEST_EQUAL(peak.second,2);
			
	//	TEST_EXCEPTION( FeatureFinderDefs::NoSuccessor , seeder.nextSeed() )		
	}
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


