// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>

///////////////////////////

START_TEST(SimpleSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

enum DimensionId
{
	RT = DimensionDescription < LCMS_Tag >::RT,
	MZ = DimensionDescription < LCMS_Tag >::MZ
};


// default ctor
SimpleSeeder* ptr = 0;
CHECK(Simple())
	ptr = new SimpleSeeder();
  TEST_EQUAL(ptr->getName(), "SimpleSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SimpleSeeder())
	delete ptr;
RESULT

MSExperiment<>::PeakType p;

// SPECTRUM 1
MSExperiment<>::SpectrumType s1;
s1.setRetentionTime(1.0);

p.setPos(500.0);
p.setIntensity(10.0);
s1.push_back(p);

p.setPos(600.0);
p.setIntensity(20.0);
s1.push_back(p);

p.setPos(800.0);
p.setIntensity(30.0);
s1.push_back(p);

p.setPos(1000.0);
p.setIntensity(40.0);
s1.push_back(p);

p.setPos(1200.0);
p.setIntensity(110.0);
s1.push_back(p);

// SPECTRUM 2
MSExperiment<>::SpectrumType s2;
s2.setRetentionTime(2.0);

p.setPos(500.0);
p.setIntensity(100.0);
s2.push_back(p);

p.setPos(600.0);
p.setIntensity(80.0);
s2.push_back(p);

p.setPos(800.0);
p.setIntensity(30.0);
s2.push_back(p);

p.setPos(1000.0);
p.setIntensity(10.0);
s2.push_back(p);

p.setPos(1200.0);
p.setIntensity(110.0);
s2.push_back(p);

CHECK(nextSeed())
  
	MSExperiment<> exp;
	exp.push_back(s1);
	exp.push_back(s2);
	
  FeaFiTraits* traits = new FeaFiTraits();
	traits->setData(exp.begin(), exp.end(),100);
	traits->getPeakFlag(make_pair(0,4)) = FeaFiTraits::INSIDE_FEATURE;
	traits->getPeakFlag(make_pair(1,4)) = FeaFiTraits::INSIDE_FEATURE;
	
	SimpleSeeder seeder;
	seeder.setTraits(traits);
	Param param;
	param.setValue("min_intensity",35);
	seeder.setParameters(param);
	FeaFiModule::IndexSet region;
	FeaFiModule::IDX peak;
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),100.0);
	TEST_EQUAL(traits->getPeakMz(peak),500.0);
	TEST_EQUAL(traits->getPeakRt(peak),2.0);
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),80.0);
	TEST_EQUAL(traits->getPeakMz(peak),600.0);
	TEST_EQUAL(traits->getPeakRt(peak),2.0);
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),40.0);
	TEST_EQUAL(traits->getPeakMz(peak),1000.0);
	TEST_EQUAL(traits->getPeakRt(peak),1.0);

	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


