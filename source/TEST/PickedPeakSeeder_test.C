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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PickedPeakSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(PickedPeakSeeder, "$Id$")

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
PickedPeakSeeder* ptr = 0;
CHECK(PickedPeakSeeder())
	ptr = new PickedPeakSeeder();
  TEST_EQUAL(ptr->getName(), "PickedPeakSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PickedPeakSeeder())
	delete ptr;
RESULT

CHECK(nextSeed())
  PickedPeakSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment<DPeak<1> > exp;
	MzDataFile().load("data/PickedPeakTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);
	
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("min_number_scans",0);
	param.setValue("min_number_peaks",3);
	seeder.setParameters(param);
	
	FeaFiModule::IndexSet  region;
	FeaFiModule::IDX peak;
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),18084);	
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),2522);					
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),100108);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),9105);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),9118);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),9902);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),4412);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),9961);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),7708);			
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),44099);			

	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


