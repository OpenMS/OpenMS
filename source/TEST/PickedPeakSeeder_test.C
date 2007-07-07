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

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(PickedPeakSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
PickedPeakSeeder* ptr = 0;
CHECK((PickedPeakSeeder()))
	ptr = new PickedPeakSeeder();
  TEST_EQUAL(ptr->getName(), "PickedPeakSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~PickedPeakSeeder()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(PickedPeakSeeder::getProductName(),"PickedPeakSeeder")
	TEST_EQUAL(PickedPeakSeeder().getName(),"PickedPeakSeeder")
RESULT

CHECK((static BaseSeeder* create()))
	TEST_NOT_EQUAL(PickedPeakSeeder::create(),0)
RESULT

CHECK((PickedPeakSeeder& operator=(const PickedPeakSeeder &rhs)))
	PickedPeakSeeder ms1;
	PickedPeakSeeder ms2;
	
	ms1 = ms2;
	
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK((PickedPeakSeeder(const PickedPeakSeeder &rhs)))
	PickedPeakSeeder ms1;
	PickedPeakSeeder ms2(ms1);
		
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK(([EXTRA]IndexSet nextSeed()))
	PRECISION(0.01)	

  PickedPeakSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment<Peak1D > exp;
	MzDataFile().load("data/PickedPeakTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),1000);	
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("min_number_scans",5);
	param.setValue("max_rt_dist_merging",0);
	param.setValue("max_mz_dist_merging",0);
	seeder.setParameters(param);
	
	// test first seeding region	
	FeaFiModule::IndexSet region = seeder.nextSeed();
	ifstream infile( "data/PickedPeakSeeder_region1");	
	DoubleReal intensity, rt, mz;
	
	FeaFiModule::IndexSet::const_iterator citer = region.begin();
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		ABORT_IF(citer == region.end())
		
		TEST_REAL_EQUAL(traits->getPeakRt(*citer),rt)
		TEST_REAL_EQUAL(traits->getPeakMz(*citer),mz)
		TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),intensity)
				
		++citer;				
	}		
	infile.close();
	
	// retrieve second region
	region = seeder.nextSeed();
 	infile.open( "data/PickedPeakSeeder_region2");	
	
	citer = region.begin();
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		ABORT_IF(citer == region.end())
		
		TEST_REAL_EQUAL(traits->getPeakRt(*citer),rt)
		TEST_REAL_EQUAL(traits->getPeakMz(*citer),mz)
		TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),intensity)
				
		++citer;				
	}		
	infile.close();
		
	// done, should be the last region !
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() );

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

// 	ofstream outfile("region1");
// 	for(FeaFiModule::IndexSet::const_iterator citer = region.begin();
// 				citer != region.end();
// 				++citer)
// 	{
// 		outfile << traits->getPeakRt(*citer) << " " << traits->getPeakMz(*citer) << " " << traits->getPeakIntensity(*citer) << endl;
// 	}				
// 	outfile.close();

