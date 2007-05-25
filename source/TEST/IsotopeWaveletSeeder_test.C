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
// $Maintainer: Ole Schulz-Trieglaff, Rene Hussong$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(IsotopeWaveletSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
IsotopeWaveletSeeder* ptr = 0;
CHECK((IsotopeWaveletSeeder()))
	ptr = new IsotopeWaveletSeeder();
  TEST_EQUAL(ptr->getName(), "IsotopeWaveletSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~IsotopeWaveletSeeder()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(IsotopeWaveletSeeder::getProductName(),"IsotopeWaveletSeeder")
	TEST_EQUAL(IsotopeWaveletSeeder().getName(),"IsotopeWaveletSeeder")
RESULT

CHECK((static BaseSeeder* create()))
	TEST_NOT_EQUAL(IsotopeWaveletSeeder::create(),0)
RESULT

CHECK((IsotopeWaveletSeeder& operator=(const IsotopeWaveletSeeder &rhs)))
	IsotopeWaveletSeeder ms1;
	IsotopeWaveletSeeder ms2;
	
	ms1 = ms2;
	
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK((IsotopeWaveletSeeder(const IsotopeWaveletSeeder &rhs)))
	IsotopeWaveletSeeder ms1;
	IsotopeWaveletSeeder ms2(ms1);
		
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK(([EXTRA]IndexSet nextSeed()))
	PRECISION(0.01)
	
  IsotopeWaveletSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment< > exp;
	MzDataFile().load("data/IsotopeWaveletTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);
	
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("min_number_scans",11);
	param.setValue("rt_tolerance_cluster",2.0);
	param.setValue("mass_tolerance_cluster",2.0);
	param.setValue("max_rt_dist_merging",0);
	param.setValue("max_mz_dist_merging",0);
	
	seeder.setParameters(param);
	
	FeaFiModule::IndexSet  region = seeder.nextSeed();
	
	ifstream infile( "data/IsotopeWaveletSeeder_region1");	
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
	
	region = seeder.nextSeed();	
	infile.open( "data/IsotopeWaveletSeeder_region2");	
	
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
	
	region = seeder.nextSeed();
	infile.open( "data/IsotopeWaveletSeeder_region3");	
	
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
	
	region = seeder.nextSeed();
	infile.open( "data/IsotopeWaveletSeeder_region4");	
	
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

	// test exception, there should be no more seeds
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )
 
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


