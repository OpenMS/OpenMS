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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummySeeder.h>

#include <OpenMS/FORMAT/MzDataFile.h>

#include<iostream>  
#include<fstream> 

using namespace OpenMS;

START_TEST(DummySeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
	
using namespace OpenMS;
using namespace std;

// default ctor
DummySeeder* ptr = 0;
CHECK((DummySeeder()))
	ptr = new DummySeeder();
  TEST_EQUAL(ptr->getName(), "DummySeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~DummySeeder()))
	delete ptr;
RESULT

CHECK((IndexSet nextSeed()))
	PRECISION(0.01)
	
  DummySeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment< > exp;
	MzDataFile().load("data/DummySeederTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);	
	seeder.setTraits(traits);
	
	Param param;
	param.setValue("min_snratio",2.0);
	param.setValue("min_number_scans",7.0);
	seeder.setParameters(param);
	
	FeaFiModule::IndexSet  region;
	
	// there should be two (seeding) regions for this data set, check the first one
	region = seeder.nextSeed();
		
	ifstream infile( "data/DummySeederTestData_region1");
	DoubleReal intensity, rt, mz;

	FeaFiModule::IndexSet::const_iterator citer = region.begin();
	
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		
		TEST_REAL_EQUAL(traits->getPeakRt(*citer),rt)
		TEST_REAL_EQUAL(traits->getPeakMz(*citer),mz)
		TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),intensity)
				
		++citer;		
	}	
	infile.close();
	
	// check second region
	region = seeder.nextSeed();		
	infile.open( "data/DummySeederTestData_region2");
	citer = region.begin();
	
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		
		TEST_REAL_EQUAL(traits->getPeakRt(*citer),rt)
		TEST_REAL_EQUAL(traits->getPeakMz(*citer),mz)
		TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),intensity)
				
		++citer;				
	}	
	infile.close();
	
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )	
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(DummySeeder::getProductName(),"DummySeeder");
RESULT

CHECK((static BaseSeeder* create()))
	BaseSeeder* base = DummySeeder::create();
	TEST_NOT_EQUAL(base,0);
	delete(base);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
