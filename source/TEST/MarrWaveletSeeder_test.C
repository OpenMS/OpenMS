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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MarrWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(MarrWaveletSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
MarrWaveletSeeder* ptr = 0;
CHECK((MarrWaveletSeeder()))
	ptr = new MarrWaveletSeeder();
  TEST_EQUAL(ptr->getName(), "MarrWaveletSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~MarrWaveletSeeder()))
	delete ptr;
RESULT

CHECK((IndexSet nextSeed()))
	
	PRECISION(0.01)
	
  MarrWaveletSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment<Peak1D > exp;
	MzDataFile().load("data/MarrWaveletTestData.mzData",exp);
	traits->setData(exp.begin(), exp.end(),100);
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("min_number_scans",4);
	param.setValue("noise_level_signal",10000);
	param.setValue("noise_level_cwt",10000);
	param.setValue("scans_to_sumup",4);
	param.setValue("cwt_scale",0.1);
	seeder.setParameters(param);
	
	// test first seeding region	
	FeaFiModule::IndexSet region = seeder.nextSeed();
	
	ifstream infile( "data/MarrWaveletSeeder_region1");	
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
	
	// retrieve second region
	region = seeder.nextSeed();
	
	infile.open( "data/MarrWaveletSeeder_region2");	
	
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
	
	// done, should be the last region !
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() );
	 
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


