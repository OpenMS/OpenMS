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
CHECK(PickedPeakSeeder())
	ptr = new PickedPeakSeeder();
  TEST_EQUAL(ptr->getName(), "PickedPeakSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PickedPeakSeeder())
	delete ptr;
RESULT

CHECK(nextSeed())
	PRECISION(0.01)	

  PickedPeakSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment<Peak1D > exp;
	MzDataFile().load("data/PickedPeakTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);	
	seeder.setTraits(traits);
	
	// test first region
	FeaFiModule::IndexSet  region = seeder.nextSeed();	
	TEST_EQUAL(region.size(),29)
	FeaFiModule::IndexSet::const_iterator citer = region.begin();

	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),619);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),695.09);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1938.13);	
	
	++citer;
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),452);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),695.44);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1938.13);	
	
	double max_intensity = 0.0;
	FeaFiModule::IDX max_peak;
	for (citer = region.begin(); citer != region.end();++citer)
	{
		if (traits->getPeakIntensity(*citer) > max_intensity)
		{
			max_intensity = traits->getPeakIntensity(*citer);
			max_peak = *citer;
		}  	
	}
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(max_peak),2139);	
	TEST_REAL_EQUAL(traits->getPeakMz(max_peak),695.082);	
	TEST_REAL_EQUAL(traits->getPeakRt(max_peak),1940.34);	
	
	citer = region.end();
	--citer;--citer;--citer;
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),336);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),694.753);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1946.98);	
	
	++citer;
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),296);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),695.065);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1946.98);	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// test second region
	region = seeder.nextSeed();	
	TEST_EQUAL(region.size(),22)
	citer = region.begin();

	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),832);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),713.811);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1932.6);	
	
	++citer;
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),576);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),714.136);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1932.6);	
	
	max_intensity = 0.0;
	for (citer = region.begin(); citer != region.end();++citer)
	{
		if (traits->getPeakIntensity(*citer) > max_intensity)
		{
			max_intensity = traits->getPeakIntensity(*citer);
			max_peak = *citer;
		}  	
	}
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(max_peak),1683);	
	TEST_REAL_EQUAL(traits->getPeakMz(max_peak),713.804);	
	TEST_REAL_EQUAL(traits->getPeakRt(max_peak),1935.92);	
	
	citer = region.end();
	--citer;--citer;--citer;
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),1044);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),713.807);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1938.13);	
	
	++citer;
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),700);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),714.144);	
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1938.13);		
	
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


