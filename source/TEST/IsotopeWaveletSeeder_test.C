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

CHECK((~IsotopeWaveletSeeder()))
	delete ptr;
RESULT

CHECK((FeaFiModule::IndexSet  nextSeed() throw(NoSuccessor)))
	PRECISION(0.01)
	
	// test single scan
  IsotopeWaveletSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment< > exp;
	MzDataFile().load("data/IsotopeWaveletTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);
	
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("rtvotes_cutoff",6);
	param.setValue("avg_intensity_factor",3.0);
	param.setValue("intensity_factor",2.0);
	param.setValue("scans_to_sumup",3);
	seeder.setParameters(param);
	
	FeaFiModule::IndexSet  region;
		
	region = seeder.nextSeed();
	
	TEST_EQUAL(region.size() == 263,true);
	
	FeaFiModule::IndexSet::const_iterator citer = region.begin();

	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),339292);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),1163.3);
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1402.77);
	
	++citer;
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),829624);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),1163.4);
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1402.77);
	
	// test point with highest intensity
	FeaFiModule::IDX max_peak;
	double max_intensity = 0.0;
	
	for (citer = region.begin(); citer != region.end(); ++citer)
	{
			if (traits->getPeakIntensity(*citer) > max_intensity)
			{
				max_peak      = *citer;
				max_intensity = traits->getPeakIntensity(*citer);
			}	
	}

	{
	// we test absolute precision only, so here we allow a lower precision
	PRECISION(10);		
	TEST_REAL_EQUAL(traits->getPeakIntensity(max_peak),2.10752e+06);	
	}
	
	TEST_REAL_EQUAL(traits->getPeakMz(max_peak),1163.6);
	TEST_REAL_EQUAL(traits->getPeakRt(max_peak),1408.41);
	
	// test last two peaks
	citer = region.end();
	--citer;
	--citer;
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),3061);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),1166.1);
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1414.08);
	
	++citer;
	
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),7342);	
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),1166.2);
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),1414.08);
	
	
	// test exception, there should be no more seeds
	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )
 
RESULT


CHECK(static const String getProductName())
	TEST_EQUAL(IsotopeWaveletSeeder::getProductName(),"IsotopeWaveletSeeder");
RESULT

CHECK(static BaseSeeder* create())
	BaseSeeder* base = IsotopeWaveletSeeder::create();
	TEST_NOT_EQUAL(base,0);
	delete(base);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


