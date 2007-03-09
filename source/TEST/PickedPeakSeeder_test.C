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

//   PickedPeakSeeder seeder;
//   FeaFiTraits* traits = new FeaFiTraits();
//  
// 	MSExperiment<Peak1D > exp;
// 	MzDataFile().load("data/PickedPeakTestData.mzData",exp);
	
// 	traits->setData(exp.begin(), exp.end(),100);	
// 	seeder.setTraits(traits);
// 	
// 	Param param;
//   param.setValue("min_number_scans",8);
// 	seeder.setParameters(param);
// 	
// 	// test first region
// 	FeaFiModule::IndexSet  region = seeder.nextSeed();	
// 	
// 	String fname("ppseeder_region");
// 	ofstream out( fname.c_str() );
// 	
// 	for (FeaFiModule::IndexSet::const_iterator citer = region.begin(); citer != region.end();++citer)
// 	{
// 	out << traits->getPeakRt(*citer) << " ";
// 	out << traits->getPeakMz(*citer) << " ";
// 	out << traits->getPeakIntensity(*citer) << endl;
// 	}
// 
// 	out.close();
// 	
// 	region = seeder.nextSeed();	
// 	
// 	String fname2("ppseeder_region2");
// 	ofstream out2( fname2.c_str() );
// 	
// 	for (FeaFiModule::IndexSet::const_iterator citer = region.begin(); citer != region.end();++citer)
// 	{
// 	out2 << traits->getPeakRt(*citer) << " ";
// 	out2 << traits->getPeakMz(*citer) << " ";
// 	out2 << traits->getPeakIntensity(*citer) << endl;
// 	}
// 	out2.close();
// 	
// 	region = seeder.nextSeed();	
// 	
// 	String fname3("ppseeder_region3");
// 	ofstream out3( fname3.c_str() );
// 	
// 	for (FeaFiModule::IndexSet::const_iterator citer = region.begin(); citer != region.end();++citer)
// 	{
// 	out3 << traits->getPeakRt(*citer) << " ";
// 	out3 << traits->getPeakMz(*citer) << " ";
// 	out3 << traits->getPeakIntensity(*citer) << endl;
// 	}
// 	out3.close();


// 	TEST_EXCEPTION( FeaFiModule::NoSuccessor , seeder.nextSeed() )

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


