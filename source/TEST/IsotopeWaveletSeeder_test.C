// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(IsotopeWaveletSeeder, "$Id: IsotopeWaveletSeeder_test.C 921 2006-11-26 20:00:44Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

enum DimensionId
{
	RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
	MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
};


// default ctor
IsotopeWaveletSeeder* ptr = 0;
CHECK(IsotopeWaveletSeeder())
	ptr = new IsotopeWaveletSeeder();
  TEST_EQUAL(ptr->getName(), "IsotopeWaveletSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~IsotopeWaveletSeeder())
	delete ptr;
RESULT

CHECK(nextSeed())
  IsotopeWaveletSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment<DPeak<1> > exp;
	MzDataFile().load("data/IsotopeWaveletTestData.mzData",exp);
	
	traits->setData(exp);
	
	seeder.setTraits(traits);
	
	Param param;
  param.setValue("rtvotes_cutoff",0);
	seeder.setParam(param);
	
	IndexSet region;
	Index peak;
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),1048);	
	
	region = seeder.nextSeed();
	peak =  *(region.begin());
	TEST_EQUAL(traits->getPeakIntensity(peak),2057);					
  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


