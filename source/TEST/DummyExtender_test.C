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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummyExtender.h>

#include <OpenMS/FORMAT/MzDataFile.h>

#include<iostream>  
#include<fstream> 

using namespace OpenMS;

START_TEST(DummyExtender, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
	
using namespace OpenMS;
using namespace std;

// default ctor
DummyExtender* ptr = 0;
CHECK((DummyExtender()))
	ptr = new DummyExtender();
  TEST_EQUAL(ptr->getName(), "DummyExtender")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~DummyExtender()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(DummyExtender::getProductName(),"DummyExtender")
	TEST_EQUAL(DummyExtender().getName(),"DummyExtender")
RESULT

CHECK(static BaseExtender* create())
	TEST_NOT_EQUAL(DummyExtender::create(),0)
RESULT

CHECK(DummyExtender& operator=(const DummyExtender &rhs))
	DummyExtender ms1;
	DummyExtender ms2;
	
	ms1 = ms2;
	
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK(DummyExtender(const DummyExtender &rhs))
	DummyExtender ms1;
	DummyExtender ms2(ms1);
		
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK((const ChargedIndexSet& extend(const ChargedIndexSet &seed_region)))
	PRECISION(0.01)
	
  DummyExtender extender;
  FeaFiTraits* traits = new FeaFiTraits();
 
	MSExperiment< > exp;
	MzDataFile().load("data/DummyExtenderTestData.mzData",exp);
	
	traits->setData(exp.begin(), exp.end(),100);	
	extender.setTraits(traits);
	
	Param param;
	param.setValue("min_intensity_contribution",0.05);
	extender.setParameters(param);
		
	FeaFiModule::ChargedIndexSet  set;
	set.insert( std::make_pair(3,10) );		// point with max. ion count		
	
	// some other points around the maximum
	set.insert( std::make_pair(3,11) );		
	set.insert( std::make_pair(3,12) );	
	set.insert( std::make_pair(3,8) );		
	set.insert( std::make_pair(3,9) );	
	set.insert( std::make_pair(2,10) );		
	set.insert( std::make_pair(4,10) );	

	// extend seeding region	
	FeaFiModule::ChargedIndexSet region = extender.extend(set);
			
	ifstream infile( "data/DummyExtender_region1");
	
	DoubleReal intensity, rt, mz;
	
	FeaFiModule::ChargedIndexSet::const_iterator citer = region.begin();
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
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
