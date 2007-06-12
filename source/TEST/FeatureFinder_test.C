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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/Param.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinder, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinder* ptr = 0;
CHECK((FeatureFinder()))
	ptr = new FeatureFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~FeatureFinder()))
	delete ptr;
RESULT

CHECK((FeatureFinder& operator = (const FeatureFinder& source)))
  FeatureFinder ff1;
    
  ff1.addSeeder("SimpleSeeder");
  ff1.addExtender("SimpleExtender");
  ff1.addFitter("SimpleModelFitter"); 
  
  FeatureFinder ff2;
  ff2 = ff1;
  
  TEST_EQUAL(ff1==ff2,true);
  
RESULT

CHECK((FeatureFinder(const FeatureFinder& source)))

  FeatureFinder ff1;
    
  ff1.addSeeder("SimpleSeeder");
  ff1.addExtender("SimpleExtender");
  ff1.addFitter("SimpleModelFitter"); 
  
  FeatureFinder ff2(ff1);
    
  TEST_EQUAL(ff1==ff2,true);
  
RESULT

CHECK((bool operator==(const FeatureFinder &rhs) const))

	 // implicitely tested above already, but let's do it again (so Marc is happy)
	 FeatureFinder ff1;
   
  ff1.addSeeder("SimpleSeeder");
  ff1.addExtender("SimpleExtender");
  ff1.addFitter("SimpleModelFitter"); 
  
  FeatureFinder ff2(ff1);
  TEST_EQUAL(ff1==ff2,true);

RESULT

CHECK((bool setParam(const Param& param)))
	
	Param p;
	p.load("data/FeaFi_testparams.ini");
	
	FeatureFinder feafi;
	
	TEST_EQUAL(feafi.setParam(p),true);
		
RESULT

CHECK((const FeatureMap& run()))

	Param p;
	p.load("data/FeaFi_testparams.ini");
	
	FeatureFinder feafi;
	feafi.setParam(p);
	
	TEST_EQUAL(feafi.setParam(p),true);
	
	FeatureMap<> features = feafi.run();
	
	TEST_EQUAL(features.size(),0);
	
RESULT

CHECK((template <class SpectrumIteratorType> void setData(const SpectrumIteratorType &begin, const SpectrumIteratorType &end, UInt buffer_size)))
  	
	DPeakArray<Peak2D> parray;
	
	for (int i=0; i<10;++i)
	{
		Peak2D p1;
		p1.getPosition()[0] = i;
		p1.getPosition()[1] = i*5;
		p1.setIntensity(i*10);
		
		parray.push_back(p1);
	}
	
	FeatureFinder ff;
	
	MSExperimentExtern<Peak1D > exp;
	exp.set2DData(parray);
	ff.setData(exp.begin(), exp.end(),100);

RESULT

CHECK(void setLogType(LogType lg) const)
	////
RESULT

CHECK((void addSeeder(const String &name, const Param *param=0)))
	// not much happening here
RESULT

CHECK((void addExtender(const String &name, const Param *param=0)))
	// not much happening here
RESULT

CHECK((void addFitter(const String &name, const Param *param=0)))
	// not much happening here
RESULT

CHECK((void removeSeeder(const String &name)))
	// not much happening here
RESULT

CHECK((void removeExtender(const String &name)))
	// not much happening here
RESULT

CHECK((void removeFitter(const String &name)))
	// not much happening here
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



