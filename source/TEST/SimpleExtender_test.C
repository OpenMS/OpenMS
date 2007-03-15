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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/DPeak.h>

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(SimpleExtender, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
SimpleExtender* ptr = 0;
CHECK((SimpleExtender()))
ptr = new SimpleExtender();
TEST_EQUAL(ptr->getName(), "SimpleExtender")
TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~SimpleExtender()))
delete ptr;
RESULT

CHECK((const IndexSet& extend(const IndexSet &seed_region)))

// this test checks the regions returned by SimpleExtender
// on one artificial data set and a picked (centroided) data set
// The test of the corresponding TOPP module performs further tests.

SimpleExtender extender;
FeaFiTraits* traits = new FeaFiTraits();
MSExperiment<Peak1D > exp;

MSExperiment<Peak1D >::SpectrumType spec;
spec.setRetentionTime(1);

double mzs[] = {1, 2, 3, 4, 5};
double its1[] = {1000, 1500,2000, 1500, 1000};

const UInt num = 5;

for (UInt i=0; i < num; i++)
{
    Peak1D p;
    p.setMZ(mzs[i]);
    p.setIntensity(its1[i]);

    spec.push_back(p);
}
exp.push_back(spec);

spec.clear();
spec.setRetentionTime(2);

double its2[] = {1000, 1500,2000, 1500, 1000};

for (UInt i=0; i < num; i++)
{
    Peak1D p;
    p.setMZ(mzs[i]);
    p.setIntensity(its2[i]);

    spec.push_back(p);
}
exp.push_back(spec);

spec.clear();
spec.setRetentionTime(3);

double its3[] = {1000, 1500,5000, 1500, 1000};

for (UInt i=0; i < num; i++)
{
    Peak1D p;
    p.setMZ(mzs[i]);
    p.setIntensity(its3[i]);

    spec.push_back(p);
}
exp.push_back(spec);

spec.clear();
spec.setRetentionTime(4);

// the last two data points should not be included (see param intensity_factor)
double its4[] = {1000, 1500,2000, 0.1, 0.1};

for (UInt i=0; i < num; i++)
{
    Peak1D p;
    p.setMZ(mzs[i]);
    p.setIntensity(its4[i]);

    spec.push_back(p);
}
exp.push_back(spec);

traits->setData(exp.begin(), exp.end(),100);

// first two points are already in some other feature region
// => check if points included in other features are ignored
traits->getPeakFlag(make_pair(0,0)) = FeaFiTraits::INSIDE_FEATURE;
traits->getPeakFlag(make_pair(0,1)) = FeaFiTraits::INSIDE_FEATURE;

extender.setTraits(traits);

FeaFiModule::IndexSet  set;

set.insert( std::make_pair(2,2) );		// start extension with point of highest intensity

FeaFiModule::IndexSet region = extender.extend(set);

// 20 points in total, 2 in another feature, 2 with too little intensity => region should be of size 16
TEST_EQUAL(region.size(),16);

FeaFiModule::IndexSet::const_iterator citer = region.begin();

// first scan
TEST_EQUAL(traits->getPeakIntensity(*citer),2000.0);
TEST_EQUAL(traits->getPeakMz(*citer),3.0);
TEST_EQUAL(traits->getPeakRt(*citer),1.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),4.0);
TEST_EQUAL(traits->getPeakRt(*citer),1.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),5.0);
TEST_EQUAL(traits->getPeakRt(*citer),1.0);

// second scan
++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),1.0);
TEST_EQUAL(traits->getPeakRt(*citer),2.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),2.0);
TEST_EQUAL(traits->getPeakRt(*citer),2.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),2000.0);
TEST_EQUAL(traits->getPeakMz(*citer),3.0);
TEST_EQUAL(traits->getPeakRt(*citer),2.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),4.0);
TEST_EQUAL(traits->getPeakRt(*citer),2.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),5.0);
TEST_EQUAL(traits->getPeakRt(*citer),2.0);


// third scan
++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),1.0);
TEST_EQUAL(traits->getPeakRt(*citer),3.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),2.0);
TEST_EQUAL(traits->getPeakRt(*citer),3.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),5000.0);
TEST_EQUAL(traits->getPeakMz(*citer),3.0);
TEST_EQUAL(traits->getPeakRt(*citer),3.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),4.0);
TEST_EQUAL(traits->getPeakRt(*citer),3.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),5.0);
TEST_EQUAL(traits->getPeakRt(*citer),3.0);


// fourth scan
++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1000.0);
TEST_EQUAL(traits->getPeakMz(*citer),1.0);
TEST_EQUAL(traits->getPeakRt(*citer),4.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),1500.0);
TEST_EQUAL(traits->getPeakMz(*citer),2.0);
TEST_EQUAL(traits->getPeakRt(*citer),4.0);

++citer;
TEST_EQUAL(traits->getPeakIntensity(*citer),2000.0);
TEST_EQUAL(traits->getPeakMz(*citer),3.0);
TEST_EQUAL(traits->getPeakRt(*citer),4.0);

// last two points have too little intensity
++ citer;

TEST_EQUAL(citer==region.end(),true)

RESULT

CHECK(([EXTRA] Extension on real-world data))

PRECISION(0.01)

SimpleExtender extender;
FeaFiTraits* traits = new FeaFiTraits();
MSExperiment<Peak1D > exp;

MzDataFile().load("data/SimpleExtender_test.mzData",exp);
traits->setData(exp.begin(),exp.end(),100);
extender.setTraits(traits);

Param param;
param.setValue("intensity_factor",0.01);
extender.setParameters(param);

FeaFiModule::IndexSet  set;
// SimpleExtender starts at maximum point
set.insert( std::make_pair(15,15) );	

FeaFiModule::IndexSet region = extender.extend(set);
		
ifstream infile( "data/SimpleExtender_region1");
	
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

RESULT

CHECK(([EXTRA] Extension on picked data))

PRECISION(0.01)

SimpleExtender extender;
FeaFiTraits* traits = new FeaFiTraits();
MSExperiment<Peak1D > exp;

MzDataFile().load("data/SimpleExtender_test2.mzData",exp);
traits->setData(exp.begin(),exp.end(),100);
extender.setTraits(traits);

FeaFiModule::IndexSet  set;
// SimpleExtender starts at maximum point
set.insert( std::make_pair(2,42) );	

FeaFiModule::IndexSet region = extender.extend(set);
		
ifstream infile( "data/SimpleExtender_region2");
	
DoubleReal intensity, rt, mz;
	
FeaFiModule::IndexSet::const_iterator citer = region.begin();
while ( infile >> rt )
{
	infile >> mz >> intensity;
		
	TEST_NOT_EQUAL(citer == region.end(),true)
		
	TEST_REAL_EQUAL(traits->getPeakRt(*citer),rt)
	TEST_REAL_EQUAL(traits->getPeakMz(*citer),mz)
	{
	PRECISION(1000)	// lower (absolute) precision for high intensities
	TEST_REAL_EQUAL(traits->getPeakIntensity(*citer),intensity)
	}		
	++citer;				
}	
infile.close();

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

