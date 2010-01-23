// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

START_TEST(SimpleExtender, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef Peak1D PeakType;
typedef Feature FeatureType;
typedef SimpleExtender<PeakType,FeatureType> ExtenderType;
typedef FeatureFinderDefs::ChargedIndexSet ChargedIndexSet;

{
	ExtenderType* ptr = 0;

	START_SECTION((SimpleExtender(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)))
		{
			MSExperiment<PeakType> map;
			FeatureMap<FeatureType> features;
			FeatureFinder ff;
			ptr = new ExtenderType(&map,&features,&ff);
			TEST_EQUAL(ptr->getName(), "SimpleExtender");
			TEST_NOT_EQUAL(ptr, 0);
		}
END_SECTION
		
		START_SECTION((virtual ~SimpleExtender()))
		delete ptr;
END_SECTION;
}
				
START_SECTION(void extend(const ChargedIndexSet &seed_region, ChargedIndexSet& result_region))
{
	// this test checks the regions returned by SimpleExtender
	// on one artificial data set and a picked (centroided) data set
	// The test of the corresponding TOPP module performs further tests.

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	MSExperiment<PeakType>::SpectrumType spec;
	spec.setRT(1);

	DoubleReal mzs[] = {1, 2, 3, 4, 5};
	DoubleReal its1[] = {1000, 1500,2000, 1500, 1000};

	const UInt num = 5;

	for (Size i=0; i < num; i++)
	{
		PeakType p;
		p.setMZ(mzs[i]);
		p.setIntensity(its1[i]);

		spec.push_back(p);
	}
	input.push_back(spec);

	spec.clear(true);
	spec.setRT(2);

	DoubleReal its2[] = {1000, 1500, 2000, 1500, 1000};

	for (Size i=0; i < num; i++)
	{
		PeakType p;
		p.setMZ(mzs[i]);
		p.setIntensity(its2[i]);

		spec.push_back(p);
	}
	input.push_back(spec);

	spec.clear(true);
	spec.setRT(3);

	DoubleReal its3[] = {1000, 1500, 5000, 1500, 1000};

	for (Size i=0; i < num; i++)
	{
		PeakType p;
		p.setMZ(mzs[i]);
		p.setIntensity(its3[i]);

		spec.push_back(p);
	}
	input.push_back(spec);

	spec.clear(true);
	spec.setRT(4);

	// the last two data points should not be included (see param intensity_factor)
	DoubleReal its4[] = {1000, 1500, 2000, 0.1, 0.1};

	for (Size i=0; i < num; i++)
	{
		PeakType p;
		p.setMZ(mzs[i]);
		p.setIntensity(its4[i]);

		spec.push_back(p);
	}
	input.push_back(spec);

	input.updateRanges();

  ff.run("none", input, features, Param(), FeatureMap<>());

	// first two points are already in some other feature region
	// => check if points included in other features are ignored
	ff.getPeakFlag(make_pair(0,0)) = FeatureFinderDefs::USED;
	ff.getPeakFlag(make_pair(0,1)) = FeatureFinderDefs::USED;

	ExtenderType extender(&input, &features, &ff);

	ChargedIndexSet index_set;

	index_set.insert( std::make_pair(2,2) );		// start extension with point of highest intensity

	ChargedIndexSet region;
  extender.extend(index_set,region);

	// We have 20 points in total, 2 in another feature, 2 with too little
	// intensity.  Hence the region should be of size 16.
	TEST_EQUAL(region.size(),16);

	ChargedIndexSet::const_iterator citer = region.begin();

	// first scan
	TEST_EQUAL(extender.getPeakIntensity(*citer),2000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),3.0);
	TEST_EQUAL(extender.getPeakRt(*citer),1.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),4.0);
	TEST_EQUAL(extender.getPeakRt(*citer),1.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),5.0);
	TEST_EQUAL(extender.getPeakRt(*citer),1.0);

	// second scan
	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),1.0);
	TEST_EQUAL(extender.getPeakRt(*citer),2.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),2.0);
	TEST_EQUAL(extender.getPeakRt(*citer),2.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),2000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),3.0);
	TEST_EQUAL(extender.getPeakRt(*citer),2.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),4.0);
	TEST_EQUAL(extender.getPeakRt(*citer),2.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),5.0);
	TEST_EQUAL(extender.getPeakRt(*citer),2.0);


	// third scan
	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),1.0);
	TEST_EQUAL(extender.getPeakRt(*citer),3.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),2.0);
	TEST_EQUAL(extender.getPeakRt(*citer),3.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),5000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),3.0);
	TEST_EQUAL(extender.getPeakRt(*citer),3.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),4.0);
	TEST_EQUAL(extender.getPeakRt(*citer),3.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),5.0);
	TEST_EQUAL(extender.getPeakRt(*citer),3.0);


	// fourth scan
	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),1.0);
	TEST_EQUAL(extender.getPeakRt(*citer),4.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),1500.0);
	TEST_EQUAL(extender.getPeakMz(*citer),2.0);
	TEST_EQUAL(extender.getPeakRt(*citer),4.0);

	++citer;
	TEST_EQUAL(extender.getPeakIntensity(*citer),2000.0);
	TEST_EQUAL(extender.getPeakMz(*citer),3.0);
	TEST_EQUAL(extender.getPeakRt(*citer),4.0);

	// last two points have too little intensity
	++ citer;

	TEST_EQUAL(citer==region.end(),true)
}
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(([EXTRA] Extension on real-world data))
{

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleExtender_test.mzData"),input);

	ExtenderType extender(&input, &features, &ff);

	Param param;
	param.setValue("intensity_factor",0.01);
	extender.setParameters(param);

	input.updateRanges();

  ff.run("none", input, features, Param(), FeatureMap<>());

	ChargedIndexSet  set;
	// SimpleExtender starts at maximum point
	set.insert( std::make_pair(15,15) );	

	ChargedIndexSet region;
  extender.extend(set,region);
		
	ifstream infile( OPENMS_GET_TEST_DATA_PATH("SimpleExtender_region1"));
	
	DoubleReal intensity, rt, mz;
	
	ChargedIndexSet::const_iterator citer = region.begin();
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		
			TEST_REAL_SIMILAR(extender.getPeakRt(*citer),rt)
			TEST_REAL_SIMILAR(extender.getPeakMz(*citer),mz)
			TEST_REAL_SIMILAR(extender.getPeakIntensity(*citer),intensity)
				
			++citer;				
	}	
	infile.close();
	TEST_EQUAL(citer == region.end(),true)

}
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(([EXTRA] Extension on picked data))
{

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	ExtenderType extender(&input, &features, &ff);

	MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleExtender_test2.mzData"),input);

	input.updateRanges();

  ff.run("none", input, features, Param(), FeatureMap<>());

	ChargedIndexSet  set;
	// SimpleExtender starts at maximum point
	set.insert( std::make_pair(2,42) );	

	ChargedIndexSet region;
	extender.extend(set,region);
	
	ifstream infile( OPENMS_GET_TEST_DATA_PATH("SimpleExtender_region2"));
	
	DoubleReal intensity, rt, mz;
	
	ChargedIndexSet::const_iterator citer = region.begin();
	while ( infile >> rt )
	{
		infile >> mz >> intensity;
		
		TEST_NOT_EQUAL(citer == region.end(),true)
		
			TEST_REAL_SIMILAR(extender.getPeakRt(*citer),rt)
			TEST_REAL_SIMILAR(extender.getPeakMz(*citer),mz)
			{
				TOLERANCE_ABSOLUTE(1000)	// lower (absolute) precision for high intensities
					TEST_REAL_SIMILAR(extender.getPeakIntensity(*citer),intensity)
					}		
		++citer;				
	}	
	infile.close();

	TEST_EQUAL(citer == region.end(),true)

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

