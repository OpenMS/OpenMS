// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Dominik Damerow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmWatershed.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

START_TEST(FeatureFinderAlgorithmWatershed, "$Id: FeatureFinderAlgorithmWatershed_test.C$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmWatershed<Peak1D,Feature> FFAW;

FFAW* ptr;
START_SECTION(FeatureFinderAlgorithmWatershed())
	ptr = new FFAW;
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((virtual ~FeatureFinderAlgorithmWatershed()))
	delete ptr;
END_SECTION

START_SECTION([EXTRA] FeatureFinderAlgorithmWatershed() - with RichPeak1D)
	FeatureFinderAlgorithmWatershed<RichPeak1D,Feature> ffa;
END_SECTION

START_SECTION(virtual void run())
  //create input
	MSExperiment<> input;
	input.reserve(500);
	for (UInt s=0; s<500; ++s)
	{
		MSExperiment<>::SpectrumType spec;
		spec.reserve(500);
		spec.setRT(s);
		spec.setMSLevel(1);
		for (UInt p=500; p<1000; ++p)
		{
			MSExperiment<>::PeakType peak;
			peak.setMZ(p);
			peak.setIntensity(1.0);
			spec.push_back(peak);
		}
		input.push_back(spec);
	}
	
	
	//create parameters
	Param param;
	param.setValue("debug","true");
	param.setValue("mz_sampling",1.0);
  param.setValue("cutoff_factor",7.0);
	
	//create dummy feature finder and run algorithm	
	FeatureMap<> output;
	FeatureFinder ff;
	FFAW ffaw;
	ffaw.setParameters(param);
	
	//--------------------------------------------------------------------
	//TEST WITH FLAT MAP
	input.updateRanges(1);
	ffaw.setData(input, output, ff);
  ffaw.run();
	TEST_EQUAL(output.size(),1);
	//MzDataFile().store(String("ff_in_")+test_number,input);
	//FeatureXMLFile().store(String("ff_out_")+(test_number++),output);
		
	//--------------------------------------------------------------------
	//TEST WITH ONE BASIN IN THE CENTER
	input[250][250].setIntensity(4.0);
	input.updateRanges(1);
	ffaw.setData(input, output, ff);
  ffaw.run();
	TEST_EQUAL(output.size(),1);
	TEST_REAL_SIMILAR(output[0].getRT(),250.0);
	TEST_REAL_SIMILAR(output[0].getMZ(),749.5);
	//MzDataFile().store(String("ff_in_")+test_number,input);
	//FeatureXMLFile().store(String("ff_out_")+(test_number++),output);
	
	//--------------------------------------------------------------------
	//TEST WITH THREE BASIN IN V-SHAPE
	input[125][125].setIntensity(2.0);
	input[125][375].setIntensity(4.0);
	input.updateRanges(1);
	ffaw.setData(input, output, ff);
  ffaw.run();
	TEST_EQUAL(output.size(),3);
	TEST_REAL_SIMILAR(output[0].getRT(),125.0);
	TEST_REAL_SIMILAR(output[0].getMZ(),874.5);
	TEST_REAL_SIMILAR(output[1].getRT(),250.0);
	TEST_REAL_SIMILAR(output[1].getMZ(),749.5);
	TEST_REAL_SIMILAR(output[2].getRT(),125.0);
	TEST_REAL_SIMILAR(output[2].getMZ(),625.5);

	//MzDataFile().store(String("ff_in_")+test_number,input);
	//FeatureXMLFile().store(String("ff_out_")+(test_number++),output);

END_SECTION

START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFAW::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(FFAW::getProductName(),"watershed")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
