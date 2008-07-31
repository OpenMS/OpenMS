// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmPicked, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmPicked<Peak1D,Feature> FFPP;

FFPP* ptr;
CHECK(FeatureFinderAlgorithmPicked())
	ptr = new FFPP;
	TEST_NOT_EQUAL(ptr,0)
	delete ptr;
RESULT

CHECK((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFPP::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(FFPP::getProductName(),"picked_peak")
RESULT

CHECK(virtual void run())
	//input and output
	MSExperiment<> input;
	MzDataFile().load("data/FeatureFinderAlgorithmPicked.mzData",input);
	input.updateRanges(1);
	FeatureMap<> output;
	
	//parameters
	Param param;
	param.load("data/FeatureFinderAlgorithmPicked.ini");
	param = param.copy("FeatureFinder:1:algorithm:",true);
	//Dummy featurefinder
	FeatureFinder ff;
	
	FFPP ffpp;
	ffpp.setParameters(param);
	ffpp.setData(input, output, ff);
	ffpp.run();
	
	TEST_EQUAL(output.size(),8);
	
	PRECISION(0.001);
	TEST_REAL_EQUAL(output[0].getOverallQuality(),0.8819);
	TEST_REAL_EQUAL(output[1].getOverallQuality(),0.8674);
	TEST_REAL_EQUAL(output[2].getOverallQuality(),0.9083);
	TEST_REAL_EQUAL(output[3].getOverallQuality(),0.9268);
	TEST_REAL_EQUAL(output[4].getOverallQuality(),0.9402);
	TEST_REAL_EQUAL(output[5].getOverallQuality(),0.9093);
	TEST_REAL_EQUAL(output[6].getOverallQuality(),0.9403);
	TEST_REAL_EQUAL(output[7].getOverallQuality(),0.9243);
	
	PRECISION(20.0);
	TEST_REAL_EQUAL(output[0].getIntensity(),51249.6);
	TEST_REAL_EQUAL(output[1].getIntensity(),44637.9);
	TEST_REAL_EQUAL(output[2].getIntensity(),34596.9);
	TEST_REAL_EQUAL(output[3].getIntensity(),19423.1);
	TEST_REAL_EQUAL(output[4].getIntensity(),12528.0);
	TEST_REAL_EQUAL(output[5].getIntensity(),8510.74);
	TEST_REAL_EQUAL(output[6].getIntensity(),7295.91);
	TEST_REAL_EQUAL(output[7].getIntensity(),5026.28);
RESULT

//remove log file
File::remove("featurefinder.log");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
