// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Oliver Kohlbacher $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/KERNEL/RichPeak1D.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmPicked, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmPicked<Peak1D,Feature> FFPP;

FFPP* ptr;
START_SECTION((FeatureFinderAlgorithmPicked()))
	ptr = new FFPP;
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((~FeatureFinderAlgorithmPicked()))
	delete ptr;
END_SECTION

START_SECTION([EXTRA] FeatureFinderAlgorithmPicked() - with RichPeak1D)
	FeatureFinderAlgorithmPicked<RichPeak1D,Feature> ffa;
	NOT_TESTABLE
END_SECTION
	
START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFPP::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(FFPP::getProductName(),"centroided")
END_SECTION

START_SECTION((virtual void run()))
	//input and output
	MSExperiment<> input;
	MzDataFile mzdata_file;
	mzdata_file.getOptions().addMSLevel(1);
	mzdata_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.mzData"),input);
	input.updateRanges(1);
	FeatureMap<> output;
	
	//parameters
	Param param;
	param.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"));
	param = param.copy("FeatureFinder:1:algorithm:",true);
	//Dummy featurefinder
	FeatureFinder ff;
	
	FFPP ffpp;
	ffpp.setParameters(param);
	ffpp.setData(input, output, ff);
	ffpp.run();
	
	TEST_EQUAL(output.size(),8);
	
	TOLERANCE_ABSOLUTE(0.001);
	TEST_REAL_SIMILAR(output[0].getOverallQuality(),0.8819);
	TEST_REAL_SIMILAR(output[1].getOverallQuality(),0.8673);
	TEST_REAL_SIMILAR(output[2].getOverallQuality(),0.9079);
	TEST_REAL_SIMILAR(output[3].getOverallQuality(),0.9271);
	TEST_REAL_SIMILAR(output[4].getOverallQuality(),0.9401);
	TEST_REAL_SIMILAR(output[5].getOverallQuality(),0.9094);
	TEST_REAL_SIMILAR(output[6].getOverallQuality(),0.9403);
	TEST_REAL_SIMILAR(output[7].getOverallQuality(),0.9243);
	
	TOLERANCE_ABSOLUTE(20.0);
	TEST_REAL_SIMILAR(output[0].getIntensity(),51260.0);
	TEST_REAL_SIMILAR(output[1].getIntensity(),44667.3);
	TEST_REAL_SIMILAR(output[2].getIntensity(),34613.3);
	TEST_REAL_SIMILAR(output[3].getIntensity(),19428.9);
	TEST_REAL_SIMILAR(output[4].getIntensity(),12513.9);
	TEST_REAL_SIMILAR(output[5].getIntensity(),8512.71);
	TEST_REAL_SIMILAR(output[6].getIntensity(),7295.3);
	TEST_REAL_SIMILAR(output[7].getIntensity(),5024.74);
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
