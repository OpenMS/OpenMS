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
// $Maintainer: Florian Zeller$
// $Authors: Florian Zeller$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/MzDataFile.h>
//#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmSH, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmSH<Peak1D,Feature> FFSH;

FFSH* ptr;
START_SECTION((FeatureFinderAlgorithmSH()))
	ptr = new FFSH;
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((~FeatureFinderAlgorithmSH()))
	delete ptr;
END_SECTION

ptr = new FeatureFinderAlgorithmSH<Peak1D, Feature>();

START_SECTION([EXTRA] FeatureFinderAlgorithmSH() - with RichPeak1D)
	FeatureFinderAlgorithmSH<RichPeak1D,Feature> ffa;
	NOT_TESTABLE
END_SECTION
	
START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFSH::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(FFSH::getProductName(),"superhirn")
END_SECTION

START_SECTION((virtual void run()))
  //input and output
  MSExperiment<> input;
  MzDataFile mzdata_file;
  mzdata_file.getOptions().addMSLevel(1);
  mzdata_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmSH_input.mzData"),input);
  input.updateRanges(1);
  FeatureMap<> output;

  //parameters
  Param param;
  //param.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"));
  //param = param.copy("FeatureFinder:1:algorithm:",true);
  //Dummy featurefinder
  FeatureFinder ff;

  FFSH ffsh;
  ffsh.setParameters(param);
  ffsh.setData(input, output, ff);
  ffsh.run();
	
	//TOLERANCE_ABSOLUTE(0.001);
	//TEST_REAL_SIMILAR(output[0].getOverallQuality(),0.8819);
	//TEST_REAL_SIMILAR(output[1].getOverallQuality(),0.8673);
	//TEST_REAL_SIMILAR(output[2].getOverallQuality(),0.9079);
	//TEST_REAL_SIMILAR(output[3].getOverallQuality(),0.9271);
	//TEST_REAL_SIMILAR(output[4].getOverallQuality(),0.9401);
	//TEST_REAL_SIMILAR(output[5].getOverallQuality(),0.9094);
	//TEST_REAL_SIMILAR(output[6].getOverallQuality(),0.9403);
	//TEST_REAL_SIMILAR(output[7].getOverallQuality(),0.9243);
	
	//TOLERANCE_ABSOLUTE(20.0);
	//TEST_REAL_SIMILAR(output[0].getIntensity(),51260.0);
	//TEST_REAL_SIMILAR(output[1].getIntensity(),44667.3);
	//TEST_REAL_SIMILAR(output[2].getIntensity(),34613.3);
	//TEST_REAL_SIMILAR(output[3].getIntensity(),19428.9);
	//TEST_REAL_SIMILAR(output[4].getIntensity(),12513.9);
	//TEST_REAL_SIMILAR(output[5].getIntensity(),8512.71);
	//TEST_REAL_SIMILAR(output[6].getIntensity(),7295.3);
	//TEST_REAL_SIMILAR(output[7].getIntensity(),5024.74);
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
