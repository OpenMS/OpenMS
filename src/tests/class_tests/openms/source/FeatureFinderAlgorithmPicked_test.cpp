// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher, Stephan Aiche $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
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

FFPP* ptr = 0;
FFPP* nullPointer = 0;
FeatureFinderAlgorithm<Peak1D,Feature>* ffA_nullPointer = 0;

START_SECTION((FeatureFinderAlgorithmPicked()))
	ptr = new FFPP;
  TEST_NOT_EQUAL(ptr,nullPointer)
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
  TEST_NOT_EQUAL(ptr2,ffA_nullPointer)
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
  ParamXMLFile paramFile;
	paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"), param);
	param = param.copy("FeatureFinder:1:algorithm:",true);
	//Dummy featurefinder
	FeatureFinder ff;
	
	FFPP ffpp;
	ffpp.setParameters(param);
	ffpp.setData(input, output, ff);
	ffpp.run();
	
	TEST_EQUAL(output.size(), 8);
	
	TOLERANCE_ABSOLUTE(0.001);
	TEST_REAL_SIMILAR(output[0].getOverallQuality(), 0.8819);
	TEST_REAL_SIMILAR(output[1].getOverallQuality(), 0.8673);
	TEST_REAL_SIMILAR(output[2].getOverallQuality(), 0.9079);
	TEST_REAL_SIMILAR(output[3].getOverallQuality(), 0.9271);
	TEST_REAL_SIMILAR(output[4].getOverallQuality(), 0.9401);
	TEST_REAL_SIMILAR(output[5].getOverallQuality(), 0.9094);
	TEST_REAL_SIMILAR(output[6].getOverallQuality(), 0.9403);
	TEST_REAL_SIMILAR(output[7].getOverallQuality(), 0.9243);
	
	TOLERANCE_ABSOLUTE(20.0);
	TEST_REAL_SIMILAR(output[0].getIntensity(), 51393.3);
	TEST_REAL_SIMILAR(output[1].getIntensity(), 44786.0);
	TEST_REAL_SIMILAR(output[2].getIntensity(), 34704.7);
	TEST_REAL_SIMILAR(output[3].getIntensity(), 19479.9);
	TEST_REAL_SIMILAR(output[4].getIntensity(), 12565.2);
	TEST_REAL_SIMILAR(output[5].getIntensity(), 8535.44);
	TEST_REAL_SIMILAR(output[6].getIntensity(), 7315.13);
	TEST_REAL_SIMILAR(output[7].getIntensity(), 5041.12);
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
