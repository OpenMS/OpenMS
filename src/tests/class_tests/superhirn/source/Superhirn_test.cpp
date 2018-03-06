// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Florian Zeller$
// $Authors: Florian Zeller$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include "test_config.h"

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>
///////////////////////////

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/MzDataFile.h>

START_TEST(FeatureFinderAlgorithmSH, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmSH FFSH;

FFSH* ptr = nullptr;
FFSH* nullPtr = nullptr;

START_SECTION((FeatureFinderAlgorithmSH()))
  ptr = new FFSH;
  TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~FeatureFinderAlgorithmSH()))
  delete ptr;
END_SECTION

ptr = new FeatureFinderAlgorithmSH();

START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
  FeatureFinderAlgorithm* ptr2 = FFSH::create();
  TEST_NOT_EQUAL(ptr2,nullPtr)
  delete ptr2;
END_SECTION

START_SECTION((static const String getProductName()))
  TEST_EQUAL(FFSH::getProductName(),"superhirn")
END_SECTION

START_SECTION((virtual void run()))
  //input and output
  PeakMap input;
  MzDataFile mzdata_file;
  mzdata_file.getOptions().addMSLevel(1);
  mzdata_file.load(SUPERHIRN_GET_TEST_DATA_PATH("FeatureFinderAlgorithmSH_input.mzData"),input);

  input.updateRanges(1);
  FeatureMap output;

  //parameters
  Param param;
  //param.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"));
  param = param.copy("FeatureFinder:1:algorithm:",true);
  //Dummy featurefinder
  FeatureFinder ff;

  FFSH ffsh;
  ffsh.setParameters(param);
  ffsh.setData(input, output, ff);
  ffsh.run();

  TEST_EQUAL(output.size(), 384);

  //TOLERANCE_ABSOLUTE(0.001);
  //TEST_REAL_SIMILAR(output[0].getOverallQuality(),0.8819);
  //TEST_REAL_SIMILAR(output[1].getOverallQuality(),0.8673);
  // ...

  //TOLERANCE_ABSOLUTE(20.0);
  TEST_REAL_SIMILAR(output[0].getIntensity(),20829);
  TEST_REAL_SIMILAR(output[1].getIntensity(),56818.6);
  // ...

  TEST_REAL_SIMILAR(output[0].getMZ(),300.060882568359);
  TEST_REAL_SIMILAR(output[1].getMZ(),300.060882568359);

  TEST_REAL_SIMILAR(output[0].getRT(),35.1000317866759);
  TEST_REAL_SIMILAR(output[1].getRT(),134.37407934271);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
