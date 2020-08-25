// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderMultiplexAlgorithm, "$Id$")

FeatureFinderMultiplexAlgorithm* ptr = 0;
FeatureFinderMultiplexAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderMultiplexAlgorithm())
{
  ptr = new FeatureFinderMultiplexAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~FeatureFinderMultiplexAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual void run()))
{
  MzMLFile mzml_file;
  MSExperiment exp;
  ConsensusMap result;
  
  mzml_file.getOptions().addMSLevel(1);
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_input.mzML"), exp);
  exp.updateRanges(1);
  
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_parameters.ini"), param);
  param = param.copy("FeatureFinderMultiplex:1:",true);
  param.remove("in");
  param.remove("out");
  param.remove("out_multiplets");
  param.remove("log");
  param.remove("debug");
  param.remove("threads");
  param.remove("no_progress");
  param.remove("force");
  param.remove("test");
  
  FeatureFinderMultiplexAlgorithm algorithm;
  algorithm.setParameters(param);
  algorithm.run(exp, true);
  result = algorithm.getConsensusMap();
  
  TEST_EQUAL(result.size(), 2);
  
  double L = result[0].getFeatures().begin()->getIntensity();
  double H = (++(result[0].getFeatures().begin()))->getIntensity();

  // Check that the HEAVY:LIGHT ratio is close to the expected 3:1 ratio
  TOLERANCE_ABSOLUTE(0.2);
  TEST_REAL_SIMILAR(H/L, 3.0);
}
END_SECTION

END_TEST
