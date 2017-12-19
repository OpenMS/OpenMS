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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderAlgorithmMRM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderAlgorithmMRM* ptr = nullptr;
FeatureFinderAlgorithmMRM* nullPointer = nullptr;
FeatureFinderAlgorithm* ffA_nullPointer = nullptr;

START_SECTION(FeatureFinderAlgorithmMRM())
{
  ptr = new FeatureFinderAlgorithmMRM();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FeatureFinderAlgorithmMRM())
{
	delete ptr;
}
END_SECTION

ptr = new FeatureFinderAlgorithmMRM();

START_SECTION((virtual void run()))
{
	FeatureFinder ff;
  ff.setLogType(ProgressLogger::NONE);

  PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmMRM_input.mzML"), exp);

	FeatureMap features, seeds;
	Param ff_param(ptr->getParameters());
	ff.run("mrm", exp, features, ff_param, seeds);

	TEST_EQUAL(exp.getChromatograms().size(), 3)

	FeatureMap new_features;
	for (Size i = 0; i != features.size(); ++i)
	{
		if (features[i].getQuality(0) > 0.99)
		{
			new_features.push_back(features[i]);
		}
	}

	TEST_EQUAL(new_features.size(), 3)

	for (Size i = 0; i != new_features.size(); ++i)
	{
		TEST_EQUAL(new_features[i].getIntensity() > 100000, true)
	}
}
END_SECTION

START_SECTION((static FeatureFinderAlgorithm<PeakType>* create()))
{
  FeatureFinderAlgorithm* ptr2 = nullptr;
  ptr2 = FeatureFinderAlgorithmMRM::create();
  TEST_NOT_EQUAL(ptr2, ffA_nullPointer)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_STRING_EQUAL(ptr->getProductName(), "mrm")
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



