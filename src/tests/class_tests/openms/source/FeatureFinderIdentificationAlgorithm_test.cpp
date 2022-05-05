// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderIdentificationAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderIdentificationAlgorithm* ptr = 0;
FeatureFinderIdentificationAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderIdentificationAlgorithm())
{
  ptr = new FeatureFinderIdentificationAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFinderIdentificationAlgorithm())
{
  delete ptr;
}
END_SECTION

// expose protected methods for testing
struct FFIdAlgoTest: public FeatureFinderIdentificationAlgorithm
{
  static auto extractTargetID(const Feature& feature, bool extract_charge = false)
  {
    return extractTargetID_(feature, extract_charge);
  }
};

START_SECTION((std::pair<String, Int> extractTargetID_(const Feature& feature, bool extract_charge = false)))
{
  Feature feat;
  feat.setMetaValue("PeptideRef", "PEP:XXXXX/2");
  auto pair1 = FFIdAlgoTest::extractTargetID(feat, false);
  TEST_EQUAL(pair1.first, "PEP:XXXXX");
  TEST_EQUAL(pair1.second, 0);
  auto pair2 = FFIdAlgoTest::extractTargetID(feat, true);
  TEST_EQUAL(pair2.first, "PEP:XXXXX");
  TEST_EQUAL(pair2.second, 2);
  feat.setMetaValue("PeptideRef", "PEP:XXXXX/2#1");
  auto pair3 = FFIdAlgoTest::extractTargetID(feat, true);
  TEST_EQUAL(pair3.first, "PEP:XXXXX");
  TEST_EQUAL(pair3.second, 2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
