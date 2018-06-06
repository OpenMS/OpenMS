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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ElutionModelFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ElutionModelFitter* ptr = nullptr;
ElutionModelFitter* null_ptr = nullptr;

START_SECTION((ElutionModelFitter()))
{
  ptr = new ElutionModelFitter();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION((~ElutionModelFitter()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void fitElutionModels(FeatureMap& features)))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionModelFitter_test.featureXML"), features);
  ABORT_IF(features.size() != 25);

  ElutionModelFitter emf;

  // symmetric model (default):
  emf.fitElutionModels(features);
  TEST_EQUAL(features.size(), 25);
  for (FeatureMap::ConstIterator it = features.begin(); it != features.end();
       ++it)
  {
    TEST_EQUAL(it->metaValueExists("model_area"), true);
    TEST_EQUAL(it->metaValueExists("model_status"), true);
    TEST_EQUAL(it->metaValueExists("raw_intensity"), true);
    TEST_NOT_EQUAL(it->getIntensity(), it->getMetaValue("raw_intensity"));
    TEST_EQUAL(it->metaValueExists("model_Gauss_sigma"), true);
  }

  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionModelFitter_test.featureXML"), features);
  ABORT_IF(features.size() != 25);

  // asymmetric model:
  Param params;
  params.setValue("asymmetric", "true");
  emf.setParameters(params);
  emf.fitElutionModels(features);
  TEST_EQUAL(features.size(), 25);
  for (FeatureMap::ConstIterator it = features.begin(); it != features.end();
       ++it)
  {
    TEST_EQUAL(it->metaValueExists("model_area"), true);
    TEST_EQUAL(it->metaValueExists("model_status"), true);
    TEST_EQUAL(it->metaValueExists("raw_intensity"), true);
    TEST_NOT_EQUAL(it->getIntensity(), it->getMetaValue("raw_intensity"));
    TEST_EQUAL(it->metaValueExists("model_EGH_tau"), true);
    TEST_EQUAL(it->metaValueExists("model_EGH_sigma"), true);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
