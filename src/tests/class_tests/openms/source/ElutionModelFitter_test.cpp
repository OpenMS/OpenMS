// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
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
  ElutionModelFitter emf;

  FeatureMap features;
  // test if exception is thrown on empty featuremap
  TEST_EXCEPTION(Exception::MissingInformation, emf.fitElutionModels(features));  

  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionModelFitter_test.featureXML"), features);
  ABORT_IF(features.size() != 25);

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
