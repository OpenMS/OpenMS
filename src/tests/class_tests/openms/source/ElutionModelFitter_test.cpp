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

START_SECTION(([EXTRA] void fitElutionModels(FeatureMap& features) - no failed fit to impute))
{
  // valid models of the same raw intensity (a linear model cannot be fitted to them), but no failed fit: nothing to
  // impute, so no imputation model is needed
  FeatureMap all;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionModelFitter_test.featureXML"), all);
  FeatureMap two;
  two.push_back(all[0]);
  two.push_back(all[1]);
  for (Feature& feature : two) feature.setIntensity(1.0e5);
  ElutionModelFitter emf;
  emf.fitElutionModels(two);
  ABORT_IF(two.size() != 2);
  TEST_EQUAL(two[0].getMetaValue("model_status"), "0 (valid)");
  TEST_EQUAL(two[1].getMetaValue("model_status"), "0 (valid)");
  TEST_NOT_EQUAL(two[0].getIntensity(), 1.0e5); // the model areas
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
