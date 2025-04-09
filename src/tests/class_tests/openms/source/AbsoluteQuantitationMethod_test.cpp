// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AbsoluteQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethod* ptr = nullptr;
AbsoluteQuantitationMethod* nullPointer = nullptr;

START_SECTION(AbsoluteQuantitationMethod())
{
  ptr = new AbsoluteQuantitationMethod();
  TEST_NOT_EQUAL(ptr, nullPointer);

  TEST_EQUAL(ptr->getComponentName().size(), 0)
  TEST_EQUAL(ptr->getFeatureName().size(), 0)
  TEST_EQUAL(ptr->getISName().size(), 0)
  TEST_REAL_SIMILAR(ptr->getLLOD(), 0.0)
  TEST_REAL_SIMILAR(ptr->getULOD(), 0.0)
  TEST_REAL_SIMILAR(ptr->getLLOQ(), 0.0)
  TEST_REAL_SIMILAR(ptr->getULOQ(), 0.0)
  TEST_EQUAL(ptr->getNPoints(), 0)
  TEST_REAL_SIMILAR(ptr->getCorrelationCoefficient(), 0.0)
  TEST_EQUAL(ptr->getConcentrationUnits().size(), 0)
  TEST_EQUAL(ptr->getTransformationModel().size(), 0)
  TEST_EQUAL(ptr->getTransformationModelParams().size(), 0)
}
END_SECTION

START_SECTION(~AbsoluteQuantitationMethod())
{
  delete ptr;
}
END_SECTION

START_SECTION(all setters and getters)
{
  AbsoluteQuantitationMethod aqm;

  aqm.setComponentName("component");
  aqm.setFeatureName("feature");
  aqm.setISName("IS");
  aqm.setLLOD(1.2);
  aqm.setULOD(3.4);
  aqm.setLLOQ(5.6);
  aqm.setULOQ(7.8);
  aqm.setNPoints(9);
  aqm.setCorrelationCoefficient(0.44);
  aqm.setConcentrationUnits("uM");
  aqm.setTransformationModel("TransformationModelLinear");
  Param params1;
  params1.setValue("slope", 1);
  aqm.setTransformationModelParams(params1);

  TEST_EQUAL(aqm.getComponentName(), "component")
  TEST_EQUAL(aqm.getFeatureName(), "feature")
  TEST_EQUAL(aqm.getISName(), "IS")
  TEST_REAL_SIMILAR(aqm.getLLOD(), 1.2)
  TEST_REAL_SIMILAR(aqm.getULOD(), 3.4)
  TEST_REAL_SIMILAR(aqm.getLLOQ(), 5.6)
  TEST_REAL_SIMILAR(aqm.getULOQ(), 7.8)
  TEST_EQUAL(aqm.getNPoints(), 9)
  TEST_REAL_SIMILAR(aqm.getCorrelationCoefficient(), 0.44)
  TEST_EQUAL(aqm.getConcentrationUnits(), "uM")
  TEST_EQUAL(aqm.getTransformationModel(), "TransformationModelLinear")
  Param params2 = aqm.getTransformationModelParams();
  TEST_EQUAL(params2.getValue("slope"), 1)
}
END_SECTION

START_SECTION(bool checkLOD(const double value) const)
{
  AbsoluteQuantitationMethod aqm;
  const double value = 2.0;
  aqm.setLLOD(0.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value), true);
  aqm.setULOD(1.0);
  TEST_EQUAL(aqm.checkLOD(value), false);
  aqm.setLLOD(3.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value), false);
}
END_SECTION

START_SECTION(bool checkLOQ(const double value) const)
{
  AbsoluteQuantitationMethod aqm;
  const double value = 2.0;
  aqm.setLLOQ(0.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value), true);
  aqm.setULOQ(1.0);
  TEST_EQUAL(aqm.checkLOQ(value), false);
  aqm.setLLOQ(3.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value), false);
}
END_SECTION

START_SECTION(inline bool operator==(const AbsoluteQuantitationMethod& other) const)
{
  AbsoluteQuantitationMethod aqm1, aqm2;
  TEST_TRUE(aqm1 == aqm2);
  aqm1.setLLOD(1.0);
  aqm2.setLLOD(1.0);
  TEST_TRUE(aqm1 == aqm2);
  aqm2.setLLOD(2.0);
  TEST_EQUAL(aqm1 == aqm2, false);
}
END_SECTION

START_SECTION(inline bool operator!=(const AbsoluteQuantitationMethod& other) const)
{
  AbsoluteQuantitationMethod aqm1, aqm2;
  TEST_EQUAL(aqm1 != aqm2, false);
  aqm1.setLLOD(1.0);
  aqm2.setLLOD(1.0);
  TEST_EQUAL(aqm1 != aqm2, false);
  aqm2.setLLOD(2.0);
  TEST_FALSE(aqm1 == aqm2);
}
END_SECTION

START_SECTION(copyConstructor)
{
  AbsoluteQuantitationMethod aqm;

  aqm.setComponentName("component");
  aqm.setFeatureName("feature");
  aqm.setISName("IS");
  aqm.setLLOD(1.2);
  aqm.setULOD(3.4);
  aqm.setLLOQ(5.6);
  aqm.setULOQ(7.8);
  aqm.setNPoints(9);
  aqm.setCorrelationCoefficient(0.44);
  aqm.setConcentrationUnits("uM");
  aqm.setTransformationModel("TransformationModelLinear");
  Param params1;
  params1.setValue("slope", 1);
  aqm.setTransformationModelParams(params1);

  const AbsoluteQuantitationMethod aqm2(aqm);

  TEST_EQUAL(aqm2.getComponentName(), "component")
  TEST_EQUAL(aqm2.getFeatureName(), "feature")
  TEST_EQUAL(aqm2.getISName(), "IS")
  TEST_REAL_SIMILAR(aqm2.getLLOD(), 1.2)
  TEST_REAL_SIMILAR(aqm2.getULOD(), 3.4)
  TEST_REAL_SIMILAR(aqm2.getLLOQ(), 5.6)
  TEST_REAL_SIMILAR(aqm2.getULOQ(), 7.8)
  TEST_EQUAL(aqm2.getNPoints(), 9)
  TEST_REAL_SIMILAR(aqm2.getCorrelationCoefficient(), 0.44)
  TEST_EQUAL(aqm2.getConcentrationUnits(), "uM")
  TEST_EQUAL(aqm2.getTransformationModel(), "TransformationModelLinear")
  Param params2 = aqm2.getTransformationModelParams();
  TEST_EQUAL(params2.getValue("slope"), 1)
}
END_SECTION

START_SECTION(copyAssignmentOperator)
{
  AbsoluteQuantitationMethod aqm;

  aqm.setComponentName("component");
  aqm.setFeatureName("feature");
  aqm.setISName("IS");
  aqm.setLLOD(1.2);
  aqm.setULOD(3.4);
  aqm.setLLOQ(5.6);
  aqm.setULOQ(7.8);
  aqm.setNPoints(9);
  aqm.setCorrelationCoefficient(0.44);
  aqm.setConcentrationUnits("uM");
  aqm.setTransformationModel("TransformationModelLinear");
  Param params1;
  params1.setValue("slope", 1);
  aqm.setTransformationModelParams(params1);

  const AbsoluteQuantitationMethod aqm2 = aqm;

  TEST_EQUAL(aqm2.getComponentName(), "component")
  TEST_EQUAL(aqm2.getFeatureName(), "feature")
  TEST_EQUAL(aqm2.getISName(), "IS")
  TEST_REAL_SIMILAR(aqm2.getLLOD(), 1.2)
  TEST_REAL_SIMILAR(aqm2.getULOD(), 3.4)
  TEST_REAL_SIMILAR(aqm2.getLLOQ(), 5.6)
  TEST_REAL_SIMILAR(aqm2.getULOQ(), 7.8)
  TEST_EQUAL(aqm2.getNPoints(), 9)
  TEST_REAL_SIMILAR(aqm2.getCorrelationCoefficient(), 0.44)
  TEST_EQUAL(aqm2.getConcentrationUnits(), "uM")
  TEST_EQUAL(aqm2.getTransformationModel(), "TransformationModelLinear")
  Param params2 = aqm2.getTransformationModelParams();
  TEST_EQUAL(params2.getValue("slope"), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
