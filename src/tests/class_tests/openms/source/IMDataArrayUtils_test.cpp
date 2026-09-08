// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS;

START_TEST(IMDataArrayUtils, "$Id$")

START_SECTION((unit round trips and compatibility entry points))
  for (const auto expected : {DriftTimeUnit::MILLISECOND, DriftTimeUnit::VSSC})
  {
    DataArrays::FloatDataArray direct, legacy;
    direct.push_back(1.25f);
    IMDataArrayUtils::setIMUnit(direct, expected);
    IMDataConverter::setIMUnit(legacy, expected);
    TEST_EQUAL(direct.getName(), legacy.getName())
    TEST_REAL_SIMILAR(direct[0], 1.25f)
    DriftTimeUnit actual = DriftTimeUnit::NONE;
    TEST_TRUE(IMDataArrayUtils::getIMUnit(direct, actual))
    TEST_TRUE(actual == expected)
    actual = DriftTimeUnit::NONE;
    TEST_TRUE(IMDataConverter::getIMUnit(direct, actual))
    TEST_TRUE(actual == expected)
  }
  DataArrays::FloatDataArray array;
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::NONE))
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE))
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::CCS))
END_SECTION

START_SECTION((ontology and vendor names retain their units))
  DataArrays::FloatDataArray array;
  DriftTimeUnit unit = DriftTimeUnit::NONE;
  array.setName("mean inverse reduced ion mobility array");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(Constants::UserParam::ION_MOBILITY);
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::MILLISECOND)
  array.setName(Constants::UserParam::INVERSE_REDUCED_ION_MOBILITY);
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(std::string(Constants::UserParam::ION_MOBILITY) + " MS:1002954");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::CCS)
  array.setName("unrelated auxiliary data");
  TEST_FALSE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::CCS)
END_SECTION

END_TEST
