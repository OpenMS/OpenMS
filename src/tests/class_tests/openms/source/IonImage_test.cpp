// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/IMAGING/IonImage.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/RangeManager.h>

using namespace OpenMS;
using namespace std;

START_TEST(IonImage, "$Id$")

/////////////////////////////////////////////////////////////

IonImage* ptr = nullptr;
IonImage* null_ptr = nullptr;

START_SECTION((IonImage()))
{
  ptr = new IonImage();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getWidth(), 0u)
  TEST_EQUAL(ptr->getHeight(), 0u)
  TEST_EQUAL(ptr->hasPixel(0, 0), false)
}
END_SECTION

START_SECTION((~IonImage()))
{
  delete ptr;
}
END_SECTION

START_SECTION((IonImage(UInt width, UInt height)))
{
  IonImage img(3, 2);
  TEST_EQUAL(img.getWidth(), 3u)
  TEST_EQUAL(img.getHeight(), 2u)
  TEST_EQUAL(img.getData().size(), 6u)
  TEST_EQUAL(img.getMask().size(), 6u)
  for (UInt y = 0; y < 2; ++y)
  {
    for (UInt x = 0; x < 3; ++x)
    {
      TEST_EQUAL(img.hasPixel(x, y), false)
      TEST_REAL_SIMILAR(img.getIntensity(x, y), 0.0)
    }
  }
}
END_SECTION

START_SECTION((void resize(UInt width, UInt height)))
{
  IonImage img;
  img.resize(2, 2);
  img.setIntensity(0, 0, 1.0);
  img.resize(4, 3);
  TEST_EQUAL(img.getWidth(), 4u)
  TEST_EQUAL(img.getHeight(), 3u)
  TEST_EQUAL(img.getData().size(), 12u)
  TEST_EQUAL(img.hasPixel(0, 0), false)
  TEST_REAL_SIMILAR(img.getIntensity(0, 0), 0.0)
}
END_SECTION

START_SECTION((void setIntensity(UInt x, UInt y, double intensity)))
{
  IonImage img(3, 2);
  img.setIntensity(1, 0, 7.5);
  TEST_EQUAL(img.hasPixel(1, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(1, 0), 7.5)
  TEST_EQUAL(img.hasPixel(2, 1), false)
  TEST_REAL_SIMILAR(img.getIntensity(2, 1), 0.0)
  // index = y * W + x = 1 * 3 + 2 = 5
  img.setIntensity(2, 1, 42.0);
  TEST_REAL_SIMILAR(img.getData()[5], 42.0)
  TEST_TRUE(img.getMask()[5])
}
END_SECTION

START_SECTION((void setMzRange(const RangeMZ& range)))
{
  IonImage img;
  img.setMzRange(RangeMZ(499.9, 500.1));
  TEST_REAL_SIMILAR(img.getMzRange().getMinMZ(), 499.9)
  TEST_REAL_SIMILAR(img.getMzRange().getMaxMZ(), 500.1)
}
END_SECTION

START_SECTION((const RangeMZ& getMzRange() const))
{
  NOT_TESTABLE // covered by setMzRange
}
END_SECTION

START_SECTION((double getIntensity(UInt x, UInt y) const))
{
  IonImage img(3, 2);
  TEST_EXCEPTION(Exception::IndexOverflow, img.getIntensity(5, 0))
  TEST_EXCEPTION(Exception::IndexOverflow, img.getIntensity(0, 5))
}
END_SECTION

START_SECTION((bool hasPixel(UInt x, UInt y) const))
{
  IonImage img(3, 2);
  TEST_EQUAL(img.hasPixel(5, 0), false) // OOB returns false without throwing
  TEST_EQUAL(img.hasPixel(0, 5), false)
}
END_SECTION

START_SECTION((UInt getWidth() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((UInt getHeight() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((const std::vector<double>& getData() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((const std::vector<bool>& getMask() const))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
