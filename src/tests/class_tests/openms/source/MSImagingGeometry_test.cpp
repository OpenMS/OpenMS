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
#include <OpenMS/IMAGING/MSImagingGeometry.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS;
using namespace std;

START_TEST(MSImagingGeometry, "$Id$")

/////////////////////////////////////////////////////////////

MSImagingGeometry* ptr = nullptr;
MSImagingGeometry* null_ptr = nullptr;

START_SECTION((MSImagingGeometry()))
{
  ptr = new MSImagingGeometry();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getWidth(), 0u)
  TEST_EQUAL(ptr->getHeight(), 0u)
  TEST_EQUAL(ptr->getNumberOfPixels(), 0u)
  TEST_REAL_SIMILAR(ptr->getPixelSizeX(), 1.0)
  TEST_REAL_SIMILAR(ptr->getPixelSizeY(), 1.0)
  TEST_EQUAL(ptr->getPixelSizeUnit(), "micrometer")
}
END_SECTION

START_SECTION((~MSImagingGeometry()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void setDimensions(UInt width, UInt height)))
{
  MSImagingGeometry g;
  g.setDimensions(10, 8);
  TEST_EQUAL(g.getWidth(), 10u)
  TEST_EQUAL(g.getHeight(), 8u)
}
END_SECTION

START_SECTION((void setPixelSize(double x, double y, const String& unit)))
{
  MSImagingGeometry g;
  g.setPixelSize(25.0, 25.0, "micrometer");
  TEST_REAL_SIMILAR(g.getPixelSizeX(), 25.0)
  TEST_REAL_SIMILAR(g.getPixelSizeY(), 25.0)
  TEST_EQUAL(g.getPixelSizeUnit(), "micrometer")
  g.setPixelSize(50.0, 75.0, "nanometer");
  TEST_REAL_SIMILAR(g.getPixelSizeX(), 50.0)
  TEST_REAL_SIMILAR(g.getPixelSizeY(), 75.0)
  TEST_EQUAL(g.getPixelSizeUnit(), "nanometer")
}
END_SECTION

START_SECTION((void addPixel(UInt x, UInt y, Size spectrum_index)))
{
  MSImagingGeometry g;
  g.addPixel(0, 0, 42);
  g.addPixel(1, 0, 7);
  TEST_EQUAL(g.getNumberOfPixels(), 2u)
  TEST_EQUAL(g.hasPixel(0, 0), true)
  TEST_EQUAL(g.hasPixel(1, 0), true)
  TEST_EQUAL(g.hasPixel(2, 0), false)
  TEST_EQUAL(g.getSpectrumIndex(0, 0), 42u)
  TEST_EQUAL(g.getSpectrumIndex(1, 0), 7u)

  // duplicate insertion
  TEST_EXCEPTION(Exception::InvalidValue, g.addPixel(0, 0, 99))

  // when dimensions are set, out-of-grid pixels are rejected
  MSImagingGeometry g2;
  g2.setDimensions(3, 2);
  g2.addPixel(2, 1, 5); // in-bounds corner
  TEST_EXCEPTION(Exception::InvalidValue, g2.addPixel(3, 0, 6))
  TEST_EXCEPTION(Exception::InvalidValue, g2.addPixel(0, 2, 7))
}
END_SECTION

START_SECTION((bool hasPixel(UInt x, UInt y) const))
{
  MSImagingGeometry g;
  g.addPixel(5, 9, 100);
  TEST_EQUAL(g.hasPixel(5, 9), true)
  TEST_EQUAL(g.hasPixel(0, 0), false)
}
END_SECTION

START_SECTION((Size getSpectrumIndex(UInt x, UInt y) const))
{
  MSImagingGeometry g;
  g.addPixel(2, 3, 17);
  TEST_EQUAL(g.getSpectrumIndex(2, 3), 17u)
  TEST_EXCEPTION(Exception::ElementNotFound, g.getSpectrumIndex(0, 0))
}
END_SECTION

START_SECTION((const std::vector<Pixel>& getPixels() const))
{
  MSImagingGeometry g;
  g.addPixel(2, 0, 10);
  g.addPixel(0, 1, 20);
  g.addPixel(1, 1, 30);
  const auto& pixels = g.getPixels();
  TEST_EQUAL(pixels.size(), 3u)
  TEST_EQUAL(pixels[0].x, 2u); TEST_EQUAL(pixels[0].y, 0u); TEST_EQUAL(pixels[0].spectrum_index, 10u)
  TEST_EQUAL(pixels[1].x, 0u); TEST_EQUAL(pixels[1].y, 1u); TEST_EQUAL(pixels[1].spectrum_index, 20u)
  TEST_EQUAL(pixels[2].x, 1u); TEST_EQUAL(pixels[2].y, 1u); TEST_EQUAL(pixels[2].spectrum_index, 30u)
}
END_SECTION

START_SECTION((Size getNumberOfPixels() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void clear()))
{
  MSImagingGeometry g;
  g.setDimensions(4, 4);
  g.setPixelSize(10, 10);
  g.addPixel(0, 0, 1);
  g.addPixel(1, 1, 2);
  g.clear();
  TEST_EQUAL(g.getWidth(), 0u)
  TEST_EQUAL(g.getHeight(), 0u)
  TEST_EQUAL(g.getNumberOfPixels(), 0u)
  TEST_EQUAL(g.hasPixel(0, 0), false)
  TEST_REAL_SIMILAR(g.getPixelSizeX(), 1.0)
  TEST_EQUAL(g.getPixelSizeUnit(), "micrometer")
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

START_SECTION((double getPixelSizeX() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((double getPixelSizeY() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((const String& getPixelSizeUnit() const))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
