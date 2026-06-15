// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann $
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
{ delete ptr; }
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
  TEST_EQUAL(pixels[0].x, 2u);
  TEST_EQUAL(pixels[0].y, 0u);
  TEST_EQUAL(pixels[0].spectrum_index, 10u)
  TEST_EQUAL(pixels[1].x, 0u);
  TEST_EQUAL(pixels[1].y, 1u);
  TEST_EQUAL(pixels[1].spectrum_index, 20u)
  TEST_EQUAL(pixels[2].x, 1u);
  TEST_EQUAL(pixels[2].y, 1u);
  TEST_EQUAL(pixels[2].spectrum_index, 30u)
}
END_SECTION

START_SECTION((Size getNumberOfPixels() const)) {NOT_TESTABLE} END_SECTION

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

START_SECTION((UInt getWidth() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((UInt getHeight() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((double getPixelSizeX() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((double getPixelSizeY() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((const String& getPixelSizeUnit() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((addRegion()))
{
  MSImagingGeometry g;
  g.setDimensions(6, 6);
  g.addPixel(0, 0, 1);
  g.addPixel(1, 0, 2);
  g.addPixel(5, 5, 3); // outside of region
  g.addRegion(MSImagingRegion::rectangle(1, "name", 0, 0, 1, 0));
  TEST_EQUAL(g.getNumberOfRegions(), 1u);

  auto spec_idx = g.getRegionSpectrumIndices(1);
  TEST_EQUAL(spec_idx.size(), 2u);
  TEST_EQUAL(spec_idx[0], 1u);
  TEST_EQUAL(spec_idx[1], 2u);

  TEST_EXCEPTION(Exception::InvalidValue, g.addRegion(MSImagingRegion::rectangle(2, "collision", 0, 0, 10, 10)));
  TEST_EQUAL(g.getNumberOfRegions(), 1u); // ensure region was not added

  TEST_EXCEPTION(Exception::InvalidValue, g.addRegion(MSImagingRegion::rectangle(1, "dup", 5, 5, 6, 6))); // duplicate id
  TEST_EXCEPTION(Exception::InvalidValue,
                 g.addRegion(MSImagingRegion::rectangle(MSImagingGeometry::NO_REGION, "bad id", 5, 5, 6, 6))); // do not set sentinel id

  g.addRegion(MSImagingRegion::rectangle(2, "far away", 100, 100, 100, 100));
  TEST_EQUAL(g.getNumberOfRegions(), 2u);
  TEST_EQUAL(g.getRegionPixels(2).size(), 0u);
}
END_SECTION

START_SECTION((regionOf()))
{
  MSImagingGeometry g;
  g.addPixel(0, 0, 0);
  g.addPixel(5, 5, 1);
  g.addRegion(MSImagingRegion::rectangle(1, "name", 0, 0, 1, 0));

  TEST_EQUAL(g.regionOf(0, 0), 1u);
  TEST_EQUAL(g.regionOf(5, 5), MSImagingGeometry::NO_REGION);   // exists but not owned
  TEST_EQUAL(g.regionOf(10, 10), MSImagingGeometry::NO_REGION); // never acquired
}
END_SECTION

START_SECTION((getRegionPixels()))
{
  MSImagingGeometry g;

  g.addPixel(0, 0, 1);
  g.addPixel(1, 0, 2);
  g.addRegion(MSImagingRegion::rectangle(1, "name", 0, 0, 1, 0));

  auto pixel_idx = g.getRegionPixels(1);
  TEST_EQUAL(pixel_idx[0], 0u); // pixel (0,0) is at position 0 in getPixels()
  TEST_EQUAL(pixel_idx[1], 1u); // pixel (1,0) is at position 1
}
END_SECTION

START_SECTION((getRegionSpectrumIndices()))
{
  MSImagingGeometry g;

  g.addPixel(0, 0, 1);
  g.addPixel(1, 0, 2);
  g.addRegion(MSImagingRegion::rectangle(1, "name", 0, 0, 1, 0));

  auto spec_idx = g.getRegionSpectrumIndices(1);
  TEST_EQUAL(spec_idx.size(), 2u);
  TEST_EQUAL(spec_idx[0], 1u);
  TEST_EQUAL(spec_idx[1], 2u);
}
END_SECTION


START_SECTION((region overlay : Mask region))
{
  MSImagingGeometry g;
  g.addPixel(10, 20, 1); // mask bit 0 (true)
  g.addPixel(11, 20, 2); // mask bit 1 (false) exists but unowned
  g.addPixel(11, 21, 3); // mask bit 3 (true)
  // 2x2 mask, origin (10,20): {true,false,false,true}
  g.addRegion(MSImagingRegion::fromMask(1, "m", 10, 20, 2, 2, std::vector<bool> {true, false, false, true}));
  TEST_EQUAL(g.regionOf(10, 20), 1u)                           // set bit -> owned
  TEST_EQUAL(g.regionOf(11, 20), MSImagingGeometry::NO_REGION) // the unowned pixel
  TEST_EQUAL(g.regionOf(11, 21), 1u)
  TEST_EQUAL(g.getRegionSpectrumIndices(1).size(), 2u) // region is 2x2 but only 2 spectra are set
}
END_SECTION

START_SECTION((removeRegion()))
{
  MSImagingGeometry g;
  g.addPixel(0, 0, 1);
  g.addPixel(5, 5, 2);
  g.addRegion(MSImagingRegion::rectangle(1, "A", 0, 0, 0, 0));
  g.addRegion(MSImagingRegion::rectangle(2, "B", 5, 5, 5, 5));
  g.removeRegion(1);
  TEST_EQUAL(g.getNumberOfRegions(), 1u);
  TEST_EQUAL(g.regionOf(0, 0), MSImagingGeometry::NO_REGION);
  TEST_EQUAL(g.regionOf(5, 5), 2u);
  TEST_EQUAL(g.getNumberOfPixels(), 2u); // removed region but NOT pixel
}
END_SECTION

START_SECTION((clearRegions()))
{
  MSImagingGeometry g;
  g.addPixel(0, 0, 1);
  g.addRegion(MSImagingRegion::rectangle(1, "A", 0, 0, 0, 0));
  g.clearRegions();
  TEST_EQUAL(g.getNumberOfRegions(), 0u);
  TEST_EQUAL(g.regionOf(0, 0), MSImagingGeometry::NO_REGION);
  TEST_EQUAL(g.getNumberOfPixels(), 1u);
}
END_SECTION

START_SECTION((getRegion()))
{
  MSImagingGeometry g;
  g.addRegion(MSImagingRegion::rectangle(1, "A", 0, 0, 1, 1));
  TEST_EQUAL(g.getRegion(1).getId(), 1u)
  TEST_EXCEPTION(Exception::ElementNotFound, g.getRegion(99))
  TEST_EXCEPTION(Exception::ElementNotFound, g.getRegionPixels(99))
  TEST_EXCEPTION(Exception::ElementNotFound, g.removeRegion(99))
}
END_SECTION

START_SECTION((membership staleness test))
{
  // we add a pixel after the region to ensure that state is live and not dependent on order of operations
  MSImagingGeometry g;
  g.addRegion(MSImagingRegion::rectangle(1, "A", 0, 0, 0, 0));
  g.addPixel(0, 0, 1);
  TEST_EQUAL(g.regionOf(0, 0), 1u);
  TEST_EQUAL(g.getRegionSpectrumIndices(1).size(), 1u);
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST