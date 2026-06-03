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
#include <OpenMS/IMAGING/MSImagingExperiment.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IMAGING/MSImagingRegion.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cmath>

using namespace OpenMS;
using namespace std;

namespace
{
  MSSpectrum makeSpec(std::initializer_list<std::pair<double, double>> peaks)
  {
    MSSpectrum s;
    for (const auto& p : peaks)
    {
      s.emplace_back(p.first, p.second);
    }
    return s;
  }

  // 2x2 grid; pixel (1,1) intentionally missing.
  MSImagingExperiment makeFixture()
  {
    MSExperiment exp;
    exp.addSpectrum(makeSpec({{499.95, 10.0}, {500.05, 20.0}}));
    exp.addSpectrum(makeSpec({{490.0, 5.0}}));
    exp.addSpectrum(makeSpec({{500.0, 100.0}}));

    MSImagingGeometry geom;
    geom.setDimensions(2, 2);
    geom.addPixel(0, 0, 0);
    geom.addPixel(1, 0, 1);
    geom.addPixel(0, 1, 2);

    MSImagingExperiment mie;
    mie.setMSExperiment(std::move(exp));
    mie.setGeometry(std::move(geom));
    return mie;
  }
}

START_TEST(MSImagingExperiment, "$Id$")

/////////////////////////////////////////////////////////////

MSImagingExperiment* ptr = nullptr;
MSImagingExperiment* null_ptr = nullptr;

START_SECTION((MSImagingExperiment()))
{
  ptr = new MSImagingExperiment();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getNumberOfPixels(), 0u)
  TEST_EQUAL(ptr->getNumberOfSpectra(), 0u)
}
END_SECTION

START_SECTION((~MSImagingExperiment()))
{
  delete ptr;
}
END_SECTION

START_SECTION((explicit MSImagingExperiment(MSExperiment exp)))
{
  MSExperiment exp;
  exp.addSpectrum(MSSpectrum());
  exp.addSpectrum(MSSpectrum());

  MSImagingExperiment mie(std::move(exp));
  TEST_EQUAL(mie.getNumberOfSpectra(), 2u)
  TEST_EQUAL(mie.getNumberOfPixels(), 0u)
}
END_SECTION

START_SECTION((MSImagingExperiment& operator=(MSExperiment exp)))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.getNumberOfSpectra(), 3u)
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)

  // assignment from a new MSExperiment resets the geometry
  MSExperiment fresh;
  fresh.addSpectrum(MSSpectrum());
  mie = std::move(fresh);
  TEST_EQUAL(mie.getNumberOfSpectra(), 1u)
  TEST_EQUAL(mie.getNumberOfPixels(), 0u)
  TEST_EQUAL(mie.getGeometry().getWidth(), 0u)
}
END_SECTION

START_SECTION((void setMSExperiment(MSExperiment exp)))
{
  MSExperiment exp;
  exp.addSpectrum(MSSpectrum());
  exp.addSpectrum(MSSpectrum());

  MSImagingExperiment mie;
  mie.setMSExperiment(std::move(exp));
  TEST_EQUAL(mie.getNumberOfSpectra(), 2u)
}
END_SECTION

START_SECTION((MSExperiment& getMSExperiment()))
{
  MSImagingExperiment mie;
  mie.getMSExperiment().addSpectrum(MSSpectrum());
  TEST_EQUAL(mie.getNumberOfSpectra(), 1u)
}
END_SECTION

START_SECTION((const MSExperiment& getMSExperiment() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setGeometry(MSImagingGeometry geom)))
{
  MSImagingGeometry g;
  g.setDimensions(3, 2);
  g.addPixel(0, 0, 0);
  MSImagingExperiment mie;
  mie.setGeometry(std::move(g));
  TEST_EQUAL(mie.getGeometry().getWidth(), 3u)
  TEST_EQUAL(mie.getNumberOfPixels(), 1u)
}
END_SECTION

START_SECTION((MSImagingGeometry& getGeometry()))
{
  MSImagingExperiment mie;
  mie.getGeometry().setDimensions(4, 4);
  TEST_EQUAL(mie.getGeometry().getWidth(), 4u)
}
END_SECTION

START_SECTION((const MSImagingGeometry& getGeometry() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size getNumberOfPixels() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size getNumberOfSpectra() const))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool hasPixel(UInt x, UInt y) const))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.hasPixel(0, 0), true)
  TEST_EQUAL(mie.hasPixel(1, 0), true)
  TEST_EQUAL(mie.hasPixel(0, 1), true)
  TEST_EQUAL(mie.hasPixel(1, 1), false)
}
END_SECTION

START_SECTION((MSSpectrum& getSpectrum(UInt x, UInt y)))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.getSpectrum(0, 1).size(), 1u)
  TEST_REAL_SIMILAR(mie.getSpectrum(0, 1)[0].getMZ(), 500.0)
  TEST_EXCEPTION(Exception::ElementNotFound, mie.getSpectrum(1, 1))

  // stale spectrum_index -> InvalidValue (not UB)
  mie.getGeometry().addPixel(1, 1, 999);
  TEST_EXCEPTION(Exception::InvalidValue, mie.getSpectrum(1, 1))
}
END_SECTION

START_SECTION((const MSSpectrum& getSpectrum(UInt x, UInt y) const))
{
  const MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.getSpectrum(0, 0).size(), 2u)
  TEST_EXCEPTION(Exception::ElementNotFound, mie.getSpectrum(1, 1))
}
END_SECTION

START_SECTION((IonImage extractIonImage(double mz, double tolerance_ppm) const))
{
  MSImagingExperiment mie = makeFixture();
  // 200 ppm @ 500 Da => dm = 0.1 => window [499.9, 500.1]
  IonImage img = mie.extractIonImage(500.0, 200.0);

  TEST_EQUAL(img.getWidth(), 2u)
  TEST_EQUAL(img.getHeight(), 2u)

  TEST_EQUAL(img.hasPixel(0, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 0), 30.0) // 10 + 20

  TEST_EQUAL(img.hasPixel(1, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(1, 0), 0.0)  // 490 outside window

  TEST_EQUAL(img.hasPixel(0, 1), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 1), 100.0)

  TEST_EQUAL(img.hasPixel(1, 1), false)

  TEST_REAL_SIMILAR(img.getMzRange().getMinMZ(), 499.9)
  TEST_REAL_SIMILAR(img.getMzRange().getMaxMZ(), 500.1)

  // invalid parameters
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(-1.0, 100.0))
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(500.0, -1.0))
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(std::nan(""), 100.0))

  // m/z window outside any peak -> all present pixels intensity 0
  IonImage far = mie.extractIonImage(1000.0, 100.0);
  TEST_EQUAL(far.hasPixel(0, 0), true)
  TEST_REAL_SIMILAR(far.getIntensity(0, 0), 0.0)
  TEST_EQUAL(far.hasPixel(0, 1), true)
  TEST_REAL_SIMILAR(far.getIntensity(0, 1), 0.0)

  // empty experiment -> empty image, no pixels visited
  MSImagingExperiment empty;
  IonImage empty_img = empty.extractIonImage(500.0, 100.0);
  TEST_EQUAL(empty_img.getWidth(), 0u)
  TEST_EQUAL(empty_img.getHeight(), 0u)
  TEST_REAL_SIMILAR(empty_img.getMzRange().getMinMZ(), 499.95)
  TEST_REAL_SIMILAR(empty_img.getMzRange().getMaxMZ(), 500.05)
}
END_SECTION

START_SECTION((IonImage extractIonImage(double mz, double tolerance_ppm, const MSImagingRegion& region) const))
{
  // Dedicated fixture with pixels at a non-zero offset to exercise local
  // coordinate translation (region bbox origin at (2,3)).
  MSExperiment exp;
  exp.addSpectrum(makeSpec({{500.0, 100.0}})); // 0 -> (2,3)
  exp.addSpectrum(makeSpec({{500.0, 50.0}}));  // 1 -> (3,3)
  exp.addSpectrum(makeSpec({{500.0, 25.0}}));  // 2 -> (2,4)
  exp.addSpectrum(makeSpec({{490.0, 5.0}}));   // 3 -> (4,4), outside any region below

  MSImagingGeometry geom;
  geom.setDimensions(6, 6);
  geom.addPixel(2, 3, 0);
  geom.addPixel(3, 3, 1);
  geom.addPixel(2, 4, 2);
  geom.addPixel(4, 4, 3);

  MSImagingExperiment mie;
  mie.setMSExperiment(std::move(exp));
  mie.setGeometry(std::move(geom));

  // Rectangle region [2,3] x [3,4]: bbox 2x2, origin (2,3).
  MSImagingRegion rect = MSImagingRegion::rectangle(1, "rect", 2, 3, 3, 4);
  IonImage img = mie.extractIonImage(500.0, 200.0, rect); // window [499.9, 500.1]

  // image is cropped to the region's bounding box
  TEST_EQUAL(img.getWidth(), 2u)
  TEST_EQUAL(img.getHeight(), 2u)

  // global (2,3) -> local (0,0)
  TEST_EQUAL(img.hasPixel(0, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 0), 100.0)
  // global (3,3) -> local (1,0)
  TEST_EQUAL(img.hasPixel(1, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(1, 0), 50.0)
  // global (2,4) -> local (0,1)
  TEST_EQUAL(img.hasPixel(0, 1), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 1), 25.0)
  // global (3,4) -> local (1,1): no pixel acquired there -> invalid
  TEST_EQUAL(img.hasPixel(1, 1), false)

  TEST_REAL_SIMILAR(img.getMzRange().getMinMZ(), 499.9)
  TEST_REAL_SIMILAR(img.getMzRange().getMaxMZ(), 500.1)

  // Mask region over the same bbox: drop (3,3) by leaving its cell unset.
  //   row 0 (y=3): (2,3)=1 (3,3)=0
  //   row 1 (y=4): (2,4)=1 (3,4)=1
  std::vector<bool> mask = {true, false,
                            true, true};
  MSImagingRegion mreg = MSImagingRegion::fromMask(2, "mask", 2, 3, 2, 2, mask);
  IonImage mimg = mie.extractIonImage(500.0, 200.0, mreg);
  TEST_EQUAL(mimg.getWidth(), 2u)
  TEST_EQUAL(mimg.getHeight(), 2u)
  TEST_EQUAL(mimg.hasPixel(0, 0), true)              // (2,3) in mask
  TEST_REAL_SIMILAR(mimg.getIntensity(0, 0), 100.0)
  TEST_EQUAL(mimg.hasPixel(1, 0), false)             // (3,3) excluded by mask
  TEST_EQUAL(mimg.hasPixel(0, 1), true)              // (2,4) in mask
  TEST_REAL_SIMILAR(mimg.getIntensity(0, 1), 25.0)
  TEST_EQUAL(mimg.hasPixel(1, 1), false)             // (3,4) in mask but no pixel acquired

  // region pixel with a spectrum but no peaks in the window -> valid, intensity 0
  MSImagingRegion far_rect = MSImagingRegion::rectangle(3, "far", 4, 4, 4, 4);
  IonImage far = mie.extractIonImage(500.0, 200.0, far_rect);
  TEST_EQUAL(far.getWidth(), 1u)
  TEST_EQUAL(far.getHeight(), 1u)
  TEST_EQUAL(far.hasPixel(0, 0), true)               // (4,4) present
  TEST_REAL_SIMILAR(far.getIntensity(0, 0), 0.0)     // 490 outside window

  // region overlapping no acquired pixels -> sized to bbox, all invalid
  MSImagingRegion none = MSImagingRegion::rectangle(4, "none", 0, 0, 1, 1);
  IonImage none_img = mie.extractIonImage(500.0, 200.0, none);
  TEST_EQUAL(none_img.getWidth(), 2u)
  TEST_EQUAL(none_img.getHeight(), 2u)
  TEST_EQUAL(none_img.hasPixel(0, 0), false)

  // invalid parameters rejected
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(-1.0, 100.0, rect))
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(500.0, -1.0, rect))
  TEST_EXCEPTION(Exception::InvalidValue, mie.extractIonImage(std::nan(""), 100.0, rect))
}
END_SECTION

START_SECTION((void validate() const))
{
  // consistent fixture passes
  MSImagingExperiment mie = makeFixture();
  mie.validate();
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)

  // append a 4th pixel that maps to a non-existent spectrum -> throws
  mie.getGeometry().addPixel(1, 1, 999);
  TEST_EXCEPTION(Exception::InvalidValue, mie.validate())
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
