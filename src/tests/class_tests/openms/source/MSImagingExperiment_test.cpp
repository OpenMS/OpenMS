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
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/IMAGING/MSImagingRegion.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <cmath>
#include <iterator>
#include <utility>

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

MSExperiment makeSpecExperiment()
{
  MSExperiment exp;
  exp.addSpectrum(makeSpec({{499.95, 10.0}, {500.05, 20.0}}));
  exp.addSpectrum(makeSpec({{490.0, 5.0}}));
  exp.addSpectrum(makeSpec({{500.0, 100.0}}));
  return exp;
}

// 2x2 grid; pixel (1,1) intentionally missing.
MSImagingExperiment makeFixture()
{
  MSExperiment exp = makeSpecExperiment();

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
} // namespace

START_TEST(MSImagingExperiment, "$Id$")

/////////////////////////////////////////////////////////////

MSImagingExperiment* ptr = nullptr;
MSImagingExperiment* null_ptr = nullptr;

START_SECTION((MSImagingExperiment()))
{
  ptr = new MSImagingExperiment();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getNumberOfPixels(), 0u)
  TEST_EQUAL(ptr->getNrSpectra(), 0u)
}
END_SECTION

START_SECTION((~MSImagingExperiment()))
{ delete ptr; }
END_SECTION

START_SECTION((explicit MSImagingExperiment(MSExperiment exp)))
{
  MSExperiment exp;
  exp.addSpectrum(MSSpectrum());
  exp.addSpectrum(MSSpectrum());

  MSImagingExperiment mie(std::move(exp));
  TEST_EQUAL(mie.getNrSpectra(), 2u)
  TEST_EQUAL(mie.getNumberOfPixels(), 0u)
}
END_SECTION

START_SECTION((MSImagingExperiment & operator=(MSExperiment exp)))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.getNrSpectra(), 3u)
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)

  // assignment from a new MSExperiment resets the geometry
  MSExperiment fresh;
  fresh.addSpectrum(MSSpectrum());
  mie = std::move(fresh);
  TEST_EQUAL(mie.getNrSpectra(), 1u)
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
  TEST_EQUAL(mie.getNrSpectra(), 2u)
}
END_SECTION

START_SECTION((MSExperiment & getMSExperiment()))
{
  MSImagingExperiment mie;
  mie.getMSExperiment().addSpectrum(MSSpectrum());
  TEST_EQUAL(mie.getNrSpectra(), 1u)
}
END_SECTION

START_SECTION((const MSExperiment& getMSExperiment() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((void setGeometry(MSImagingGeometry geom)))
{
  MSImagingGeometry g;
  g.setDimensions(3, 2);
  g.addPixel(0, 0, 0);
  MSImagingExperiment mie;
  mie.setGeometry(std::move(g)); // no spectra yet: nothing to annotate, pixel is kept
  TEST_EQUAL(mie.getGeometry().getWidth(), 3u)
  TEST_EQUAL(mie.getNumberOfPixels(), 1u)

  // with spectra present, the bound spectra get their (1-based) pixel coordinate recorded
  MSImagingExperiment fx = makeFixture();
  for (Size i = 0; i < 3; ++i)
  {
    TEST_EQUAL(fx[i].metaValueExists(MSImagingExperiment::META_PIXEL_X), true)
    TEST_EQUAL(fx[i].metaValueExists(MSImagingExperiment::META_PIXEL_Y), true)
    TEST_EQUAL((Int)fx[i].getMetaValue(MSImagingExperiment::META_PIXEL_Z), 1)
  }
  TEST_EQUAL((Int)fx[1].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2) // pixel (1,0) -> spectrum 1
  TEST_EQUAL((Int)fx[1].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 1)
  TEST_EQUAL((Int)fx[2].getMetaValue(MSImagingExperiment::META_PIXEL_X), 1) // pixel (0,1) -> spectrum 2
  TEST_EQUAL((Int)fx[2].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 2)

  // the geometry passed in is authoritative: a spectrum's previous coordinate is overwritten
  MSImagingGeometry moved;
  moved.setDimensions(2, 2);
  moved.addPixel(1, 1, 1);
  fx.setGeometry(std::move(moved));
  TEST_EQUAL((Int)fx[1].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2)
  TEST_EQUAL((Int)fx[1].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 2)
  fx.validate();

  // a geometry binding one spectrum to two pixels is rejected before anything is touched
  MSImagingGeometry twice;
  twice.setDimensions(2, 2);
  twice.addPixel(0, 0, 1);
  twice.addPixel(0, 1, 1);
  TEST_EXCEPTION(Exception::InvalidValue, fx.setGeometry(std::move(twice)))
  TEST_EQUAL(fx.getNumberOfPixels(), 1u)
  TEST_EQUAL(fx.hasPixel(1, 1), true)
  TEST_EQUAL((Int)fx[1].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2)
  fx.validate();
}
END_SECTION

START_SECTION((MSImagingGeometry & getGeometry()))
{
  MSImagingExperiment mie;
  mie.getGeometry().setDimensions(4, 4);
  TEST_EQUAL(mie.getGeometry().getWidth(), 4u)
}
END_SECTION

START_SECTION((const MSImagingGeometry& getGeometry() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((Size getNumberOfPixels() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((Size getNrSpectra() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((Size size() const))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.size(), 3u)
  TEST_EQUAL(MSImagingExperiment().size(), 0u)
}
END_SECTION

START_SECTION((bool empty() const))
{
  TEST_EQUAL(MSImagingExperiment().empty(), true)
  TEST_EQUAL(makeFixture().empty(), false)
}
END_SECTION

START_SECTION((MSSpectrum& getSpectrum(Size index)))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.getSpectrum(Size(0)).size(), 2u)
  mie.getSpectrum(Size(1)).setRT(42.0); // a reference: edits land
  TEST_REAL_SIMILAR(mie.getMSExperiment()[1].getRT(), 42.0)
  TEST_EXCEPTION(Exception::IndexOverflow, mie.getSpectrum(Size(3)))
}
END_SECTION

START_SECTION((const MSSpectrum& getSpectrum(Size index) const))
{
  const MSImagingExperiment mie = makeFixture();
  TEST_REAL_SIMILAR(mie.getSpectrum(Size(2))[0].getMZ(), 500.0)
  TEST_EXCEPTION(Exception::IndexOverflow, mie.getSpectrum(Size(99)))
}
END_SECTION

START_SECTION((MSSpectrum& operator[](Size index)))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie[0].size(), 2u)
  TEST_EXCEPTION(Exception::IndexOverflow, mie[3])
}
END_SECTION

START_SECTION((const MSSpectrum& operator[](Size index) const))
{
  const MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie[1].size(), 1u)
  TEST_EXCEPTION(Exception::IndexOverflow, mie[3])
}
END_SECTION

START_SECTION((Iterator begin()))
{
  MSImagingExperiment mie = makeFixture();
  Size n = 0;
  for (MSSpectrum& s : mie)
  {
    s.setRT(double(n));
    ++n;
  }
  TEST_EQUAL(n, 3u)
  TEST_REAL_SIMILAR(mie[2].getRT(), 2.0)
}
END_SECTION

START_SECTION((Iterator end())) {NOT_TESTABLE} END_SECTION

START_SECTION((ConstIterator begin() const))
{
  const MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(std::distance(mie.begin(), mie.end()), 3)
}
END_SECTION

START_SECTION((ConstIterator end() const)) {NOT_TESTABLE} END_SECTION

  START_SECTION((bool hasPixel(UInt x, UInt y) const))
{
  MSImagingExperiment mie = makeFixture();
  TEST_EQUAL(mie.hasPixel(0, 0), true)
  TEST_EQUAL(mie.hasPixel(1, 0), true)
  TEST_EQUAL(mie.hasPixel(0, 1), true)
  TEST_EQUAL(mie.hasPixel(1, 1), false)
}
END_SECTION

START_SECTION((void bindPixel(UInt x, UInt y, Size spectrum_index)))
{
  MSImagingExperiment mie(makeSpecExperiment());
  mie.getGeometry().setDimensions(2, 2);
  mie.bindPixel(1, 1, 2);
  TEST_EQUAL(mie.hasPixel(1, 1), true)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(1, 1), 2u)
  TEST_EQUAL((Int)mie[2].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2) // recorded 1-based
  TEST_EQUAL((Int)mie[2].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 2)
  TEST_EQUAL((Int)mie[2].getMetaValue(MSImagingExperiment::META_PIXEL_Z), 1)
  mie.validate();

  TEST_EXCEPTION(Exception::IndexOverflow, mie.bindPixel(0, 0, 3)) // no such spectrum
  TEST_EXCEPTION(Exception::InvalidValue, mie.bindPixel(1, 1, 0))  // duplicate pixel
  TEST_EXCEPTION(Exception::InvalidValue, mie.bindPixel(2, 0, 0))  // outside the 2x2 grid
  TEST_EQUAL(mie[0].metaValueExists(MSImagingExperiment::META_PIXEL_X), false) // nothing written on failure

  // a spectrum has one acquisition location: a second binding of spectrum 2 is rejected
  // and leaves the single existing binding intact
  TEST_EXCEPTION(Exception::InvalidValue, mie.bindPixel(0, 0, 2))
  TEST_EQUAL(mie.getNumberOfPixels(), 1u)
  TEST_EQUAL(mie.hasPixel(0, 0), false)
  TEST_EQUAL((Int)mie[2].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2)
  TEST_EQUAL((Int)mie[2].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 2)
  mie.validate();
  // rebinding to the very same pixel is a duplicate-pixel error, as before
  TEST_EXCEPTION(Exception::InvalidValue, mie.bindPixel(1, 1, 2))
  // a spectrum that merely carries a coordinate (no pixel claims it) can be bound anywhere
  MSImagingExperiment::setPixelCoordinate(mie[0], 1, 0);
  mie.bindPixel(0, 0, 0);
  TEST_EQUAL((Int)mie[0].getMetaValue(MSImagingExperiment::META_PIXEL_X), 1)
  TEST_EQUAL((Int)mie[0].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 1)
  mie.validate();
}
END_SECTION

START_SECTION((void unbindPixel(UInt x, UInt y)))
{
  MSImagingExperiment mie = makeFixture();
  mie.getGeometry().addRegion(MSImagingRegion::rectangle(1, "col0", 0, 0, 0, 1));
  mie.unbindPixel(1, 0); // spectrum 1
  TEST_EQUAL(mie.getNumberOfPixels(), 2u)
  TEST_EQUAL(mie.hasPixel(1, 0), false)
  TEST_EQUAL(mie[1].metaValueExists(MSImagingExperiment::META_PIXEL_X), false)
  TEST_EQUAL(mie[1].metaValueExists(MSImagingExperiment::META_PIXEL_Y), false)
  TEST_EQUAL(mie[1].metaValueExists(MSImagingExperiment::META_PIXEL_Z), false)
  TEST_EQUAL(mie.getNrSpectra(), 3u) // the spectrum itself stays
  TEST_EQUAL(mie.getGeometry().getNumberOfRegions(), 1u)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(0, 1), 2u) // remaining lookups still resolve
  mie.validate();
  TEST_EXCEPTION(Exception::ElementNotFound, mie.unbindPixel(1, 0))

  // ... which is how a binding is moved
  mie.bindPixel(1, 1, 1);
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(1, 1), 1u)
  TEST_EQUAL((Int)mie[1].getMetaValue(MSImagingExperiment::META_PIXEL_X), 2)
  TEST_EQUAL((Int)mie[1].getMetaValue(MSImagingExperiment::META_PIXEL_Y), 2)
  mie.validate();
  mie.rebuildGeometry();
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(1, 1), 1u)
}
END_SECTION

START_SECTION((MSSpectrum & getSpectrum(UInt x, UInt y)))
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
  TEST_REAL_SIMILAR(img.getIntensity(1, 0), 0.0) // 490 outside window

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

START_SECTION((void validate() const))
{
  // consistent fixture passes
  MSImagingExperiment mie = makeFixture();
  mie.validate();
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)

  // append a 4th pixel that maps to a non-existent spectrum -> throws
  mie.getGeometry().addPixel(1, 1, 999);
  TEST_EXCEPTION(Exception::InvalidValue, mie.validate())

  // every index in range, but the spectrum list was reordered: the coordinate each spectrum
  // carries no longer matches the pixel it is bound to -> throws (an index-only check passes)
  MSImagingExperiment reordered = makeFixture();
  std::swap(reordered.getMSExperiment().getSpectra()[0], reordered.getMSExperiment().getSpectra()[2]);
  TEST_EXCEPTION(Exception::InvalidValue, reordered.validate())

  // spectra without coordinates are accepted as-is (hand-built geometry over plain spectra)
  MSImagingExperiment plain(makeSpecExperiment());
  plain.getGeometry().setDimensions(2, 2);
  plain.getGeometry().addPixel(0, 0, 0);
  plain.validate();

  // a stale index that happens to be in range on a replaced experiment is caught too
  MSImagingExperiment replaced = makeFixture();
  MSExperiment other = replaced.getMSExperiment();
  other.getSpectra().erase(other.getSpectra().begin()); // drop spectrum 0: every index shifts by one
  replaced.setMSExperiment(std::move(other));
  TEST_EXCEPTION(Exception::InvalidValue, replaced.validate())
}
END_SECTION

START_SECTION((void rebuildGeometry()))
{
  MSImagingExperiment mie = makeFixture();
  mie.getGeometry().setPixelSize(20.0, 20.0, "micrometer");
  mie.getGeometry().addRegion(MSImagingRegion::rectangle(1, "col0", 0, 0, 0, 1));

  // reorder the spectra behind the geometry's back, then rebuild the index from the
  // coordinates the spectra carry: dims, pixel size and regions survive, indices follow
  MSExperiment& exp = mie.getMSExperiment();
  std::swap(exp.getSpectra()[0], exp.getSpectra()[2]);
  TEST_EXCEPTION(Exception::InvalidValue, mie.validate())
  mie.rebuildGeometry();
  mie.validate();
  TEST_EQUAL(mie.getNumberOfPixels(), 3u)
  TEST_EQUAL(mie.getGeometry().getWidth(), 2u)
  TEST_EQUAL(mie.getGeometry().getHeight(), 2u)
  TEST_REAL_SIMILAR(mie.getGeometry().getPixelSizeX(), 20.0)
  TEST_EQUAL(mie.getGeometry().getNumberOfRegions(), 1u)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(0, 0), 2u) // the 500.0-peak spectrum moved to index 0 ... 
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(0, 1), 0u) // ... and pixel (0,1)'s spectrum is now at 0
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(1, 0), 1u)
  TEST_REAL_SIMILAR(mie.getSpectrum(0, 0)[0].getMZ(), 499.95) // pixel (0,0) still yields its own spectrum
  TEST_REAL_SIMILAR(mie.getSpectrum(0, 1)[0].getMZ(), 500.0)

  // erasing a spectrum drops its pixel; the others are re-indexed
  exp.getSpectra().erase(exp.getSpectra().begin() + 1); // spectrum 1 = pixel (1,0)
  mie.rebuildGeometry();
  mie.validate();
  TEST_EQUAL(mie.getNumberOfPixels(), 2u)
  TEST_EQUAL(mie.hasPixel(1, 0), false)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(0, 1), 0u)
  TEST_EQUAL(mie.getGeometry().getSpectrumIndex(0, 0), 1u)

  // no declared grid: dimensions grow to the bounding box of the mapped pixels
  MSImagingExperiment grown(makeSpecExperiment());
  MSImagingExperiment::setPixelCoordinate(grown[0], 4, 1);
  MSImagingExperiment::setPixelCoordinate(grown[2], 0, 3);
  grown.rebuildGeometry();
  TEST_EQUAL(grown.getNumberOfPixels(), 2u)
  TEST_EQUAL(grown.getGeometry().getWidth(), 5u)
  TEST_EQUAL(grown.getGeometry().getHeight(), 4u)
  TEST_EQUAL(grown.getGeometry().getSpectrumIndex(4, 1), 0u)

  // duplicates: first spectrum per pixel wins; out-of-grid and non-z1 spectra are skipped
  MSImagingExperiment dup(makeSpecExperiment());
  dup.getGeometry().setDimensions(2, 2);
  MSImagingExperiment::setPixelCoordinate(dup[0], 0, 0);
  MSImagingExperiment::setPixelCoordinate(dup[1], 0, 0);
  MSImagingExperiment::setPixelCoordinate(dup[2], 5, 5);
  dup.rebuildGeometry();
  TEST_EQUAL(dup.getNumberOfPixels(), 1u)
  TEST_EQUAL(dup.getGeometry().getSpectrumIndex(0, 0), 0u)
  dup[2].setMetaValue(MSImagingExperiment::META_PIXEL_X, 2);
  dup[2].setMetaValue(MSImagingExperiment::META_PIXEL_Y, 2);
  dup[2].setMetaValue(MSImagingExperiment::META_PIXEL_Z, 2);
  dup.rebuildGeometry();
  TEST_EQUAL(dup.hasPixel(1, 1), false)
}
END_SECTION

START_SECTION((static void bindPixelsFromSpectra(const MSExperiment& exp, MSImagingGeometry& geom)))
{
  MSImagingExperiment fx = makeFixture(); // spectra annotated by setGeometry()
  MSImagingGeometry geom;
  MSImagingExperiment::bindPixelsFromSpectra(fx.getMSExperiment(), geom);
  TEST_EQUAL(geom.getNumberOfPixels(), 3u)
  TEST_EQUAL(geom.getWidth(), 2u)
  TEST_EQUAL(geom.getHeight(), 2u)
  TEST_EQUAL(geom.getSpectrumIndex(0, 1), 2u)
}
END_SECTION

START_SECTION((static bool getPixelCoordinate(const MSSpectrum& spec, UInt& x, UInt& y)))
{
  MSSpectrum s;
  UInt x = 7, y = 7;
  TEST_EQUAL(MSImagingExperiment::getPixelCoordinate(s, x, y), false)
  s.setMetaValue(MSImagingExperiment::META_PIXEL_X, 3);
  s.setMetaValue(MSImagingExperiment::META_PIXEL_Y, 1);
  TEST_EQUAL(MSImagingExperiment::getPixelCoordinate(s, x, y), true)
  TEST_EQUAL(x, 2u) // 1-based meta value -> 0-based geometry coordinate
  TEST_EQUAL(y, 0u)
  s.setMetaValue(MSImagingExperiment::META_PIXEL_Z, 2);
  TEST_EQUAL(MSImagingExperiment::getPixelCoordinate(s, x, y), false) // only plane 1 is modeled
  s.setMetaValue(MSImagingExperiment::META_PIXEL_Z, 1);
  s.setMetaValue(MSImagingExperiment::META_PIXEL_X, 0);
  TEST_EQUAL(MSImagingExperiment::getPixelCoordinate(s, x, y), false) // non-conformant (< 1)
}
END_SECTION

START_SECTION((static void setPixelCoordinate(MSSpectrum& spec, UInt x, UInt y)))
{
  MSSpectrum s;
  MSImagingExperiment::setPixelCoordinate(s, 0, 4);
  TEST_EQUAL((Int)s.getMetaValue(MSImagingExperiment::META_PIXEL_X), 1)
  TEST_EQUAL((Int)s.getMetaValue(MSImagingExperiment::META_PIXEL_Y), 5)
  TEST_EQUAL((Int)s.getMetaValue(MSImagingExperiment::META_PIXEL_Z), 1)
  UInt x = 0, y = 0;
  TEST_EQUAL(MSImagingExperiment::getPixelCoordinate(s, x, y), true)
  TEST_EQUAL(x, 0u)
  TEST_EQUAL(y, 4u)
}
END_SECTION

START_SECTION((std::vector<Size> getRegionSpectrumIndices(Size region_id) const))
{
  MSImagingExperiment mie = makeFixture();
  // region covering (0,0) and (0,1) -> spectrum indices 0 and 2
  mie.getGeometry().addRegion(MSImagingRegion::rectangle(1, "col0", 0, 0, 0, 1));
  auto spec = mie.getRegionSpectrumIndices(1);
  TEST_EQUAL(spec.size(), 2u)
  TEST_EQUAL(spec[0], 0u) // pixel (0,0) -> spectrum 0
  TEST_EQUAL(spec[1], 2u) // pixel (0,1) -> spectrum 2
  TEST_EXCEPTION(Exception::ElementNotFound, mie.getRegionSpectrumIndices(99))
}
END_SECTION


START_SECTION((IonImage extractIonImage(double mz, double tolerance_ppm, Size region_id) const))
{
  MSImagingExperiment mie = makeFixture();
  mie.getGeometry().addRegion(MSImagingRegion::rectangle(1, "col0", 0, 0, 0, 1)); // (0,0) and (0,1) only
  IonImage img = mie.extractIonImage(500.0, 200.0, 1);

  // global dims preserved
  TEST_EQUAL(img.getWidth(), 2u)
  TEST_EQUAL(img.getHeight(), 2u)

  // region pixels populated, same values as whole-image
  TEST_EQUAL(img.hasPixel(0, 0), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 0), 30.0)
  TEST_EQUAL(img.hasPixel(0, 1), true)
  TEST_REAL_SIMILAR(img.getIntensity(0, 1), 100.0)

  // THE KEY ASSERT: (1,0) is acquired but OUTSIDE the region -> must stay masked/invalid
  TEST_EQUAL(img.hasPixel(1, 0), false)

  // unknown region id
  TEST_EXCEPTION(Exception::ElementNotFound, mie.extractIonImage(500.0, 200.0, 99))
}
END_SECTION
/////////////////////////////////////////////////////////////
END_TEST
