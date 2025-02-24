// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MobilityPeak1D.h>

///////////////////////////

START_TEST(MobilityPeak1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MobilityPeak1D* d10_ptr = nullptr;
MobilityPeak1D* d10_nullPointer = nullptr;

static_assert(std::is_trivially_destructible<MobilityPeak1D> {});
//static_assert(std::is_trivially_default_constructible<MobilityPeak1D> {});
static_assert(std::is_trivially_copy_constructible<MobilityPeak1D> {});
static_assert(std::is_trivially_copy_assignable<MobilityPeak1D> {});
static_assert(std::is_trivially_move_constructible<MobilityPeak1D> {});
static_assert(std::is_nothrow_move_constructible<MobilityPeak1D> {});
static_assert(std::is_trivially_move_assignable<MobilityPeak1D> {});

START_SECTION((MobilityPeak1D()))
  d10_ptr = new MobilityPeak1D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~MobilityPeak1D()))
  delete d10_ptr;
END_SECTION

START_SECTION((IntensityType getIntensity() const))
  TEST_REAL_SIMILAR(MobilityPeak1D().getIntensity(), 0.0)
END_SECTION

START_SECTION((PositionType const& getPosition() const))
  TEST_REAL_SIMILAR(MobilityPeak1D().getPosition()[0], 0.0)
END_SECTION

START_SECTION((CoordinateType getMobility() const))
  TEST_REAL_SIMILAR(MobilityPeak1D().getMobility(), 0.0)
END_SECTION

START_SECTION((CoordinateType getPos() const))
  TEST_REAL_SIMILAR(MobilityPeak1D().getPos(), 0.0)
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
  MobilityPeak1D p;
  p.setIntensity(17.8f);
  TEST_REAL_SIMILAR(p.getIntensity(), 17.8)
END_SECTION

START_SECTION((void setPosition(PositionType const& position)))
  MobilityPeak1D::PositionType pos;
  pos[0] = 1.0;
  MobilityPeak1D p;
  p.setPosition(pos);
  TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((PositionType & getPosition()))
  MobilityPeak1D::PositionType pos;
  pos[0] = 1.0;
  MobilityPeak1D p;
  p.getPosition() = pos;
  TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((void setMobility(CoordinateType mb)))
  MobilityPeak1D p;
  p.setMobility(5.0);
  TEST_REAL_SIMILAR(p.getMobility(), 5.0)
END_SECTION

START_SECTION((void setPos(CoordinateTypepos)))
  MobilityPeak1D p;
  p.setPos(5.0);
  TEST_REAL_SIMILAR(p.getPos(), 5.0)
END_SECTION

START_SECTION((MobilityPeak1D(const MobilityPeak1D& p)))
  MobilityPeak1D::PositionType pos;
  pos[0] = 21.21;
  MobilityPeak1D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  MobilityPeak1D::PositionType pos2;
  MobilityPeak1D::IntensityType i2;

  MobilityPeak1D copy_of_p(p);

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((MobilityPeak1D & operator=(const MobilityPeak1D& rhs)))
  MobilityPeak1D::PositionType pos;
  pos[0] = 21.21;
  MobilityPeak1D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  MobilityPeak1D::PositionType pos2;
  MobilityPeak1D::IntensityType i2;

  MobilityPeak1D copy_of_p;
  copy_of_p = p;

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((bool operator==(const MobilityPeak1D& rhs) const))
  MobilityPeak1D p1;
  MobilityPeak1D p2(p1);
  TEST_TRUE(p1 == p2)

  p1.setIntensity(5.0f);
  TEST_FALSE(p1 == p2)
  p2.setIntensity(5.0f);
  TEST_TRUE(p1 == p2)

  p1.getPosition()[0] = 5;
  TEST_FALSE(p1 == p2)
  p2.getPosition()[0] = 5;
  TEST_TRUE(p1 == p2)
END_SECTION

START_SECTION((bool operator!=(const MobilityPeak1D& rhs) const))
  MobilityPeak1D p1;
  MobilityPeak1D p2(p1);
  TEST_FALSE(p1 != p2)

  p1.setIntensity(5.0f);
  TEST_TRUE(p1 != p2)
  p2.setIntensity(5.0f);
  TEST_FALSE(p1 != p2)

  p1.getPosition()[0] = 5;
  TEST_TRUE(p1 != p2)
  p2.getPosition()[0] = 5;
  TEST_FALSE(p1 != p2)
END_SECTION


/////////////////////////////////////////////////////////////
// Nested stuff
/////////////////////////////////////////////////////////////

MobilityPeak1D p1;
p1.setIntensity(10.0);
p1.setMobility(10.0);
MobilityPeak1D p2;
p2.setIntensity(12.0);
p2.setMobility(12.0);

// IntensityLess
START_SECTION(([MobilityPeak1D::IntensityLess] bool operator()(MobilityPeak1D const& left, MobilityPeak1D const& right) const))

  std::vector<MobilityPeak1D> v;
  MobilityPeak1D p;

  p.setIntensity(2.5f);
  v.push_back(p);

  p.setIntensity(3.5f);
  v.push_back(p);

  p.setIntensity(1.5f);
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak1D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  v[0] = v[2];
  v[2] = p;
  std::sort(v.begin(), v.end(), MobilityPeak1D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  // some more
  TEST_TRUE(MobilityPeak1D::IntensityLess()(p1, p2))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2, p1))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2, p2))

END_SECTION

START_SECTION(([MobilityPeak1D::IntensityLess] bool operator()(MobilityPeak1D const& left, IntensityType right) const))
  TEST_TRUE(MobilityPeak1D::IntensityLess()(p1, p2.getIntensity()))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2, p1.getIntensity()))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2, p2.getIntensity()))
END_SECTION

START_SECTION(([MobilityPeak1D::IntensityLess] bool operator()(IntensityType left, MobilityPeak1D const& right) const))
  TEST_TRUE(MobilityPeak1D::IntensityLess()(p1.getIntensity(), p2))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2.getIntensity(), p1))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2.getIntensity(), p2))
END_SECTION

START_SECTION(([MobilityPeak1D::IntensityLess] bool operator()(IntensityType left, IntensityType right) const))
  TEST_TRUE(MobilityPeak1D::IntensityLess()(p1.getIntensity(), p2.getIntensity()))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2.getIntensity(), p1.getIntensity()))
  TEST_FALSE(MobilityPeak1D::IntensityLess()(p2.getIntensity(), p2.getIntensity()))
END_SECTION

// MobilityLess
START_SECTION(([MobilityPeak1D::MobilityLess] bool operator()(const MobilityPeak1D& left, const MobilityPeak1D& right) const))

  std::vector<MobilityPeak1D> v;
  MobilityPeak1D p;

  p.setMobility(3.0);
  v.push_back(p);

  p.setMobility(2.0);
  v.push_back(p);

  p.setMobility(1.0);
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak1D::MobilityLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  //
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak1D::MobilityLess] bool operator()(MobilityPeak1D const& left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p1, p2.getMobility()), true)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2, p1.getMobility()), false)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2, p2.getMobility()), false)

END_SECTION

START_SECTION(([MobilityPeak1D::MobilityLess] bool operator()(CoordinateType left, MobilityPeak1D const& right) const))

  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p1.getMobility(), p2), true)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2.getMobility(), p1), false)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2.getMobility(), p2), false)

END_SECTION

START_SECTION(([MobilityPeak1D::MobilityLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p1.getMobility(), p2.getMobility()), true)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2.getMobility(), p1.getMobility()), false)
  TEST_EQUAL(MobilityPeak1D::MobilityLess()(p2.getMobility(), p2.getMobility()), false)

END_SECTION

// PositionLess
START_SECTION(([MobilityPeak1D::PositionLess] bool operator()(const MobilityPeak1D& left, const MobilityPeak1D& right) const))

  std::vector<MobilityPeak1D> v;
  MobilityPeak1D p;

  p.getPosition()[0] = 3.0;
  v.push_back(p);

  p.getPosition()[0] = 2.0;
  v.push_back(p);

  p.getPosition()[0] = 1.0;
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak1D::PositionLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  //
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak1D::PositionLess] bool operator()(const MobilityPeak1D& left, const PositionType& right) const))
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p1, p2.getPosition()), true)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2, p1.getPosition()), false)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2, p2.getPosition()), false)
END_SECTION

START_SECTION(([MobilityPeak1D::PositionLess] bool operator()(const PositionType& left, const MobilityPeak1D& right) const))
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p1.getPosition(), p2), true)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2.getPosition(), p1), false)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2.getPosition(), p2), false)
END_SECTION

START_SECTION(([MobilityPeak1D::PositionLess] bool operator()(const PositionType& left, const PositionType& right) const))
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p1.getPosition(), p2.getPosition()), true)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2.getPosition(), p1.getPosition()), false)
  TEST_EQUAL(MobilityPeak1D::PositionLess()(p2.getPosition(), p2.getPosition()), false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
