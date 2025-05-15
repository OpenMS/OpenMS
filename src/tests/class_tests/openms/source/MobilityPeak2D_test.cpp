// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MobilityPeak2D.h>

///////////////////////////

START_TEST(MobilityPeak2D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

static_assert(std::is_trivially_destructible<MobilityPeak2D> {});
// static_assert(std::is_trivially_default_constructible<MobilityPeak2D> {});
static_assert(std::is_trivially_copy_constructible<MobilityPeak2D> {});
static_assert(std::is_trivially_copy_assignable<MobilityPeak2D> {});
static_assert(std::is_trivially_move_constructible<MobilityPeak2D> {});
static_assert(std::is_nothrow_move_constructible<MobilityPeak2D> {});
static_assert(std::is_trivially_move_assignable<MobilityPeak2D> {});


MobilityPeak2D* d10_ptr = nullptr;
MobilityPeak2D* d10_nullPointer = nullptr;

START_SECTION((MobilityPeak2D()))
{
  d10_ptr = new MobilityPeak2D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
}
END_SECTION

START_SECTION((~MobilityPeak2D()))
{
  delete d10_ptr;
}
END_SECTION


START_SECTION((MobilityPeak2D(const MobilityPeak2D& p)))
{
  MobilityPeak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  MobilityPeak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  MobilityPeak2D::PositionType pos2;
  MobilityPeak2D::IntensityType i2;

  MobilityPeak2D copy_of_p(p);

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
}
END_SECTION

START_SECTION((MobilityPeak2D(MobilityPeak2D && rhs)))
{
  // Ensure that MobilityPeak2D has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(MobilityPeak2D(std::declval<MobilityPeak2D&&>())), true)

  MobilityPeak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  MobilityPeak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  MobilityPeak2D::PositionType pos2;
  MobilityPeak2D::IntensityType i2;

  MobilityPeak2D copy_of_p(std::move(p));

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
}
END_SECTION

START_SECTION((explicit MobilityPeak2D(const PositionType& pos, const IntensityType in)))
{
  MobilityPeak2D p(MobilityPeak2D::PositionType(21.21, 22.22), 123.456f);
  MobilityPeak2D copy_of_p(p);
  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[0], 21.21)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[1], 22.22)
}
END_SECTION

START_SECTION((MobilityPeak2D & operator=(const MobilityPeak2D& rhs)))
  MobilityPeak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  MobilityPeak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  MobilityPeak2D::PositionType pos2;
  MobilityPeak2D::IntensityType i2;

  MobilityPeak2D copy_of_p;
  copy_of_p = p;

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
END_SECTION

START_SECTION((IntensityType getIntensity() const))
  TEST_REAL_SIMILAR(MobilityPeak2D().getIntensity(), 0.0)
END_SECTION

START_SECTION((PositionType const& getPosition() const))
  const MobilityPeak2D p {};
  TEST_REAL_SIMILAR(p.getPosition()[0], 0.0)
  TEST_REAL_SIMILAR(p.getPosition()[1], 0.0)
END_SECTION

START_SECTION((CoordinateType getMobility() const))
  TEST_REAL_SIMILAR(MobilityPeak2D().getMobility(), 0.0)
END_SECTION

START_SECTION((CoordinateType getMZ() const))
  TEST_REAL_SIMILAR(MobilityPeak2D().getMZ(), 0.0)
END_SECTION

START_SECTION((void setMobility(CoordinateTypecoordinate)))
  MobilityPeak2D p0;
  p0.setMobility(12345.0);
  TEST_REAL_SIMILAR(p0.getMobility(), 12345.0)
END_SECTION

START_SECTION((void setMZ(CoordinateTypecoordinate)))
  MobilityPeak2D p0;
  p0.setMZ(12345.0);
  TEST_REAL_SIMILAR(p0.getMZ(), 12345.0)
END_SECTION

START_SECTION((void setPosition(const PositionType& position)))
  DPosition<2> p;
  p[0] = 876;
  p[1] = 12345.0;
  MobilityPeak2D p1;
  p1.setPosition(p);
  TEST_REAL_SIMILAR(p1.getPosition()[0], 876)
  TEST_REAL_SIMILAR(p1.getPosition()[1], 12345.0)
END_SECTION

START_SECTION((PositionType & getPosition()))
DPosition<2> p;
p[0] = 876;
p[1] = 12345.0;
MobilityPeak2D p1;
p1.getPosition() = p;
TEST_REAL_SIMILAR(p1.getPosition()[0], 876)
TEST_REAL_SIMILAR(p1.getPosition()[1], 12345.0)
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
MobilityPeak2D p;
p.setIntensity(17.8f);
TEST_REAL_SIMILAR(p.getIntensity(), 17.8)
END_SECTION


START_SECTION((bool operator==(const MobilityPeak2D& rhs) const))
MobilityPeak2D p1;
MobilityPeak2D p2(p1);
TEST_TRUE(p1 == p2)

p1.setIntensity(5.0f);
TEST_EQUAL(p1 == p2, false)
p2.setIntensity(5.0f);
TEST_TRUE(p1 == p2)

p1.getPosition()[0] = 5;
TEST_EQUAL(p1 == p2, false)
p2.getPosition()[0] = 5;
TEST_TRUE(p1 == p2)
END_SECTION

START_SECTION((bool operator!=(const MobilityPeak2D& rhs) const))
MobilityPeak2D p1;
MobilityPeak2D p2(p1);
TEST_EQUAL(p1 != p2, false)

p1.setIntensity(5.0f);
TEST_FALSE(p1 == p2)
p2.setIntensity(5.0f);
TEST_EQUAL(p1 != p2, false)

p1.getPosition()[0] = 5;
TEST_FALSE(p1 == p2)
p2.getPosition()[0] = 5;
TEST_EQUAL(p1 != p2, false)
END_SECTION

START_SECTION(([EXTRA] enum value MobilityPeak2D::IM))
{
  TEST_EQUAL(MobilityPeak2D::IM, 0);
}
END_SECTION

START_SECTION(([EXTRA] enum value MobilityPeak2D::MZ))
{
  TEST_EQUAL(MobilityPeak2D::MZ, 1);
}
END_SECTION

START_SECTION(([EXTRA] enum value MobilityPeak2D::DIMENSION))
{
  TEST_EQUAL(MobilityPeak2D::DIMENSION, 2);
}
END_SECTION

START_SECTION(([EXTRA] enum MobilityPeak2D::DimensionId))
{
  MobilityPeak2D::DimensionDescription dim;
  dim = MobilityPeak2D::IM;
  TEST_EQUAL(dim, MobilityPeak2D::IM);
  dim = MobilityPeak2D::MZ;
  TEST_EQUAL(dim, MobilityPeak2D::MZ);
  dim = MobilityPeak2D::DIMENSION;
  TEST_EQUAL(dim, MobilityPeak2D::DIMENSION);
}
END_SECTION

START_SECTION((static char const* shortDimensionName(UInt const dim)))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionName(MobilityPeak2D::IM), "IM");
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionName(MobilityPeak2D::MZ), "MZ");
}
END_SECTION

START_SECTION((static char const* shortDimensionNameIM()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionNameIM(), "IM");
}
END_SECTION

START_SECTION((static char const* shortDimensionNameMZ()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionNameMZ(), "MZ");
}
END_SECTION

START_SECTION((static char const* fullDimensionName(UInt const dim)))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionName(MobilityPeak2D::IM), "ion mobility");
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionName(MobilityPeak2D::MZ), "mass-to-charge");
}
END_SECTION

START_SECTION((static char const* fullDimensionNameIM()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionNameIM(), "ion mobility");
}
END_SECTION

START_SECTION((static char const* fullDimensionNameMZ()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionNameMZ(), "mass-to-charge");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnit(UInt const dim)))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionUnit(MobilityPeak2D::IM), "?");
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionUnit(MobilityPeak2D::MZ), "Th");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnitIM()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionUnitIM(), "?");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnitMZ()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::shortDimensionUnitMZ(), "Th");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnit(UInt const dim)))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionUnit(MobilityPeak2D::IM), "?");
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionUnit(MobilityPeak2D::MZ), "Thomson");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnitIM()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionUnitIM(), "?");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnitMZ()))
{
  TEST_STRING_EQUAL(MobilityPeak2D::fullDimensionUnitMZ(), "Thomson");
}
END_SECTION

/////////////////////////////////////////////////////////////
// Nested Stuff
/////////////////////////////////////////////////////////////

MobilityPeak2D p1;
p1.setIntensity(10.0);
p1.setMZ(10.0);
p1.setMobility(10.0);
MobilityPeak2D p2;
p2.setIntensity(12.0);
p2.setMZ(12.0);
p2.setMobility(12.0);

// IntensityLess
START_SECTION(([MobilityPeak2D::IntensityLess] bool operator()(const MobilityPeak2D& left, const MobilityPeak2D& right) const))

  std::vector<MobilityPeak2D> v;
  MobilityPeak2D p;

  p.setIntensity(2.5f);
  v.push_back(p);

  p.setIntensity(3.5f);
  v.push_back(p);

  p.setIntensity(1.5f);
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak2D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  v[0] = v[2];
  v[2] = p;
  std::sort(v.begin(), v.end(), MobilityPeak2D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  //
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IntensityLess] bool operator()(const MobilityPeak2D& left, IntensityType right) const))

  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p1, p2.getIntensity()), true)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p1.getIntensity()), false)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p2.getIntensity()), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IntensityLess] bool operator()(IntensityType left, const MobilityPeak2D& right) const))

  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p1.getIntensity(), p2), true)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2.getIntensity(), p1), false)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2.getIntensity(), p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IntensityLess] bool operator()(IntensityType left, IntensityType right) const))

  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p1, p2.getIntensity()), true)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p1.getIntensity()), false)
  TEST_EQUAL(MobilityPeak2D::IntensityLess()(p2, p2.getIntensity()), false)

END_SECTION

// IMLess
START_SECTION(([MobilityPeak2D::IMLess] bool operator()(const MobilityPeak2D& left, const MobilityPeak2D& right) const))

  std::vector<MobilityPeak2D> v;
  MobilityPeak2D p;

  p.getPosition()[0] = 3.0;
  p.getPosition()[1] = 2.5;
  v.push_back(p);

  p.getPosition()[0] = 2.0;
  p.getPosition()[1] = 3.5;
  v.push_back(p);

  p.getPosition()[0] = 1.0;
  p.getPosition()[1] = 1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak2D::IMLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  TEST_EQUAL(MobilityPeak2D::IMLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IMLess] bool operator()(const MobilityPeak2D& left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak2D::IMLess()(p1, p2.getMobility()), true)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2, p1.getMobility()), false)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2, p2.getMobility()), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IMLess] bool operator()(CoordinateType left, const MobilityPeak2D& right) const))

  TEST_EQUAL(MobilityPeak2D::IMLess()(p1.getMobility(), p2), true)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2.getMobility(), p1), false)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2.getMobility(), p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::IMLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak2D::IMLess()(p1.getMobility(), p2.getMobility()), true)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2.getMobility(), p1.getMobility()), false)
  TEST_EQUAL(MobilityPeak2D::IMLess()(p2.getMobility(), p2.getMobility()), false)

END_SECTION

// PositionLess
START_SECTION(([MobilityPeak2D::PositionLess] bool operator()(const MobilityPeak2D& left, const MobilityPeak2D& right) const))

  std::vector<MobilityPeak2D> v;
  MobilityPeak2D p;

  p.getPosition()[0] = 3.0;
  p.getPosition()[1] = 2.5;
  v.push_back(p);

  p.getPosition()[0] = 2.0;
  p.getPosition()[1] = 3.5;
  v.push_back(p);

  p.getPosition()[0] = 1.0;
  p.getPosition()[1] = 1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak2D::PositionLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 3.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 2.5)

  std::sort(v.begin(), v.end(), MobilityPeak2D::MZLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 2.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 3.5)
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 3.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 2.0)

  //
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::PositionLess] bool operator()(const MobilityPeak2D& left, const PositionType& right) const))

  TEST_EQUAL(MobilityPeak2D::PositionLess()(p1, p2.getPosition()), true)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2, p1.getPosition()), false)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2, p2.getPosition()), false)

END_SECTION

START_SECTION(([MobilityPeak2D::PositionLess] bool operator()(const PositionType& left, const MobilityPeak2D& right) const))

  TEST_EQUAL(MobilityPeak2D::PositionLess()(p1.getPosition(), p2), true)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2.getPosition(), p1), false)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2.getPosition(), p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::PositionLess] bool operator()(const PositionType& left, const PositionType& right) const))

  TEST_EQUAL(MobilityPeak2D::PositionLess()(p1.getPosition(), p2.getPosition()), true)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2.getPosition(), p1.getPosition()), false)
  TEST_EQUAL(MobilityPeak2D::PositionLess()(p2.getPosition(), p2.getPosition()), false)

END_SECTION

// MZLess
START_SECTION(([MobilityPeak2D::MZLess] bool operator()(const MobilityPeak2D& left, const MobilityPeak2D& right) const))

  std::vector<MobilityPeak2D> v;
  MobilityPeak2D p;

  p.getPosition()[0] = 3.0;
  p.getPosition()[1] = 2.5;
  v.push_back(p);

  p.getPosition()[0] = 2.0;
  p.getPosition()[1] = 3.5;
  v.push_back(p);

  p.getPosition()[0] = 1.0;
  p.getPosition()[1] = 1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), MobilityPeak2D::MZLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 2.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 3.5)

  TEST_EQUAL(MobilityPeak2D::MZLess()(p1, p2), true)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2, p1), false)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2, p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::MZLess] bool operator()(const MobilityPeak2D& left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak2D::MZLess()(p1, p2.getMZ()), true)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2, p1.getMZ()), false)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2, p2.getMZ()), false)

END_SECTION

START_SECTION(([MobilityPeak2D::MZLess] bool operator()(CoordinateType left, const MobilityPeak2D& right) const))

  TEST_EQUAL(MobilityPeak2D::MZLess()(p1.getMZ(), p2), true)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2.getMZ(), p1), false)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2.getMZ(), p2), false)

END_SECTION

START_SECTION(([MobilityPeak2D::MZLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(MobilityPeak2D::MZLess()(p1.getMZ(), p2.getMZ()), true)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2.getMZ(), p1.getMZ()), false)
  TEST_EQUAL(MobilityPeak2D::MZLess()(p2.getMZ(), p2.getMZ()), false)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
