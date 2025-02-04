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

#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

START_TEST(Peak1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

Peak1D* d10_ptr = nullptr;
Peak1D* d10_nullPointer = nullptr;

static_assert(std::is_trivially_destructible<Peak1D> {});
// static_assert(std::is_trivially_default_constructible<Peak1D> {});
static_assert(std::is_trivially_copy_constructible<Peak1D> {});
static_assert(std::is_trivially_copy_assignable<Peak1D> {});
static_assert(std::is_trivially_move_constructible<Peak1D> {});
static_assert(std::is_nothrow_move_constructible<Peak1D> {});
static_assert(std::is_trivially_move_assignable<Peak1D> {});

START_SECTION((Peak1D()))
	d10_ptr = new Peak1D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~Peak1D()))
	delete d10_ptr;
END_SECTION

START_SECTION((IntensityType getIntensity() const))
	TEST_REAL_SIMILAR(Peak1D().getIntensity(), 0.0)
END_SECTION

START_SECTION((PositionType const& getPosition() const))
	TEST_REAL_SIMILAR(Peak1D().getPosition()[0], 0.0)
END_SECTION

START_SECTION((CoordinateType getMZ() const))
	TEST_REAL_SIMILAR(Peak1D().getMZ(), 0.0)
END_SECTION

START_SECTION((CoordinateType getPos() const))
	TEST_REAL_SIMILAR(Peak1D().getPos(), 0.0)
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
	Peak1D p;
 	p.setIntensity(17.8f);
 	TEST_REAL_SIMILAR(p.getIntensity(), 17.8)
END_SECTION

START_SECTION((void setPosition(PositionType const &position)))
	Peak1D::PositionType pos;
	pos[0] = 1.0;
	Peak1D p;
	p.setPosition(pos);
	TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((PositionType& getPosition()))
	Peak1D::PositionType pos;
	pos[0] = 1.0;
	Peak1D p;
	p.getPosition() = pos;
	TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((void setMZ(CoordinateTypemz)))
	Peak1D p;
	p.setMZ(5.0);
	TEST_REAL_SIMILAR(p.getMZ(), 5.0)
END_SECTION

START_SECTION((void setPos(CoordinateTypepos)))
	Peak1D p;
	p.setPos(5.0);
	TEST_REAL_SIMILAR(p.getPos(), 5.0)
END_SECTION

START_SECTION((Peak1D(const Peak1D& p)))
	Peak1D::PositionType pos;
	pos[0] = 21.21;
	Peak1D p;
	p.setIntensity(123.456f);
	p.setPosition(pos);
	Peak1D::PositionType pos2;
	Peak1D::IntensityType i2;

	Peak1D copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_SIMILAR(i2, 123.456)

	TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((Peak1D& operator = (const Peak1D& rhs)))
	Peak1D::PositionType pos;
	pos[0] = 21.21;
	Peak1D p;
	p.setIntensity(123.456f);
	p.setPosition(pos);
	Peak1D::PositionType pos2;
	Peak1D::IntensityType i2;

	Peak1D copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_SIMILAR(i2, 123.456)

	TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((bool operator == (const Peak1D& rhs) const))
	Peak1D p1;
	Peak1D p2(p1);
	TEST_TRUE(p1 == p2)
	
	p1.setIntensity(5.0f);
	TEST_EQUAL(p1==p2, false)
	p2.setIntensity(5.0f);
	TEST_TRUE(p1 == p2)
	
	p1.getPosition()[0]=5;
	TEST_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_TRUE(p1 == p2)
END_SECTION

START_SECTION((bool operator != (const Peak1D& rhs) const))
	Peak1D p1;
	Peak1D p2(p1);
	TEST_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5.0f);
	TEST_FALSE(p1 == p2)
	p2.setIntensity(5.0f);
	TEST_EQUAL(p1!=p2, false)
	
	p1.getPosition()[0]=5;
	TEST_FALSE(p1 == p2)
	p2.getPosition()[0]=5;
	TEST_EQUAL(p1!=p2, false)
END_SECTION


/////////////////////////////////////////////////////////////
// Nested stuff
/////////////////////////////////////////////////////////////

Peak1D p1;
p1.setIntensity(10.0);
p1.setMZ(10.0);
Peak1D p2;
p2.setIntensity(12.0);
p2.setMZ(12.0);

// IntensityLess
START_SECTION(([Peak1D::IntensityLess] bool operator()(Peak1D const &left, Peak1D const &right) const))

  std::vector<Peak1D > v;
  Peak1D p;

  p.setIntensity(2.5f);
  v.push_back(p);

  p.setIntensity(3.5f);
  v.push_back(p);

  p.setIntensity(1.5f);
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak1D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  v[0]=v[2];
  v[2]=p;
  std::sort(v.begin(), v.end(), Peak1D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  // some more
  TEST_EQUAL(Peak1D::IntensityLess()(p1,p2), true)
  TEST_EQUAL(Peak1D::IntensityLess()(p2,p1), false)
  TEST_EQUAL(Peak1D::IntensityLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak1D::IntensityLess] bool operator()(Peak1D const &left, IntensityType right) const))

  TEST_EQUAL(Peak1D::IntensityLess()(p1,p2.getIntensity()), true)
  TEST_EQUAL(Peak1D::IntensityLess()(p2,p1.getIntensity()), false)
  TEST_EQUAL(Peak1D::IntensityLess()(p2,p2.getIntensity()), false)

END_SECTION

START_SECTION(([Peak1D::IntensityLess] bool operator()(IntensityType left, Peak1D const &right) const))

  TEST_EQUAL(Peak1D::IntensityLess()(p1.getIntensity(),p2), true)
  TEST_EQUAL(Peak1D::IntensityLess()(p2.getIntensity(),p1), false)
  TEST_EQUAL(Peak1D::IntensityLess()(p2.getIntensity(),p2), false)

END_SECTION

START_SECTION(([Peak1D::IntensityLess] bool operator()(IntensityType left, IntensityType right) const))

  TEST_EQUAL(Peak1D::IntensityLess()(p1.getIntensity(),p2.getIntensity()), true)
  TEST_EQUAL(Peak1D::IntensityLess()(p2.getIntensity(),p1.getIntensity()), false)
  TEST_EQUAL(Peak1D::IntensityLess()(p2.getIntensity(),p2.getIntensity()), false)

END_SECTION

// MZLess
START_SECTION(([Peak1D::MZLess] bool operator()(const Peak1D &left, const Peak1D &right) const))

  std::vector<Peak1D > v;
  Peak1D p;

  p.setMZ(3.0);
  v.push_back(p);

  p.setMZ(2.0);
  v.push_back(p);

  p.setMZ(1.0);
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak1D::MZLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  //
  TEST_EQUAL(Peak1D::MZLess()(p1,p2), true)
  TEST_EQUAL(Peak1D::MZLess()(p2,p1), false)
  TEST_EQUAL(Peak1D::MZLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak1D::MZLess] bool operator()(Peak1D const &left, CoordinateType right) const))

  TEST_EQUAL(Peak1D::MZLess()(p1,p2.getMZ()), true)
  TEST_EQUAL(Peak1D::MZLess()(p2,p1.getMZ()), false)
  TEST_EQUAL(Peak1D::MZLess()(p2,p2.getMZ()), false)

END_SECTION

START_SECTION(([Peak1D::MZLess] bool operator()(CoordinateType left, Peak1D const &right) const))

  TEST_EQUAL(Peak1D::MZLess()(p1.getMZ(),p2), true)
  TEST_EQUAL(Peak1D::MZLess()(p2.getMZ(),p1), false)
  TEST_EQUAL(Peak1D::MZLess()(p2.getMZ(),p2), false)

END_SECTION

START_SECTION(([Peak1D::MZLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(Peak1D::MZLess()(p1.getMZ(),p2.getMZ()), true)
  TEST_EQUAL(Peak1D::MZLess()(p2.getMZ(),p1.getMZ()), false)
  TEST_EQUAL(Peak1D::MZLess()(p2.getMZ(),p2.getMZ()), false)

END_SECTION

// PositionLess
START_SECTION(([Peak1D::PositionLess] bool operator()(const Peak1D &left, const Peak1D &right) const))

  std::vector<Peak1D > v;
  Peak1D p;

  p.getPosition()[0]=3.0;
  v.push_back(p);

  p.getPosition()[0]=2.0;
  v.push_back(p);

  p.getPosition()[0]=1.0;
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak1D::PositionLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  //
  TEST_EQUAL(Peak1D::PositionLess()(p1,p2), true)
  TEST_EQUAL(Peak1D::PositionLess()(p2,p1), false)
  TEST_EQUAL(Peak1D::PositionLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak1D::PositionLess] bool operator()(const Peak1D &left, const PositionType &right) const))

  TEST_EQUAL(Peak1D::PositionLess()(p1,p2.getPosition()), true)
  TEST_EQUAL(Peak1D::PositionLess()(p2,p1.getPosition()), false)
  TEST_EQUAL(Peak1D::PositionLess()(p2,p2.getPosition()), false)

END_SECTION

START_SECTION(([Peak1D::PositionLess] bool operator()(const PositionType &left, const Peak1D &right) const))

  TEST_EQUAL(Peak1D::PositionLess()(p1.getPosition(),p2), true)
  TEST_EQUAL(Peak1D::PositionLess()(p2.getPosition(),p1), false)
  TEST_EQUAL(Peak1D::PositionLess()(p2.getPosition(),p2), false)

END_SECTION

START_SECTION(([Peak1D::PositionLess] bool operator()(const PositionType &left, const PositionType &right) const))

  TEST_EQUAL(Peak1D::PositionLess()(p1.getPosition(),p2.getPosition()), true)
  TEST_EQUAL(Peak1D::PositionLess()(p2.getPosition(),p1.getPosition()), false)
  TEST_EQUAL(Peak1D::PositionLess()(p2.getPosition(),p2.getPosition()), false)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
