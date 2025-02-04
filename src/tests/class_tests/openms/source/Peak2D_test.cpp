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

#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(Peak2D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

static_assert(std::is_trivially_destructible<Peak2D> {});
// static_assert(std::is_trivially_default_constructible<Peak2D> {});
static_assert(std::is_trivially_copy_constructible<Peak2D> {});
static_assert(std::is_trivially_copy_assignable<Peak2D> {});
static_assert(std::is_trivially_move_constructible<Peak2D> {});
static_assert(std::is_nothrow_move_constructible<Peak2D> {});
static_assert(std::is_trivially_move_assignable<Peak2D> {});


Peak2D* d10_ptr = nullptr;
Peak2D* d10_nullPointer = nullptr;

START_SECTION((Peak2D()))
{
  d10_ptr = new Peak2D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
}
END_SECTION

START_SECTION((~Peak2D()))
{
  delete d10_ptr;
}
END_SECTION


START_SECTION((Peak2D(const Peak2D &p)))
{
  Peak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Peak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  Peak2D::PositionType pos2;
  Peak2D::IntensityType i2;

  Peak2D copy_of_p(p);

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
}
END_SECTION

START_SECTION((Peak2D(Peak2D &&rhs)))
{
  // Ensure that Peak2D has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(Peak2D(std::declval<Peak2D&&>())), true)

  Peak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Peak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  Peak2D::PositionType pos2;
  Peak2D::IntensityType i2;

  Peak2D copy_of_p(std::move(p));

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
}
END_SECTION

START_SECTION((explicit Peak2D(const PositionType& pos, const IntensityType in)))
{
  Peak2D p(Peak2D::PositionType(21.21, 22.22), 123.456f);
  Peak2D copy_of_p(p);
  TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[0], 21.21)
  TEST_REAL_SIMILAR(copy_of_p.getPosition()[1], 22.22)
}
END_SECTION

START_SECTION((Peak2D& operator=(const Peak2D &rhs)))
  Peak2D::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  Peak2D p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  Peak2D::PositionType pos2;
  Peak2D::IntensityType i2;

  Peak2D copy_of_p;
  copy_of_p = p;

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
END_SECTION

START_SECTION((IntensityType getIntensity() const))
	TEST_REAL_SIMILAR(Peak2D().getIntensity(), 0.0)
END_SECTION

START_SECTION((PositionType const& getPosition() const))
	const Peak2D	p;
	TEST_REAL_SIMILAR(p.getPosition()[0], 0.0)
	TEST_REAL_SIMILAR(p.getPosition()[1], 0.0)
END_SECTION

START_SECTION((CoordinateType getRT() const))
	TEST_REAL_SIMILAR(Peak2D().getRT(), 0.0)
END_SECTION

START_SECTION((CoordinateType getMZ() const))
	TEST_REAL_SIMILAR(Peak2D().getMZ(), 0.0)
END_SECTION

START_SECTION((void setRT(CoordinateTypecoordinate)))
	Peak2D p0;
  	p0.setRT(12345.0);
	TEST_REAL_SIMILAR(p0.getRT(), 12345.0)
END_SECTION

START_SECTION((void setMZ(CoordinateTypecoordinate)))
	Peak2D p0;
  	p0.setMZ(12345.0);
	TEST_REAL_SIMILAR(p0.getMZ(), 12345.0)
END_SECTION

START_SECTION((void setPosition(const PositionType &position)))
	DPosition<2> p;
  	p[0] = 876;
  	p[1] = 12345.0;
	Peak2D p1;
  	p1.setPosition(p);
	TEST_REAL_SIMILAR(p1.getPosition()[0], 876)
	TEST_REAL_SIMILAR(p1.getPosition()[1], 12345.0)
END_SECTION

START_SECTION((PositionType& getPosition()))
	DPosition<2> p;
  	p[0] = 876;
  	p[1] = 12345.0;
	Peak2D p1;
  	p1.getPosition() = p;
	TEST_REAL_SIMILAR(p1.getPosition()[0], 876)
	TEST_REAL_SIMILAR(p1.getPosition()[1], 12345.0)
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
	Peak2D p;
  p.setIntensity(17.8f);
  TEST_REAL_SIMILAR(p.getIntensity(), 17.8)
END_SECTION


START_SECTION((bool operator == (const Peak2D& rhs) const))
	Peak2D p1;
	Peak2D p2(p1);
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

START_SECTION((bool operator != (const Peak2D& rhs) const))
	Peak2D p1;
	Peak2D p2(p1);
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

START_SECTION(([EXTRA]enum value Peak2D::RT))
{
	TEST_EQUAL(Peak2D::RT,0);
}
END_SECTION

START_SECTION(([EXTRA]enum value Peak2D::MZ))
{
	TEST_EQUAL(Peak2D::MZ,1);
}
END_SECTION

START_SECTION(([EXTRA]enum value Peak2D::DIMENSION))
{
	TEST_EQUAL(Peak2D::DIMENSION,2);
}
END_SECTION

START_SECTION(([EXTRA]enum Peak2D::DimensionId))
{
	Peak2D::DimensionDescription dim;
	dim = Peak2D::RT;
	TEST_EQUAL(dim,Peak2D::RT);
	dim = Peak2D::MZ;
	TEST_EQUAL(dim,Peak2D::MZ);
	dim = Peak2D::DIMENSION;
	TEST_EQUAL(dim,Peak2D::DIMENSION);
}
END_SECTION

START_SECTION((static char const* shortDimensionName(UInt const dim)))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionName(Peak2D::RT),"RT");
	TEST_STRING_EQUAL(Peak2D::shortDimensionName(Peak2D::MZ),"MZ");
}
END_SECTION

START_SECTION((static char const* shortDimensionNameRT()))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionNameRT(),"RT");
}
END_SECTION

START_SECTION((static char const* shortDimensionNameMZ()))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionNameMZ(),"MZ");
}
END_SECTION

START_SECTION((static char const* fullDimensionName(UInt const dim)))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionName(Peak2D::RT),"retention time");
	TEST_STRING_EQUAL(Peak2D::fullDimensionName(Peak2D::MZ),"mass-to-charge");
}
END_SECTION

START_SECTION((static char const* fullDimensionNameRT()))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionNameRT(),"retention time");
}
END_SECTION

START_SECTION((static char const* fullDimensionNameMZ()))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionNameMZ(),"mass-to-charge");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnit(UInt const dim)))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionUnit(Peak2D::RT),"sec");
	TEST_STRING_EQUAL(Peak2D::shortDimensionUnit(Peak2D::MZ),"Th");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnitRT()))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionUnitRT(),"sec");
}
END_SECTION

START_SECTION((static char const* shortDimensionUnitMZ()))
{
	TEST_STRING_EQUAL(Peak2D::shortDimensionUnitMZ(),"Th");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnit(UInt const dim)))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionUnit(Peak2D::RT),"Seconds");
	TEST_STRING_EQUAL(Peak2D::fullDimensionUnit(Peak2D::MZ),"Thomson");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnitRT()))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionUnitRT(),"Seconds");
}
END_SECTION

START_SECTION((static char const* fullDimensionUnitMZ()))
{
	TEST_STRING_EQUAL(Peak2D::fullDimensionUnitMZ(),"Thomson");
}
END_SECTION

/////////////////////////////////////////////////////////////
// Nested Stuff
/////////////////////////////////////////////////////////////

Peak2D p1;
p1.setIntensity(10.0);
p1.setMZ(10.0);
p1.setRT(10.0);
Peak2D p2;
p2.setIntensity(12.0);
p2.setMZ(12.0);
p2.setRT(12.0);

// IntensityLess
START_SECTION(([Peak2D::IntensityLess] bool operator()(const Peak2D &left, const Peak2D &right) const))

  std::vector<Peak2D > v;
  Peak2D p;

  p.setIntensity(2.5f);
  v.push_back(p);

  p.setIntensity(3.5f);
  v.push_back(p);

  p.setIntensity(1.5f);
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak2D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  v[0]=v[2];
  v[2]=p;
  std::sort(v.begin(), v.end(), Peak2D::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  //
  TEST_EQUAL(Peak2D::IntensityLess()(p1,p2), true)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p1), false)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak2D::IntensityLess] bool operator()(const Peak2D &left, IntensityType right) const))

  TEST_EQUAL(Peak2D::IntensityLess()(p1,p2.getIntensity()), true)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p1.getIntensity()), false)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p2.getIntensity()), false)

END_SECTION

START_SECTION(([Peak2D::IntensityLess] bool operator()(IntensityType left, const Peak2D &right) const))

  TEST_EQUAL(Peak2D::IntensityLess()(p1.getIntensity(),p2), true)
  TEST_EQUAL(Peak2D::IntensityLess()(p2.getIntensity(),p1), false)
  TEST_EQUAL(Peak2D::IntensityLess()(p2.getIntensity(),p2), false)

END_SECTION

START_SECTION(([Peak2D::IntensityLess] bool operator()(IntensityType left, IntensityType right) const))

  TEST_EQUAL(Peak2D::IntensityLess()(p1,p2.getIntensity()), true)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p1.getIntensity()), false)
  TEST_EQUAL(Peak2D::IntensityLess()(p2,p2.getIntensity()), false)

END_SECTION

// RTLess
START_SECTION(([Peak2D::RTLess] bool operator()(const Peak2D &left, const Peak2D &right) const))

  std::vector<Peak2D > v;
  Peak2D p;

  p.getPosition()[0]=3.0;
  p.getPosition()[1]=2.5;
  v.push_back(p);

  p.getPosition()[0]=2.0;
  p.getPosition()[1]=3.5;
  v.push_back(p);

  p.getPosition()[0]=1.0;
  p.getPosition()[1]=1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak2D::RTLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)

  TEST_EQUAL(Peak2D::RTLess()(p1,p2), true)
  TEST_EQUAL(Peak2D::RTLess()(p2,p1), false)
  TEST_EQUAL(Peak2D::RTLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak2D::RTLess] bool operator()(const Peak2D &left, CoordinateType right) const))

  TEST_EQUAL(Peak2D::RTLess()(p1,p2.getRT()), true)
  TEST_EQUAL(Peak2D::RTLess()(p2,p1.getRT()), false)
  TEST_EQUAL(Peak2D::RTLess()(p2,p2.getRT()), false)

END_SECTION

START_SECTION(([Peak2D::RTLess] bool operator()(CoordinateType left, const Peak2D &right) const))

  TEST_EQUAL(Peak2D::RTLess()(p1.getRT(),p2), true)
  TEST_EQUAL(Peak2D::RTLess()(p2.getRT(),p1), false)
  TEST_EQUAL(Peak2D::RTLess()(p2.getRT(),p2), false)

END_SECTION

START_SECTION(([Peak2D::RTLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(Peak2D::RTLess()(p1.getRT(),p2.getRT()), true)
  TEST_EQUAL(Peak2D::RTLess()(p2.getRT(),p1.getRT()), false)
  TEST_EQUAL(Peak2D::RTLess()(p2.getRT(),p2.getRT()), false)

END_SECTION

// PositionLess
START_SECTION(([Peak2D::PositionLess] bool operator()(const Peak2D &left, const Peak2D &right) const))

  std::vector<Peak2D > v;
  Peak2D p;

  p.getPosition()[0]=3.0;
  p.getPosition()[1]=2.5;
  v.push_back(p);

  p.getPosition()[0]=2.0;
  p.getPosition()[1]=3.5;
  v.push_back(p);

  p.getPosition()[0]=1.0;
  p.getPosition()[1]=1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak2D::PositionLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 3.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 2.5)

  std::sort(v.begin(), v.end(), Peak2D::MZLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 2.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 3.5)
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 3.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 2.0)

  //
  TEST_EQUAL(Peak2D::PositionLess()(p1,p2), true)
  TEST_EQUAL(Peak2D::PositionLess()(p2,p1), false)
  TEST_EQUAL(Peak2D::PositionLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak2D::PositionLess] bool operator()(const Peak2D &left, const PositionType &right) const))

  TEST_EQUAL(Peak2D::PositionLess()(p1,p2.getPosition()), true)
  TEST_EQUAL(Peak2D::PositionLess()(p2,p1.getPosition()), false)
  TEST_EQUAL(Peak2D::PositionLess()(p2,p2.getPosition()), false)

END_SECTION

START_SECTION(([Peak2D::PositionLess] bool operator()(const PositionType &left, const Peak2D &right) const))

  TEST_EQUAL(Peak2D::PositionLess()(p1.getPosition(),p2), true)
  TEST_EQUAL(Peak2D::PositionLess()(p2.getPosition(),p1), false)
  TEST_EQUAL(Peak2D::PositionLess()(p2.getPosition(),p2), false)

END_SECTION

START_SECTION(([Peak2D::PositionLess] bool operator()(const PositionType &left, const PositionType &right) const))

  TEST_EQUAL(Peak2D::PositionLess()(p1.getPosition(),p2.getPosition()), true)
  TEST_EQUAL(Peak2D::PositionLess()(p2.getPosition(),p1.getPosition()), false)
  TEST_EQUAL(Peak2D::PositionLess()(p2.getPosition(),p2.getPosition()), false)

END_SECTION

// MZLess
START_SECTION(([Peak2D::MZLess] bool operator()(const Peak2D &left, const Peak2D &right) const))

  std::vector<Peak2D > v;
  Peak2D p;

  p.getPosition()[0]=3.0;
  p.getPosition()[1]=2.5;
  v.push_back(p);

  p.getPosition()[0]=2.0;
  p.getPosition()[1]=3.5;
  v.push_back(p);

  p.getPosition()[0]=1.0;
  p.getPosition()[1]=1.5;
  v.push_back(p);

  std::sort(v.begin(), v.end(), Peak2D::MZLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[1], 1.5)
  TEST_REAL_SIMILAR(v[1].getPosition()[1], 2.5)
  TEST_REAL_SIMILAR(v[2].getPosition()[1], 3.5)

  TEST_EQUAL(Peak2D::MZLess()(p1,p2), true)
  TEST_EQUAL(Peak2D::MZLess()(p2,p1), false)
  TEST_EQUAL(Peak2D::MZLess()(p2,p2), false)

END_SECTION

START_SECTION(([Peak2D::MZLess] bool operator()(const Peak2D &left, CoordinateType right) const))

  TEST_EQUAL(Peak2D::MZLess()(p1,p2.getMZ()), true)
  TEST_EQUAL(Peak2D::MZLess()(p2,p1.getMZ()), false)
  TEST_EQUAL(Peak2D::MZLess()(p2,p2.getMZ()), false)

END_SECTION

START_SECTION(([Peak2D::MZLess] bool operator()(CoordinateType left, const Peak2D &right) const))

  TEST_EQUAL(Peak2D::MZLess()(p1.getMZ(),p2), true)
  TEST_EQUAL(Peak2D::MZLess()(p2.getMZ(),p1), false)
  TEST_EQUAL(Peak2D::MZLess()(p2.getMZ(),p2), false)

END_SECTION

START_SECTION(([Peak2D::MZLess] bool operator()(CoordinateType left, CoordinateType right) const))

  TEST_EQUAL(Peak2D::MZLess()(p1.getMZ(),p2.getMZ()), true)
  TEST_EQUAL(Peak2D::MZLess()(p2.getMZ(),p1.getMZ()), false)
  TEST_EQUAL(Peak2D::MZLess()(p2.getMZ(),p2.getMZ()), false)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
