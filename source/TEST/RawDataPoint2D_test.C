// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RawDataPoint2D.h>

///////////////////////////

START_TEST(RawDataPoint2D<D>, "$Id: DRawDataPoint_test.C 1300 2007-01-18 07:27:04Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RawDataPoint2D* d10_ptr = 0;
CHECK(RawDataPoint2D())
	d10_ptr = new RawDataPoint2D;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~RawDataPoint2D())
	delete d10_ptr;
RESULT

CHECK(const IntensityType& getIntensity() const)
	const RawDataPoint2D p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(const PositionType& getPos() const)
	const RawDataPoint2D	p;
	TEST_REAL_EQUAL(p.getPos()[0], 0.0)
	TEST_REAL_EQUAL(p.getPos()[1], 0.0)
RESULT


CHECK(void getRT(CoordinateType const& coordinate))
	RawDataPoint2D p0;
	const RawDataPoint2D p1(p0);
	TEST_REAL_EQUAL(p1.getRT(), 0.0)
RESULT

CHECK(void setRT(Size const i, const CoordinateType& coordinate))
	RawDataPoint2D p0;
  p0.setRT(12345.0);
	const RawDataPoint2D p1(p0);
	TEST_REAL_EQUAL(p1.getRT(), 12345.0)
RESULT

CHECK(void getMZ(CoordinateType const& coordinate))
	RawDataPoint2D p0;
	const RawDataPoint2D p1(p0);
	TEST_REAL_EQUAL(p1.getMZ(), 0.0)
RESULT

CHECK(void setMZ(Size const i, const CoordinateType& coordinate))
	RawDataPoint2D p0;
  p0.setMZ(12345.0);
	const RawDataPoint2D p1(p0);
	TEST_REAL_EQUAL(p1.getMZ(), 12345.0)
RESULT

CHECK(void setPos(PositionType const& position))
	RawDataPoint2D p0;
  p0.getPos()[0] = 876;
  p0.getPos()[1] = 12345.0;
	RawDataPoint2D p1;
  p1.setPos(p0.getPos());
	TEST_REAL_EQUAL(p1.getPos()[0], 876)
	TEST_REAL_EQUAL(p1.getPos()[1], 12345.0)
RESULT

CHECK(CoordinateType const& getMZ() const)
	RawDataPoint2D p0;
  p0.setPos(12345.0);
	const RawDataPoint2D p1(p0);
	TEST_REAL_EQUAL(p1.getPos()[0], 12345.0)
RESULT

CHECK(IntensityType& getIntensity())
	RawDataPoint2D p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.setIntensity(123.456);
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.setIntensity(-0.12345);
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.setIntensity(0.0);
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(void setIntensity(const IntensityType& intensity))
	RawDataPoint2D p;
  p.setIntensity(17.8);
  TEST_REAL_EQUAL(p.getIntensity(), 17.8)
RESULT


CHECK(PositionType& getPos())
	RawDataPoint2D::PositionType pos;
	RawDataPoint2D p;
	pos = p.getPos();
	TEST_REAL_EQUAL(pos[0], 0.0)
	TEST_REAL_EQUAL(pos[1], 0.0)
	pos[0] = 1.0;
	pos[1] = 2.0;
	p.setPos(pos);
	RawDataPoint2D::PositionType pos2(p.getPos());
	TEST_REAL_EQUAL(pos2[0], 1.0)
	TEST_REAL_EQUAL(pos2[1], 2.0)
RESULT

CHECK(RawDataPoint2D(const DRawDataPoint& p))
	RawDataPoint2D::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	RawDataPoint2D p;
	p.setIntensity(123.456);
	p.setPos(pos);
	RawDataPoint2D::PositionType pos2;
	RawDataPoint2D::IntensityType i2;

	RawDataPoint2D copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
RESULT

CHECK(RawDataPoint2D& operator = (const DRawDataPoint& rhs))
	RawDataPoint2D::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	RawDataPoint2D p;
	p.setIntensity(123.456);
	p.setPos(pos);
	RawDataPoint2D::PositionType pos2;
	RawDataPoint2D::IntensityType i2;

	RawDataPoint2D copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
RESULT

CHECK(bool operator == (const RawDataPoint2D& rhs) const)
	RawDataPoint2D p1;
	RawDataPoint2D p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
RESULT

CHECK(bool operator != (const RawDataPoint2D& rhs) const)
	RawDataPoint2D p1;
	RawDataPoint2D p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
RESULT

CHECK([EXTRA] class PositionLess)
	std::vector<RawDataPoint2D > v;
	RawDataPoint2D p;
	
	p.getPos()[0]=3.0;
	p.getPos()[1]=2.5;
	v.push_back(p);

	p.getPos()[0]=2.0;
	p.getPos()[1]=3.5;
	v.push_back(p);

	p.getPos()[0]=1.0;
	p.getPos()[1]=1.5;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), RawDataPoint2D::PositionLess());
	TEST_REAL_EQUAL(v[0].getPos()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPos()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPos()[0], 3.0)
	TEST_REAL_EQUAL(v[0].getPos()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPos()[1], 3.5)
	TEST_REAL_EQUAL(v[2].getPos()[1], 2.5)

	std::sort(v.begin(), v.end(), RawDataPoint2D::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPos()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPos()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPos()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPos()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPos()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPos()[0], 2.0)
RESULT

CHECK([EXTRA] struct NthPositionLess)

	std::vector<RawDataPoint2D > v;
	RawDataPoint2D p;
	
	p.getPos()[0]=3.0;
	p.getPos()[1]=2.5;
	v.push_back(p);

	p.getPos()[0]=2.0;
	p.getPos()[1]=3.5;
	v.push_back(p);

	p.getPos()[0]=1.0;
	p.getPos()[1]=1.5;
	v.push_back(p);

	std::sort(v.begin(), v.end(), RawDataPoint2D::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPos()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPos()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPos()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPos()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPos()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPos()[0], 2.0)
RESULT

CHECK([EXTRA] struct IntensityLess)
	std::vector<RawDataPoint2D > v;
	RawDataPoint2D p;
	
	p.setIntensity(2.5);
	v.push_back(p);

	p.setIntensity(3.5);
	v.push_back(p);

	p.setIntensity(1.5);
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), RawDataPoint2D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)

	v[0]=v[2];
	v[2]=p;
	std::sort(v.begin(), v.end(), RawDataPoint2D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)
RESULT

CHECK(enum value RawDataPoint2D::RT)
{
	TEST_EQUAL(RawDataPoint2D::RT,0);
}
RESULT

CHECK(enum value RawDataPoint2D::MZ)
{
	TEST_EQUAL(RawDataPoint2D::MZ,1);
}
RESULT

CHECK(enum value RawDataPoint2D::DIMENSION)
{
	TEST_EQUAL(RawDataPoint2D::DIMENSION,2);
}
RESULT

CHECK(enum RawDataPoint2D::DimensionId)
{
	RawDataPoint2D::DimensionDescription dim;
	dim = RawDataPoint2D::RT;
	TEST_EQUAL(dim,RawDataPoint2D::RT);
	dim = RawDataPoint2D::MZ;
	TEST_EQUAL(dim,RawDataPoint2D::MZ);
	dim = RawDataPoint2D::DIMENSION;
	TEST_EQUAL(dim,RawDataPoint2D::DIMENSION);
}
RESULT

CHECK(static char const * const shortDimensionName(DimensionId const dim))
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionName(RawDataPoint2D::RT),"RT");
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionName(RawDataPoint2D::MZ),"MZ");
}
RESULT

CHECK(static char const * const shortDimensionNameRT())
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionNameRT(),"RT");
}
RESULT

CHECK(static char const * const shortDimensionNameMZ())
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionNameMZ(),"MZ");
}
RESULT

CHECK(static char const * const fullDimensionName(DimensionId const dim))
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionName(RawDataPoint2D::RT),"retention time");
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionName(RawDataPoint2D::MZ),"mass-to-charge");
}
RESULT

CHECK(static char const * const fullDimensionNameRT())
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionNameRT(),"retention time");
}
RESULT

CHECK(static char const * const fullDimensionNameMZ())
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionNameMZ(),"mass-to-charge");
}
RESULT

CHECK(static char const * const shortDimensionUnit(DimensionId const dim))
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::RT),"sec");
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::MZ),"Th");
}
RESULT

CHECK(static char const * const shortDimensionUnitRT())
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionUnitRT(),"sec");
}
RESULT

CHECK(static char const * const shortDimensionUnitMZ())
{
	TEST_STRING_EQUAL(RawDataPoint2D::shortDimensionUnitMZ(),"Th");
}
RESULT

CHECK(static char const * const fullDimensionUnit(DimensionId const dim))
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionUnit(RawDataPoint2D::RT),"Seconds");
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionUnit(RawDataPoint2D::MZ),"Thomson");
}
RESULT

CHECK(static char const * const fullDimensionUnitRT())
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionUnitRT(),"Seconds");
}
RESULT

CHECK(static char const * const fullDimensionUnitMZ())
{
	TEST_STRING_EQUAL(RawDataPoint2D::fullDimensionUnitMZ(),"Thomson");
}
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
