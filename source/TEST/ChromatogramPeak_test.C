// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/ChromatogramPeak.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramPeak, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramPeak* ptr = 0;
ChromatogramPeak* nullPointer = 0;
START_SECTION(ChromatogramPeak())
{
	ptr = new ChromatogramPeak();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~ChromatogramPeak())
{
	delete ptr;
}
END_SECTION

START_SECTION((IntensityType getIntensity() const))
  TEST_REAL_SIMILAR(ChromatogramPeak().getIntensity(), 0.0)
END_SECTION

START_SECTION((PositionType const& getPosition() const))
  TEST_REAL_SIMILAR(ChromatogramPeak().getPosition()[0], 0.0)
END_SECTION

START_SECTION((CoordinateType getRT() const))
  TEST_REAL_SIMILAR(ChromatogramPeak().getRT(), 0.0)
END_SECTION

START_SECTION((CoordinateType getPos() const))
  TEST_REAL_SIMILAR(ChromatogramPeak().getPos(), 0.0)
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
  ChromatogramPeak p;
  p.setIntensity(17.8f);
  TEST_REAL_SIMILAR(p.getIntensity(), 17.8)
END_SECTION

START_SECTION((void setPosition(PositionType const &position)))
  ChromatogramPeak::PositionType pos;
  pos[0] = 1.0;
  ChromatogramPeak p;
  p.setPosition(pos);
  TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((PositionType& getPosition()))
  ChromatogramPeak::PositionType pos;
  pos[0] = 1.0;
  ChromatogramPeak p;
  p.getPosition() = pos;
  TEST_REAL_SIMILAR(p.getPosition()[0], 1.0)
END_SECTION

START_SECTION((void setRT(CoordinateType rt)))
  ChromatogramPeak p;
  p.setRT(5.0);
  TEST_REAL_SIMILAR(p.getRT(), 5.0)
END_SECTION

START_SECTION((void setPos(CoordinateTypepos)))
  ChromatogramPeak p;
  p.setPos(5.0);
  TEST_REAL_SIMILAR(p.getPos(), 5.0)
END_SECTION

START_SECTION((ChromatogramPeak(const ChromatogramPeak& p)))
  ChromatogramPeak::PositionType pos;
  pos[0] = 21.21;
  ChromatogramPeak p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  ChromatogramPeak::PositionType pos2;
  ChromatogramPeak::IntensityType i2;

  ChromatogramPeak copy_of_p(p);

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((ChromatogramPeak& operator = (const ChromatogramPeak& rhs)))
  ChromatogramPeak::PositionType pos;
  pos[0] = 21.21;
  ChromatogramPeak p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  ChromatogramPeak::PositionType pos2;
  ChromatogramPeak::IntensityType i2;

  ChromatogramPeak copy_of_p;
  copy_of_p = p;

  i2 = copy_of_p.getIntensity();
  pos2 = copy_of_p.getPosition();
  TEST_REAL_SIMILAR(i2, 123.456)

  TEST_REAL_SIMILAR(pos2[0], 21.21)
END_SECTION

START_SECTION((bool operator == (const ChromatogramPeak& rhs) const))
  ChromatogramPeak p1;
  ChromatogramPeak p2(p1);
  TEST_EQUAL(p1==p2, true)

  p1.setIntensity(5.0f);
  TEST_EQUAL(p1==p2, false)
  p2.setIntensity(5.0f);
  TEST_EQUAL(p1==p2, true)

  p1.getPosition()[0]=5;
  TEST_EQUAL(p1==p2, false)
  p2.getPosition()[0]=5;
  TEST_EQUAL(p1==p2, true)
END_SECTION

START_SECTION((bool operator != (const ChromatogramPeak& rhs) const))
  ChromatogramPeak p1;
  ChromatogramPeak p2(p1);
  TEST_EQUAL(p1!=p2, false)

  p1.setIntensity(5.0f);
  TEST_EQUAL(p1!=p2, true)
  p2.setIntensity(5.0f);
  TEST_EQUAL(p1!=p2, false)

  p1.getPosition()[0]=5;
  TEST_EQUAL(p1!=p2, true)
  p2.getPosition()[0]=5;
  TEST_EQUAL(p1!=p2, false)
END_SECTION

START_SECTION([EXTRA] class PositionLess)
  std::vector<ChromatogramPeak > v;
  ChromatogramPeak p;

  p.getPosition()[0]=3.0;
  v.push_back(p);

  p.getPosition()[0]=2.0;
  v.push_back(p);

  p.getPosition()[0]=1.0;
  v.push_back(p);

  std::sort(v.begin(), v.end(), ChromatogramPeak::PositionLess());
  TEST_REAL_SIMILAR(v[0].getPosition()[0], 1.0)
  TEST_REAL_SIMILAR(v[1].getPosition()[0], 2.0)
  TEST_REAL_SIMILAR(v[2].getPosition()[0], 3.0)
END_SECTION

START_SECTION([EXTRA] struct IntensityLess)
  std::vector<ChromatogramPeak > v;
  ChromatogramPeak p;

  p.setIntensity(2.5f);
  v.push_back(p);

  p.setIntensity(3.5f);
  v.push_back(p);

  p.setIntensity(1.5f);
  v.push_back(p);

  std::sort(v.begin(), v.end(), ChromatogramPeak::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)

  v[0]=v[2];
  v[2]=p;
  std::sort(v.begin(), v.end(), ChromatogramPeak::IntensityLess());
  TEST_REAL_SIMILAR(v[0].getIntensity(), 1.5)
  TEST_REAL_SIMILAR(v[1].getIntensity(), 2.5)
  TEST_REAL_SIMILAR(v[2].getIntensity(), 3.5)
END_SECTION

START_SECTION(([ChromatogramPeak::IntensityLess] bool operator()(ChromatogramPeak const &left, ChromatogramPeak const &right) const))
{
  ChromatogramPeak left,right;
  left.setIntensity(10.0);
  right.setIntensity(20.0);

  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::IntensityLess] bool operator()(ChromatogramPeak const &left, IntensityType right) const))
{
  ChromatogramPeak left;
  left.setIntensity(10.0);

  ChromatogramPeak::IntensityType right;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::IntensityLess] bool operator()(IntensityType left, ChromatogramPeak const &right) const))
{
  ChromatogramPeak::IntensityType left;
  left = 10.0;

  ChromatogramPeak right;
  right.setIntensity(20.0);

  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::IntensityLess] bool operator()(IntensityType left, IntensityType right) const))
{
  ChromatogramPeak::IntensityType left,right;
  left = 10.0;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::IntensityLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::PositionLess] bool operator()(const ChromatogramPeak &left, const ChromatogramPeak &right) const))
{
  ChromatogramPeak left,right;
  left.setPosition(10.0);
  right.setPosition(20.0);

  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::PositionLess] bool operator()(const ChromatogramPeak &left, const PositionType &right) const))
{
  ChromatogramPeak left;
  left.setPosition(10.0);

  ChromatogramPeak::PositionType right;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::PositionLess] bool operator()(const PositionType &left, const ChromatogramPeak &right) const))
{
  ChromatogramPeak::PositionType left;
  left = 10.0;
  ChromatogramPeak right;
  right.setPosition(20.0);

  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::PositionLess] bool operator()(const PositionType &left, const PositionType &right) const))
{
  ChromatogramPeak::PositionType left,right;
  left = 10.0;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::PositionLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::RTLess] bool operator()(const ChromatogramPeak &left, const ChromatogramPeak &right) const))
{
  ChromatogramPeak left;
  left.setRT(10.0);

  ChromatogramPeak right;
  right.setRT(20.0);

  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::RTLess] bool operator()(ChromatogramPeak const &left, CoordinateType right) const))
{
  ChromatogramPeak left;
  left.setRT(10.0);

  ChromatogramPeak::CoordinateType right;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::RTLess] bool operator()(CoordinateType left, ChromatogramPeak const &right) const))
{
  ChromatogramPeak::CoordinateType left;
  left = 10.0;

  ChromatogramPeak right;
  right.setRT(20.0);

  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,left), false)
}
END_SECTION

START_SECTION(([ChromatogramPeak::RTLess] bool operator()(CoordinateType left, CoordinateType right) const))
{
  ChromatogramPeak::CoordinateType left;
  left = 10.0;

  ChromatogramPeak::CoordinateType right;
  right = 20.0;

  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,right), true)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(right,left), false)
  TEST_EQUAL(ChromatogramPeak::RTLess().operator ()(left,left), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



