// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
START_SECTION(ChromatogramPeak())
{
	ptr = new ChromatogramPeak();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~ChromatogramPeak())
{
	delete ptr;
}
END_SECTION

START_SECTION((ChromatogramPeak(const ChromatogramPeak &p)))
{
  // TODO
}
END_SECTION

START_SECTION((IntensityType getIntensity() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIntensity(IntensityType intensity)))
{
  ChromatogramPeak p;
	TEST_REAL_SIMILAR(p.getIntensity(), 0)
	p.setIntensity(0.35f);
	TEST_REAL_SIMILAR(p.getIntensity(), 0.35)
}
END_SECTION

START_SECTION((CoordinateType getRT() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setRT(CoordinateType rt)))
{
  ChromatogramPeak p;
	TEST_REAL_SIMILAR(p.getRT(), 0)
	p.setRT(1.5);
	TEST_REAL_SIMILAR(p.getRT(), 1.5)
}
END_SECTION

START_SECTION((CoordinateType getPos() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPos(CoordinateType pos)))
{
  ChromatogramPeak p;
	TEST_REAL_SIMILAR(p.getPos(), 0)
	p.setPos(2.5);
	TEST_REAL_SIMILAR(p.getPos(), 2.5)
}
END_SECTION

START_SECTION((PositionType const& getPosition() const ))
{
  ChromatogramPeak p;
	//TEST_REAL_SIMILAR(p.getPosition(), 2.5)
}
END_SECTION

START_SECTION((PositionType& getPosition()))
{
  ChromatogramPeak p;
	p.getPosition() = 3.8;
}
END_SECTION

START_SECTION((void setPosition(PositionType const &position)))
{
  // TODO
}
END_SECTION

START_SECTION((ChromatogramPeak& operator=(const ChromatogramPeak &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const ChromatogramPeak &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const ChromatogramPeak &rhs) const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



