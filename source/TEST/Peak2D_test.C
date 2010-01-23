// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(Peak2D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

Peak2D* d10_ptr = 0;
START_SECTION((Peak2D()))
	d10_ptr = new Peak2D;
	TEST_NOT_EQUAL(d10_ptr, 0)
END_SECTION

START_SECTION((~Peak2D()))
	delete d10_ptr;
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

START_SECTION((Peak2D(const Peak2D &p)))
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

START_SECTION((bool operator == (const Peak2D& rhs) const))
	Peak2D p1;
	Peak2D p2(p1);
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

START_SECTION((bool operator != (const Peak2D& rhs) const))
	Peak2D p1;
	Peak2D p2(p1);
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
END_SECTION

START_SECTION([EXTRA] struct MZLess)

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
END_SECTION

START_SECTION([EXTRA] struct RTLess)

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
END_SECTION

START_SECTION([EXTRA] struct IntensityLess)
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
/////////////////////////////////////////////////////////////
END_TEST
