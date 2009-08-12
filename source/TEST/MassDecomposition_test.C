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
#include <OpenMS/ANALYSIS/DENOVO/MassDecomposition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassDecomposition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassDecomposition* ptr = 0;
START_SECTION(MassDecomposition())
{
	ptr = new MassDecomposition();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~MassDecomposition())
{
	delete ptr;
}
END_SECTION

START_SECTION((MassDecomposition(const MassDecomposition &deco)))
{
  // TODO
}
END_SECTION

START_SECTION((MassDecomposition(const String &deco)))
{
  // TODO
}
END_SECTION

START_SECTION((MassDecomposition& operator=(const MassDecomposition &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((MassDecomposition& operator+=(const MassDecomposition &d)))
{
  // TODO
}
END_SECTION

START_SECTION((String toString() const ))
{
  // TODO
}
END_SECTION

START_SECTION((String toExpandedString() const ))
{
  // TODO
}
END_SECTION

START_SECTION((MassDecomposition operator+(const MassDecomposition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((UInt getNumberOfMaxAA() const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator<(const MassDecomposition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const String &deco) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool containsTag(const String &tag) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool compatible(const MassDecomposition &deco) const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



