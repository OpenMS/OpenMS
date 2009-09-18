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
#include <OpenMS/KERNEL/MSChromatogram.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSChromatogram, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSChromatogram<ChromatogramPeak>* ptr = 0;
START_SECTION(MSChromatogram())
{
	ptr = new MSChromatogram<ChromatogramPeak>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~MSChromatogram())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const FloatDataArrays& getFloatDataArrays() const ))
{
  // TODO
}
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
{
  // TODO
}
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const ))
{
  // TODO
}
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
{
  // TODO
}
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const ))
{
  // TODO
}
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
{
  // TODO
}
END_SECTION

START_SECTION((void sortByIntensity(bool reverse=false)))
{
  // TODO
}
END_SECTION

START_SECTION((void sortByPosition()))
{
  // TODO
}
END_SECTION

START_SECTION((bool isSorted() const ))
{
  // TODO
}
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZEnd(CoordinateType mz)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  // TODO
}
END_SECTION

START_SECTION((MSChromatogram(const MSChromatogram &source)))
{
  // TODO
}
END_SECTION

START_SECTION((~MSChromatogram()))
{
  // TODO
}
END_SECTION

START_SECTION((MSChromatogram& operator=(const MSChromatogram &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const MSChromatogram &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const MSChromatogram &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void updateRanges()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



