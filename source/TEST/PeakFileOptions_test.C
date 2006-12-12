// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/PeakFileOptions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

START_TEST(PeakFileOptions, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakFileOptions* ptr = 0;
CHECK(PeakFileOptions())
	ptr = new PeakFileOptions();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PeakFileOptions())
	delete ptr;
RESULT

CHECK(void setMetadataOnly(bool only))
	PeakFileOptions tmp;
	tmp.setMetadataOnly(true);
	TEST_EQUAL(tmp.getMetadataOnly(), true);
RESULT

CHECK(bool getMetadataOnly() const)
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMetadataOnly(), false);
RESULT

CHECK(void setRTRange(const DRange<1>& range))
	PeakFileOptions tmp;
	tmp.setRTRange(makeRange(2, 4));
	TEST_EQUAL(tmp.hasRTRange(), true);
	TEST_EQUAL(tmp.getRTRange(), makeRange(2, 4));
RESULT

CHECK(bool hasRTRange())
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasRTRange(), false);
RESULT

CHECK(const DRange<1>& getRTRange() const)
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getRTRange(), DRange<1>());
RESULT

CHECK(void setMZRange(const DRange<1>& range))
	PeakFileOptions tmp;
	tmp.setMZRange(makeRange(3, 5));
	TEST_EQUAL(tmp.hasMZRange(), true);
	TEST_EQUAL(tmp.getMZRange(), makeRange(3, 5));
RESULT

CHECK(bool hasMZRange())
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasMZRange(), false);
RESULT

CHECK(const DRange<1>& getMZRange() const)
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMZRange(), DRange<1>());
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
