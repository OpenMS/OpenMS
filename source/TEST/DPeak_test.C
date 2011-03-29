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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/DPeak.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DPeak<1>::Type* ptr1 = 0;
DPeak<1>::Type* nullPointer1 = 0;
START_SECTION(DPeak())
{
	ptr1 = new DPeak<1>::Type();
  TEST_NOT_EQUAL(ptr1, nullPointer1);
}
END_SECTION

START_SECTION(~DPeak())
{
	delete ptr1;
}
END_SECTION

DPeak<2>::Type* ptr2 = 0;
DPeak<2>::Type* nullPointer2 = 0;
START_SECTION([EXTRA]DPeak())
{
	ptr2 = new DPeak<2>::Type();
  TEST_NOT_EQUAL(ptr2, nullPointer2);
}
END_SECTION

START_SECTION([EXTRA]~DPeak())
{
	delete ptr2;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
