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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LowessSmoothing* ptr = 0;
LowessSmoothing* null_ptr = 0;
START_SECTION(LowessSmoothing())
{
	ptr = new LowessSmoothing();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~LowessSmoothing())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~LowessSmoothing()))
{
  // TODO
}
END_SECTION

START_SECTION((void smoothData(const DoubleVector &, const DoubleVector &, DoubleVector &)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



