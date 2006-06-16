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
// $Id: Histogram_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

///////////////////////////

START_TEST(Distribution, "$Id: Histogram_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

CHECK(inline static double ceil_decimal(double x, int decPow))
	TEST_REAL_EQUAL(ceil_decimal(12345.671,-2),12345.68)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,-1),12345.7)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,0),12346.0)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,1),12350.0)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,2),12400.0)
RESULT

CHECK(inline static double round_decimal(double x, int decPow))
	TEST_REAL_EQUAL(round_decimal(12345.671,-2),12345.67)
	TEST_REAL_EQUAL(round_decimal(12345.67,-1),12345.7)
	TEST_REAL_EQUAL(round_decimal(12345.67,0),12346.0)
	TEST_REAL_EQUAL(round_decimal(12345.67,1),12350.0)
	TEST_REAL_EQUAL(round_decimal(12345.67,2),12300.0)
RESULT

CHECK(inline static double intervalTransformation(double x,double left1,double right1,double left2,double right2))
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.0,1.0,0.0,100.0),50.0)
RESULT 

CHECK(inline double linear2log(double x, bool is_percent=false, double max=0))

RESULT

CHECK(inline double log2linear(double x, bool is_percent=false, double max=0))

RESULT

CHECK(inline bool isOdd(UnsignedInt x))

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
