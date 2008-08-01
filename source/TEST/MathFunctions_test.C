// -*- mode: C++; tab-width: 2; -*-
//  vi: set ts=2:
// 
//  --------------------------------------------------------------------------
//                    OpenMS Mass Spectrometry Framework 
//  --------------------------------------------------------------------------
//   Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
// 
//   This library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License as published by the Free Software Foundation; either
//   version 2.1 of the License, or (at your option) any later version.
// 
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
// 
//   You should have received a copy of the GNU Lesser General Public
//   License along with this library; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
//  --------------------------------------------------------------------------
//  $Maintainer: Marc Sturm $
//  --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

// /////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

// /////////////////////////

START_TEST(Distribution, "$Id$")

// ///////////////////////////////////////////////////////////

CHECK((ceil_decimal))
	TEST_REAL_EQUAL(ceil_decimal(12345.671,-2),12345.68)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,-1),12345.7)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,0),12346.0)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,1),12350.0)
	TEST_REAL_EQUAL(ceil_decimal(12345.67,2),12400.0)
RESULT

CHECK((round_decimal))
	TEST_REAL_EQUAL(round_decimal(12345.671,-2),12345.67)
	TEST_REAL_EQUAL(round_decimal(12345.67,-1),12345.7)
	TEST_REAL_EQUAL(round_decimal(12345.67,0),12346.0)
	TEST_REAL_EQUAL(round_decimal(12345.67,1),12350.0)
	TEST_REAL_EQUAL(round_decimal(12345.67,2),12300.0)
RESULT

CHECK((intervalTransformation))
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.0,1.0,0.0,600.0),300.0)
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.25,1.0,0.0,600.0),200.0)
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.0,0.75,0.0,600.0),400.0)
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.0,1.0,150.0,600.0),375.0)
	TEST_REAL_EQUAL(intervalTransformation(0.5,0.0,1.0,0.0,450.0),225.0)
RESULT 

CHECK((linear2log))
	TEST_REAL_EQUAL(linear2log(0.0),0.0)
	TEST_REAL_EQUAL(linear2log(9.0),1.0)
	TEST_REAL_EQUAL(linear2log(99.0),2.0)
	TEST_REAL_EQUAL(linear2log(999.0),3.0)
RESULT

CHECK((log2linear))
	TEST_REAL_EQUAL(log2linear(0.0),0.0)
	TEST_REAL_EQUAL(log2linear(1.0),9.0)
	TEST_REAL_EQUAL(log2linear(2.0),99.0)
	TEST_REAL_EQUAL(log2linear(3.0),999.0)
RESULT

CHECK((isOdd))
	TEST_EQUAL(isOdd(0),false)
	TEST_EQUAL(isOdd(1),true)
	TEST_EQUAL(isOdd(2),false)
	TEST_EQUAL(isOdd(3),true)
RESULT

/////////////////////////////////////////////////////////////);
/////////////////////////////////////////////////////////////
END_TEST
