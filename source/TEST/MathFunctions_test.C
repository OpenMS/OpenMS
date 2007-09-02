// -*- Mode: C++; tab-width: 2; -*-
//  vi: set ts=2:
// 
//  --------------------------------------------------------------------------
//                    OpenMS Mass Spectrometry Framework 
//  --------------------------------------------------------------------------
//   Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

PRECISION(0.01)

CHECK(Real pearsonCorrelation(std::vector<Real> model, std::vector<Real> data))
	vector<Real> v1,v2;
	v1.push_back(1);
	v1.push_back(2);
	v1.push_back(3);
	v1.push_back(4);
	v1.push_back(5);
	
	v2.push_back(1);
	v2.push_back(2);
	v2.push_back(3);
	v2.push_back(4);
	v2.push_back(5);
	
	TEST_REAL_EQUAL(pearsonCorrelation(v1,v2),1)

	v2.clear();
	v2.push_back(-1);
	v2.push_back(-2);
	v2.push_back(-3);
	v2.push_back(-4);
	v2.push_back(-5);


	TEST_REAL_EQUAL(pearsonCorrelation(v1,v2),-1)

	
	v1.clear();
	v2.clear();
		
	v1.push_back(0.3716803);
	v1.push_back(0.2778111);
	v1.push_back(0.8152372);
	v1.push_back(0.7715097);
	v1.push_back(0.0163179);
	v1.push_back(-0.4898738);
	v1.push_back(-0.6060137);
	v1.push_back(-0.8882970);
	v1.push_back(0.2913591);
	v1.push_back(-0.3661791);
	v1.push_back(0.1320750);
	v1.push_back(0.2637229);
	v1.push_back(-0.7390226);
	v1.push_back(-0.0395929);
	v1.push_back(0.3387334);
	v1.push_back(0.8598541);
	v1.push_back(0.7388236);
	v1.push_back(-0.5928083);
	v1.push_back(0.9226006);
	v1.push_back(-0.3571427);
	
	v2.push_back(0.6396969);
	v2.push_back(0.7942405);
	v2.push_back(-0.6364473);
	v2.push_back(-0.6845633);
	v2.push_back(-0.6908862);
	v2.push_back(-0.5034169);
	v2.push_back(0.5745298);
	v2.push_back(-0.1247591);
	v2.push_back(-0.5129564);
	v2.push_back(0.0745857);
	v2.push_back(0.0733665);
	v2.push_back(-0.0118882);
	v2.push_back(0.1763471);
	v2.push_back(0.1027599);
	v2.push_back(-0.9737805);
	v2.push_back(0.8747677);
	v2.push_back(0.9479392);
	v2.push_back(0.0843604);
	v2.push_back(-0.3518961);
	v2.push_back(-0.3034039);

	TEST_REAL_EQUAL(pearsonCorrelation(v1,v2),0)
	
	v1.clear();
	v2.clear();
	
	v1.push_back(-0.1833341);
	v1.push_back(0.6564449);
	v1.push_back(0.8725039);
	v1.push_back(0.3610921);
	v1.push_back(0.7926144);
	v1.push_back(0.1833341);
	v1.push_back(-0.6564449);
	v1.push_back(-0.4141061);
	v1.push_back(-0.8725039);
	v1.push_back(0.8269985);
	v1.push_back(-0.5878715);
	v1.push_back(-0.2950443);
	v1.push_back(-0.3610921);
	v1.push_back(-0.8269985);
	v1.push_back(-0.0470327);
	v1.push_back(0.4141061);
	v1.push_back(0.0470327);
	v1.push_back(0.2950443);
	v1.push_back(-0.7926144);
	v1.push_back(0.5878715);
	
	v2.push_back(0.0336114);
	v2.push_back(0.4309199);
	v2.push_back(0.7612631);
	v2.push_back(0.1303875);
	v2.push_back(0.6282377);
	v2.push_back(0.0336114);
	v2.push_back(0.4309199);
	v2.push_back(0.1714839);
	v2.push_back(0.7612631);
	v2.push_back(0.6839264);
	v2.push_back(0.3455929);
	v2.push_back(0.0870511);
	v2.push_back(0.1303875);
	v2.push_back(0.6839264);
	v2.push_back(0.0022121);
	v2.push_back(0.1714839);
	v2.push_back(0.0022121);
	v2.push_back(0.0870511);
	v2.push_back(0.6282377);
	v2.push_back(0.3455929);

	TEST_REAL_EQUAL(pearsonCorrelation(v1,v2),0)

RESULT

/////////////////////////////////////////////////////////////);
/////////////////////////////////////////////////////////////
END_TEST
