// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

// This one is going to be tested.
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


///////////////////////////

// More headers

// #include <algorithm>
// #include <functional>
// #include <iostream>
// #include <iterator>
// #include <vector>
// #include <string>


///////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( StatisticFunctions, "$Id$" );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType meanSquareError( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::list<DoubleReal> numbers1(20, 1.5);
	std::list<DoubleReal> numbers2(20, 1.3);
	DoubleReal result = 0;

	TOLERANCE_ABSOLUTE(0.000001);
	result = Math::meanSquareError(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_SIMILAR(result, 0.04);
}
END_SECTION

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType classificationRate( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::vector<DoubleReal> numbers1(20, 1);
	std::vector<DoubleReal> numbers2(20, 1);
	DoubleReal result = 0;

	numbers1.resize(40, -1);
	numbers2.resize(40, -1);

	numbers1[2] = -1;
	numbers1[7] = -1;
	numbers1[11] = -1;
	numbers1[15] = -1;
	numbers1[17] = -1;
	numbers1[25] = 1;
	numbers1[27] = 1;
	numbers1[29] = 1;
	numbers1[31] = 1;
	numbers1[37] = 1;

	result = Math::classificationRate(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_SIMILAR(result, 0.75);
}
END_SECTION

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType pearsonCorrelationCoefficient( const IteratorType1 begin_a, const IteratorType1 end_a, const IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::vector<DoubleReal> numbers1(20, 1.5);
	std::vector<DoubleReal> numbers2(20, 1.3);
	DoubleReal result = 0;

	numbers1[0] = 0.1;
	numbers2[0] = 0.5;
	numbers1[1] = 0.2;
	numbers2[1] = 0.7;
	numbers1[2] = 0.01;
	numbers2[2] = 0.03;
	numbers1[3] = 1.7;
	numbers2[3] = 1.0;
	numbers1[4] = 3.2;
	numbers2[4] = 4.0;
	
	result = Math::pearsonCorrelationCoefficient(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_SIMILAR(result, 0.897811);
	
// ************ TEST for nan *****************	
	std::vector<Real> vv1,vv2;
	vv1.push_back(1);
	vv1.push_back(1);
	vv1.push_back(1);
	vv1.push_back(1);
	vv1.push_back(1);

	vv2.push_back(1);
	vv2.push_back(2);
	vv2.push_back(3);
	vv2.push_back(4);
	vv2.push_back(5);
	
	result = Math::pearsonCorrelationCoefficient(vv1.begin(), vv1.end(), vv2.begin(), vv2.end());
	if (isnan(result) ) result = -1.0;

	TEST_REAL_SIMILAR(result, -1.0);
// ************ TEST for nan *****************	

	std::vector<Real> v1,v2;
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

	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),1);

	v2.clear();
	v2.push_back(-1);
	v2.push_back(-2);
	v2.push_back(-3);
	v2.push_back(-4);
	v2.push_back(-5);


	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),-1);


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

	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);

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

	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);
}
END_SECTION

START_SECTION([EXTRA](static void computeRank(std::vector<DoubleReal>& w)))
{
  std::vector<DoubleReal> numbers1(10, 1.5);

  numbers1[0] = 1.4;
  numbers1[1] = 0.2;
  numbers1[2] = 0.01;
  numbers1[3] = 1.7;
  numbers1[4] = 3.2;
  numbers1[5] = 2.2;
  
  TEST_REAL_SIMILAR(numbers1[0], 1.4);
  TEST_REAL_SIMILAR(numbers1[5], 2.2);
  
  Math::computeRank(numbers1);
  
  TEST_REAL_SIMILAR(numbers1[0], 0);
  TEST_REAL_SIMILAR(numbers1[1], 1);
  TEST_REAL_SIMILAR(numbers1[2], 2);
  TEST_REAL_SIMILAR(numbers1[3], 3);
  TEST_REAL_SIMILAR(numbers1[4], 4);
  TEST_REAL_SIMILAR(numbers1[5], 5);
}
END_SECTION        
    
START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType rankCorrelationCoefficient( const IteratorType1 begin_a, const IteratorType1 end_a, const IteratorType2 begin_b, const IteratorType2 end_b )))
{
  std::vector<DoubleReal> numbers1(10, 1.5);
  std::vector<DoubleReal> numbers2(10, 1.3);
  DoubleReal result = 0;

  numbers1[0] = 0.4;
  numbers2[0] = 0.5;
  numbers1[1] = 0.2;
  numbers2[1] = 0.7;
  numbers1[2] = 0.01;
  numbers2[2] = 0.03;
  numbers1[3] = 1.7;
  numbers2[3] = 1.0;
  numbers1[4] = 3.2;
  numbers2[4] = 4.0;
  numbers1[5] = 2.2;
  numbers2[5] = 3.0;
  
  result = Math::rankCorrelationCoefficient(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
  TEST_REAL_SIMILAR(result, 0.953125);
}
END_SECTION    


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
