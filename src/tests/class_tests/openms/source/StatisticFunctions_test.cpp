// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Johannes Junker, Mathias Walzer, Chris Bielow $
// --------------------------------------------------------------------------

///////////////////////////
// This one is going to be tested.
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <boost/math/special_functions/fpclassify.hpp>
///////////////////////////

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <list>

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( StatisticFunctions, "$Id$" );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA](template <typename IteratorType> static double sum(IteratorType begin, IteratorType end)))
{
	int x[] = {-1, 0, 1, 2, 3};
	TEST_EQUAL(int(Math::sum(x, x + 5)), 5);
	TEST_EQUAL(int(Math::sum(x, x)), 0);

	DoubleList y;
	y.push_back(-1.0);
  y.push_back(-0.5);
  y.push_back(0.0);
  y.push_back(0.5);
  y.push_back(1.0);
  y.push_back(1.5);
  y.push_back(2.0);
	TEST_REAL_SIMILAR(Math::sum(y.begin(), y.end()), 3.5);
}
END_SECTION

START_SECTION([EXTRA](template <typename IteratorType> static double mean(IteratorType begin, IteratorType end)))
{
	int x[] = {-1, 0, 1, 2, 3};
	TEST_EQUAL(Math::mean(x, x + 5), 1);
	TEST_EXCEPTION(Exception::InvalidRange, Math::mean(x, x));

	DoubleList y = ListUtils::create<double>("-1.0,-0.5,0.0,0.5,1.0,1.5,2.0");
	TEST_REAL_SIMILAR(Math::mean(y.begin(), y.end()), 0.5);
}
END_SECTION

START_SECTION([EXTRA](template <typename IteratorType> static double median(IteratorType begin, IteratorType end, bool sorted = false)))
{
	int x[] = {-1, 0, 1, 2, 3};
	TEST_REAL_SIMILAR(Math::median(x, x + 5, true), 1.0);
  int x2[] = {-1, 0, 1, 2, 3, 4}; // (1+2)/2
  TEST_REAL_SIMILAR(Math::median(x2, x2 + 6, true), 1.5);
	TEST_EXCEPTION(Exception::InvalidRange, Math::median(x, x));

  // unsorted
	DoubleList y = ListUtils::create<double>("1.0,-0.5,2.0,0.5,-1.0,1.5,0.0");
	TEST_REAL_SIMILAR(Math::median(y.begin(), y.end()), 0.5);
	y.push_back(-1.5); // even length
	TEST_REAL_SIMILAR(Math::median(y.begin(), y.end()), 0.25);

  // sorted
  DoubleList z_odd = ListUtils::create<double>("-1.0,-0.5,0.0,0.5,1.0,1.5,2.0");
  TEST_REAL_SIMILAR(Math::median(z_odd.begin(), z_odd.end(), true), 0.5);
  DoubleList z_even = ListUtils::create<double>("-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0");
  TEST_REAL_SIMILAR(Math::median(z_even.begin(), z_even.end(), true), 0.25);
}
END_SECTION

START_SECTION([EXTRA](template <typename IteratorType> double MAD(IteratorType begin, IteratorType end, double median_of_numbers)))
{
  int x[] = {-1, 0, 1, 2, 3};
  TEST_EQUAL(Math::MAD(x, x + 5, 1), 1);   // median{2, 1, 0, 1, 2}
  int x2[] = {-1, 0, 1, 2, 3, 4}; // median = 1.5 --> median{2.5, 1.5, 0.5, 0.5, 1.5, 2.5}
  TEST_REAL_SIMILAR(Math::MAD(x2, x2 + 6, true), 1.5);
  
  DoubleList z_odd = ListUtils::create<double>("-1.0,-0.5,0.0,0.5,1.0,1.5,2.0"); // median{1.5, 1, 0.5, 0, 0.5, 1, 1.5} == median{0, 0.5, 0.5, 1, 1, 1.5 ,1.5}
  TEST_REAL_SIMILAR(Math::MAD(z_odd.begin(), z_odd.end(), 0.5), 1);
  DoubleList z_even = ListUtils::create<double>("-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0"); // median{2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5} == median{0, 0.5, 0.5, 1, 1, 1.5 , 1.5, 2}
  TEST_REAL_SIMILAR(Math::MAD(z_even.begin(), z_even.end(), 0.5), 1);
}
END_SECTION

START_SECTION([EXTRA](template <typename IteratorType> double MeanAbsoluteDeviation(IteratorType begin, IteratorType end, double mean_of_numbers)))
{
  int x1[] = {-1, 0, 1, 2, 3}; // mean = 1
  TEST_EQUAL(Math::MeanAbsoluteDeviation(x1, x1 + 5, 1), 1.2);
  
  int x2[] = {-1, 0, 1, 2, 3, 4}; // mean = 1.5
  TEST_REAL_SIMILAR(Math::MeanAbsoluteDeviation(x2, x2 + 6, 1.5), 1.5);
  
  // single element test
  int x3[] = {-1}; // mean = -1
  TEST_REAL_SIMILAR(Math::MeanAbsoluteDeviation(x3, x3 + 1, -1.0), 0.0);
  
  // empty range	
  int x4[] = {1}; // mean = 1
  TEST_EQUAL(std::isnan(Math::MeanAbsoluteDeviation(x4, x4, 1)), true); // NaN
  
  DoubleList z_odd = ListUtils::create<double>("-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0"); // mean = 0.5
  TEST_REAL_SIMILAR(Math::MeanAbsoluteDeviation(z_odd.begin(), z_odd.end(), 0.5), 0.857142);
  
  DoubleList z_even = ListUtils::create<double>("-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0"); // mean = 0.25
  TEST_REAL_SIMILAR(Math::MeanAbsoluteDeviation(z_even.begin(), z_even.end(), 0.25), 1);
}
END_SECTION

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType meanSquareError( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::list<double> numbers1(20, 1.5);
	std::list<double> numbers2(20, 1.3);
	double result = 0;

	TOLERANCE_ABSOLUTE(0.000001);
	result = Math::meanSquareError(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_SIMILAR(result, 0.04);
}
END_SECTION

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType classificationRate( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::vector<double> numbers1(20, 1);
	std::vector<double> numbers2(20, 1);
	double result = 0;

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
	std::vector<double> numbers1(20, 1.5);
	std::vector<double> numbers2(20, 1.3);
	double result = 0;

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
	std::vector<float> vv1,vv2;
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
	if (std::isnan(result) ) result = -1.0;

	TEST_REAL_SIMILAR(result, -1.0);
// ************ TEST for nan *****************

	std::vector<float> v1,v2;
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

	v1.push_back(0.3716803f);
	v1.push_back(0.2778111f);
	v1.push_back(0.8152372f);
	v1.push_back(0.7715097f);
	v1.push_back(0.0163179f);
	v1.push_back(-0.4898738f);
	v1.push_back(-0.6060137f);
	v1.push_back(-0.8882970f);
	v1.push_back(0.2913591f);
	v1.push_back(-0.3661791f);
	v1.push_back(0.1320750f);
	v1.push_back(0.2637229f);
	v1.push_back(-0.7390226f);
	v1.push_back(-0.0395929f);
	v1.push_back(0.3387334f);
	v1.push_back(0.8598541f);
	v1.push_back(0.7388236f);
	v1.push_back(-0.5928083f);
	v1.push_back(0.9226006f);
	v1.push_back(-0.3571427f);

	v2.push_back(0.6396969f);
	v2.push_back(0.7942405f);
	v2.push_back(-0.6364473f);
	v2.push_back(-0.6845633f);
	v2.push_back(-0.6908862f);
	v2.push_back(-0.5034169f);
	v2.push_back(0.5745298f);
	v2.push_back(-0.1247591f);
	v2.push_back(-0.5129564f);
	v2.push_back(0.0745857f);
	v2.push_back(0.0733665f);
	v2.push_back(-0.0118882f);
	v2.push_back(0.1763471f);
	v2.push_back(0.1027599f);
	v2.push_back(-0.9737805f);
	v2.push_back(0.8747677f);
	v2.push_back(0.9479392f);
	v2.push_back(0.0843604f);
	v2.push_back(-0.3518961f);
	v2.push_back(-0.3034039f);

	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);

	v1.clear();
	v2.clear();

	v1.push_back(-0.1833341f);
	v1.push_back(0.6564449f);
	v1.push_back(0.8725039f);
	v1.push_back(0.3610921f);
	v1.push_back(0.7926144f);
	v1.push_back(0.1833341f);
	v1.push_back(-0.6564449f);
	v1.push_back(-0.4141061f);
	v1.push_back(-0.8725039f);
	v1.push_back(0.8269985f);
	v1.push_back(-0.5878715f);
	v1.push_back(-0.2950443f);
	v1.push_back(-0.3610921f);
	v1.push_back(-0.8269985f);
	v1.push_back(-0.0470327f);
	v1.push_back(0.4141061f);
	v1.push_back(0.0470327f);
	v1.push_back(0.2950443f);
	v1.push_back(-0.7926144f);
	v1.push_back(0.5878715f);

	v2.push_back(0.0336114f);
	v2.push_back(0.4309199f);
	v2.push_back(0.7612631f);
	v2.push_back(0.1303875f);
	v2.push_back(0.6282377f);
	v2.push_back(0.0336114f);
	v2.push_back(0.4309199f);
	v2.push_back(0.1714839f);
	v2.push_back(0.7612631f);
	v2.push_back(0.6839264f);
	v2.push_back(0.3455929f);
	v2.push_back(0.0870511f);
	v2.push_back(0.1303875f);
	v2.push_back(0.6839264f);
	v2.push_back(0.0022121f);
	v2.push_back(0.1714839f);
	v2.push_back(0.0022121f);
	v2.push_back(0.0870511f);
	v2.push_back(0.6282377f);
	v2.push_back(0.3455929f);

	TEST_REAL_SIMILAR(Math::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);
}
END_SECTION

START_SECTION([EXTRA](static void computeRank(std::vector<double>& w)))
{
  std::vector<double> numbers1(10, 1.5);

  numbers1[0] = 1.4;
  numbers1[1] = 0.2;
  numbers1[2] = 0.01;
  numbers1[3] = 1.7;
  numbers1[4] = 3.2;
  numbers1[5] = 2.2;

  TEST_REAL_SIMILAR(numbers1[0], 1.4);
  TEST_REAL_SIMILAR(numbers1[5], 2.2);
  TEST_REAL_SIMILAR(numbers1[6], 1.5);
  TEST_REAL_SIMILAR(numbers1[9], 1.5);

  Math::computeRank(numbers1);

  TEST_REAL_SIMILAR(numbers1[0], 3);
  TEST_REAL_SIMILAR(numbers1[1], 2);
  TEST_REAL_SIMILAR(numbers1[2], 1);
  TEST_REAL_SIMILAR(numbers1[3], 8);
  TEST_REAL_SIMILAR(numbers1[4], 10);
  TEST_REAL_SIMILAR(numbers1[5], 9);
  TEST_REAL_SIMILAR(numbers1[6], 5.5);
  TEST_REAL_SIMILAR(numbers1[7], 5.5);
  TEST_REAL_SIMILAR(numbers1[8], 5.5);
  TEST_REAL_SIMILAR(numbers1[9], 5.5);
}
END_SECTION

START_SECTION([EXTRA](template< typename IteratorType1, typename IteratorType2 > static RealType rankCorrelationCoefficient( const IteratorType1 begin_a, const IteratorType1 end_a, const IteratorType2 begin_b, const IteratorType2 end_b )))
{
  std::vector<double> numbers1(10, 1.5);
  std::vector<double> numbers2(10, 1.3);
  std::vector<double> numbers3(10, 0.42);
  std::vector<double> numbers4(10, 0.0);
  double result = 0;

  for (Size i = 0; i < numbers4.size(); ++i)
  {
    numbers4[i] = (double)(i+1);
  }

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
  TEST_REAL_SIMILAR(result, 0.858064516129032);

  result = Math::rankCorrelationCoefficient(numbers1.begin(), numbers1.end(),
																						numbers2.rbegin(), numbers2.rend());
  TEST_REAL_SIMILAR(result, 0.303225806451613);

  result = Math::rankCorrelationCoefficient(numbers3.begin(), numbers3.end(), numbers4.begin(), numbers4.end());
  TEST_REAL_SIMILAR(result, 0.0);

  result = Math::rankCorrelationCoefficient(numbers3.begin(), numbers3.end(), numbers3.begin(), numbers3.end());
  TEST_REAL_SIMILAR(result, 0.0);

  result = Math::rankCorrelationCoefficient(numbers4.begin(), numbers4.end(), numbers4.begin(), numbers4.end());
  TEST_REAL_SIMILAR(result, 1.0);

  result = Math::rankCorrelationCoefficient(numbers4.begin(), numbers4.end(), numbers4.rbegin(), numbers4.rend());
  TEST_REAL_SIMILAR(result, -1.0);
}
END_SECTION

START_SECTION([EXTRA](template <typename IteratorType> static double quantile(IteratorType begin, IteratorType end, UInt quantile, bool sorted = false) ))
{
  std::vector<int> x = {3,6,7,8,8,10,13,15,16,20};
  std::vector<int> y = {3,6,7,8,8,10,13,15,16};

  TEST_REAL_SIMILAR(Math::quantile1st(x.begin(), x.end(), true), 6.5);
  TEST_REAL_SIMILAR(Math::median(x.begin(), x.end(), true), 9.0);
  TEST_REAL_SIMILAR(Math::quantile3rd(x.begin(), x.end(), true), 15.5);
  TEST_REAL_SIMILAR(Math::quantile1st(y.begin(), y.end(), true),6.5);
  TEST_REAL_SIMILAR(Math::median(y.begin(), y.end(), true), 8.0);
  TEST_REAL_SIMILAR(Math::quantile3rd(y.begin(), y.end(), true), 14.0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
