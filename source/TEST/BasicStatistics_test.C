// -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

///////////////////////////

// More headers

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>


///////////////////////////

namespace OpenMS
{

	// extra stuff required for this test
	double dvector_data[] =
	{
		82.70033, 18.53697, 130.43985, 71.42455, 50.63099, 20.31581, 30.19521,
		36.79161, 135.08596, 84.68491, 124.30681, 71.33620, 126.07538, 73.61598,
		130.07241, 88.97545, 112.80919, 81.12736, 170.80468, 74.20200, 29.40524,
		44.20175, 124.63237, 84.51534, 165.35688, 79.33067, 68.44432, 18.62523,
		112.01351, 77.03597, 29.93905, 49.71414, 30.82335, 61.01894, 113.46661,
		78.16001, 162.25406, 89.78833, 158.70900, 74.51220, 73.57289, 124.63237,
		84.51534, 165.35688, 79.33067, 68.44432, 18.62523, 112.01351, 77.03597,
		29.93905, 49.71414, 30.82335, 61.01894, 113.46661, 78.16001, 162.25406,
		89.78833, 158.70900, 74.51220, 73.57289, 17.14514, 130.14515, 83.68410,
		29.89634, 47.08373, 76.58917, 29.00928, 57.22767, 22.04459, 108.34564,
		79.49656, 140.83229, 67.81030, 28.82848, 78.72329, 31.32767, 62.28604,
		29.48579, 76.01188, 142.99623, 71.69667, 140.45532, 78.81924, 57.99051,
		19.66125, 29.71268, 63.73135, 65.07940, 27.78494, 127.22279, 67.27982,
		29.50484, 142.99623, 71.69667, 140.45532, 78.81924, 57.99051, 19.66125,
		29.71268, 63.73135, 65.07940, 27.78494, 127.22279, 67.27982, 29.50484,
		142.99623, 71.69667, 140.45532, 78.81924, 57.99051, 19.66125, 29.71268,
		63.73135, 65.07940, 27.78494, 127.22279, 67.27982, 29.50484, 54.54108,
		30.53517, 86.44319, 67.76178, 18.95834, 123.73745, 77.66034, 30.29570,
		60.94120, 142.92731, 82.77405, 141.99247, 76.17666, 157.02459, 78.28177,
		96.25540, 19.82469, 27.72561, 53.91157, 29.91151, 60.05424, 61.35466,
		16.14011, 163.18400, 77.86948, 153.28102, 91.43451, 29.32177, 83.93723,
		111.66644, 80.25561, 129.31559, 90.71809, 107.97381, 75.83463, 147.61897,
		78.47707, 29.93856, 68.92398, 177.78189, 81.44311, 68.58626, 24.30645,
		132.16980, 79.22136, 28.12488, 78.71920, 151.88722, 83.39256, 29.69833,
		71.72692, 52.76207, 15.71214, 116.18279, 75.74875, 115.52147, 91.14405,
		127.02429, 95.27849, 67.42286, 20.34733, 102.67339, 93.84615, 128.95366,
		69.28015, 138.62953, 94.72963, 129.24376, 66.28535, 27.90273, 58.98529,
		29.84631, 47.59564, 118.73823, 77.77458, 72.75859, 18.41622
	};

	size_t num_numbers = sizeof (dvector_data) / sizeof (*dvector_data);

} // namespace OpenMS

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( BasicStatistics, "$Id$" );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
CHECK((BasicStatistics()))
{

	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	PRECISION(0.1);
	STATUS( stats );
	TEST_REAL_EQUAL( stats.sum(), 15228.2 );
	TEST_REAL_EQUAL( stats.mean(), 96.4639 );
	TEST_REAL_EQUAL( stats.variance(), 3276.51 );

	float fvector_coord[195];

	for ( int i = 0; i < 195; fvector_coord[i] = 1000 - i, ++i ) ;

	BasicStatistics < double > stats2;
	stats2.update( &*dvector_data, dvector_data + num_numbers, &*fvector_coord );

	STATUS( stats );
	TEST_REAL_EQUAL( stats2.sum(), stats.sum() );
	TEST_REAL_EQUAL( stats2.mean(), 1000. - stats.mean() );
	TEST_REAL_EQUAL( stats2.variance(), 3276.51 );

}
RESULT
//-----------------------------------------------------------
CHECK((BasicStatistics(BasicStatistics const &arg)))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	PRECISION(0.1);

	BasicStatistics<> const & stats_cref = stats;

	BasicStatistics<> stats_copy (stats_cref);

	TEST_REAL_EQUAL( stats_copy.sum(), stats.sum() );
	TEST_REAL_EQUAL( stats_copy.mean(), stats.mean() );
	TEST_REAL_EQUAL( stats_copy.variance(), stats.variance() );
}
RESULT
//-----------------------------------------------------------
CHECK((BasicStatistics& operator=(BasicStatistics const &arg)))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	PRECISION(0.1);

	BasicStatistics<> const & stats_cref = stats;

	BasicStatistics<> stats_copy;
	stats_copy = stats_cref;

	TEST_REAL_EQUAL( stats_copy.sum(), stats.sum() );
	TEST_REAL_EQUAL( stats_copy.mean(), stats.mean() );
	TEST_REAL_EQUAL( stats_copy.variance(), stats.variance() );

}
RESULT
//-----------------------------------------------------------
CHECK((void clear()))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	PRECISION(0.1);
	TEST_REAL_EQUAL( stats.sum(), 15228.2 );
	TEST_REAL_EQUAL( stats.mean(), 96.4639 );
	TEST_REAL_EQUAL( stats.variance(), 3276.51 );
	stats.clear();
	TEST_REAL_EQUAL( stats.sum(), 0. );
	TEST_REAL_EQUAL( stats.mean(), 0. );
	TEST_REAL_EQUAL( stats.variance(), 0. );
}
RESULT
//-----------------------------------------------------------
CHECK((template <typename ProbabilityIterator> void update(ProbabilityIterator probability_begin, ProbabilityIterator const probability_end)))
{
	BasicStatistics<> stats;
	TEST_REAL_EQUAL( stats.sum(), 0. );
	TEST_REAL_EQUAL( stats.mean(), 0. );
	TEST_REAL_EQUAL( stats.variance(), 0. );
	stats.update( &*dvector_data, dvector_data + num_numbers );
	PRECISION(0.1);
	TEST_REAL_EQUAL( stats.sum(), 15228.2 );
	TEST_REAL_EQUAL( stats.mean(), 96.4639 );
	TEST_REAL_EQUAL( stats.variance(), 3276.51 );
}
RESULT
//-----------------------------------------------------------
CHECK((template <typename ProbabilityIterator, typename CoordinateIterator> void update(ProbabilityIterator const probability_begin, ProbabilityIterator const probability_end, CoordinateIterator const coordinate_begin)))
{
	BasicStatistics<> stats;
	TEST_REAL_EQUAL( stats.sum(), 0. );
	TEST_REAL_EQUAL( stats.mean(), 0. );
	TEST_REAL_EQUAL( stats.variance(), 0. );

	float fvector_coord[195];
	for ( int i = 0; i < 195; fvector_coord[i] = 1000 - i, ++i ) ;

	stats.update( &*dvector_data, dvector_data + num_numbers, &*fvector_coord );

	PRECISION(0.1);
	TEST_REAL_EQUAL( stats.sum(), 15228.2 );
	TEST_REAL_EQUAL( stats.mean(), 1000.-96.4639 );
	TEST_REAL_EQUAL( stats.variance(), 3276.51 );
}
RESULT;
//-----------------------------------------------------------

BasicStatistics<double> bid;

CHECK((RealType mean() const))
{
	TEST_EQUAL(bid.mean(),0.);
	// continued below
}
RESULT

CHECK((void setMean(RealType const &mean)))
{
	TEST_EQUAL(bid.mean(),0.);
	bid.setMean(17.);
	TEST_EQUAL(bid.mean(),17.);
}
RESULT

CHECK((RealType variance() const))
{
	TEST_EQUAL(bid.variance(),0.);
	// continued below
}
RESULT

CHECK((void setVariance(RealType const &variance)))
{
	TEST_EQUAL(bid.variance(),0.);
	bid.setVariance(18.);
	TEST_EQUAL(bid.variance(),18.);
}
RESULT

CHECK((RealType sum() const))
{
	TEST_EQUAL(bid.sum(),0.);
	// continued below
}
RESULT

CHECK((void setSum(RealType const &sum)))
{
	TEST_EQUAL(bid.sum(),0.);
	bid.setSum(19.);
	TEST_EQUAL(bid.sum(),19.);
}
RESULT

//-----------------------------------------------------------

CHECK((static RealType sqrt2pi()))
{
	TEST_EQUAL(BasicStatistics<>::sqrt2pi(),2.50662827463100050240);
}
RESULT

CHECK((RealType normalDensity_sqrt2pi(RealType coordinate) const))
{
	bid.clear();
	bid.setMean(10.);
	bid.setVariance(3.);
	PRECISION(.0001);
	TEST_REAL_EQUAL(bid.normalDensity_sqrt2pi(10.),1.);
	TEST_REAL_EQUAL(bid.normalDensity_sqrt2pi(7.),.22313016014842982893);
	TEST_REAL_EQUAL(bid.normalDensity_sqrt2pi(9.),.84648172489061407405);
	TEST_REAL_EQUAL(bid.normalDensity_sqrt2pi(11.),.84648172489061407405);
}
RESULT

CHECK((RealType normalDensity(RealType const coordinate) const))
{
	bid.clear();
	bid.setMean(10.);
	bid.setVariance(3.);
	PRECISION(.0001);
	TEST_REAL_EQUAL(bid.normalDensity(10.),1./2.50662827463100050240);
	TEST_REAL_EQUAL(bid.normalDensity(7.),.22313016014842982893/2.50662827463100050240);
	TEST_REAL_EQUAL(bid.normalDensity(9.),.84648172489061407405/2.50662827463100050240);
	TEST_REAL_EQUAL(bid.normalDensity(11.),.84648172489061407405/2.50662827463100050240);
}
RESULT

//-----------------------------------------------------------
CHECK((void normalApproximation(probability_container &probability, typename probability_container::size_type const size)))
{

	double dvector2_data[] =
	{
		0,
		1,
		3,
		2,
		0
	};

	size_t num2_numbers = sizeof (dvector2_data) / sizeof (*dvector2_data);

	BasicStatistics <> stats;
	stats.update( &*dvector2_data, dvector2_data + num2_numbers );

	TEST_EQUAL(num2_numbers,5);
	PRECISION(0.1);
	STATUS( stats );
	TEST_REAL_EQUAL( stats.sum(), 6. );
	TEST_REAL_EQUAL( stats.mean(), 13. / 6. );
	TEST_REAL_EQUAL( stats.variance(), 17. / 36. );

	PRECISION(0.000001);
	TEST_REAL_EQUAL( stats.sqrt2pi(), 2.50662827463100050240 );
	PRECISION(0.1);

	for ( double pos = -2.; pos < 6.; pos += .25 )
	{
		STATUS( std::showpoint << "pos:" << pos << "  density:" << stats.normalDensity_sqrt2pi ( pos ) / stats.sqrt2pi() );
	}

	std::vector < double > probs;
	stats.normalApproximation ( probs, 6 );

	double good_probs [] =
	{
		0.0241689,
		0.824253,
		3.38207,
		1.66963,
		0.0991695,
		0.000708684
	};

	for ( UInt i = 0; i < probs.size(); ++i )
	{
		// STATUS( std::showpoint << "i:" << i << "  probs[i]:" << probs[i] << '\n');
		TEST_REAL_EQUAL( probs[i], good_probs[i] );
	}

	// testing CHECK((void normalApproximation(probability_container &probability)))
	std::vector < double > probs2(6);
	stats.normalApproximation ( probs2 );

	for ( UInt i = 0; i < probs2.size(); ++i )
	{
		// STATUS( std::showpoint << "i:" << i << "  probs[i]:" << probs[i] << '\n');
		TEST_REAL_EQUAL( probs2[i], good_probs[i] );
	}

}
RESULT
//-----------------------------------------------------------
CHECK((void normalApproximation(probability_container &probability)))
{
	// already tested in CHECK((void normalApproximation(probability_container &probability, typename probability_container::size_type const size)))
}
RESULT
//-----------------------------------------------------------
CHECK((void normalApproximation(probability_container &probability, coordinate_container const &coordinate)))
{

	double magic1 = 200, magic2 = 100;

	std::vector < double > data ( (UInt) magic1 );
	std::copy ( &*dvector_data, dvector_data + num_numbers, std::back_inserter ( data ) );

	BasicStatistics <> stats;
	stats.update( data.begin(), data.end() );
	std::vector < double > fit;
	stats.normalApproximation ( fit, UInt ( data.size() + magic1 ) );
	BasicStatistics <> stats2;
	stats2.update( fit.begin(), fit.end() );
	STATUS( stats );
	STATUS( stats2 );

	PRECISION(0.1);
	TEST_REAL_EQUAL( stats.sum(), stats2.sum() );
	TEST_REAL_EQUAL( stats.mean(), stats2.mean() );
	TEST_REAL_EQUAL( stats.variance(), stats2.variance() );

	std::vector < double > pos2;
	for ( double i = 0; i < fit.size(); pos2.push_back(i), i += 1./magic2 ) ;
	std::vector < double > fit2;
	stats.normalApproximation ( fit2, pos2 );
	stats2.update ( fit2.begin(), fit2.end() );
	STATUS( stats  );
	STATUS( stats2 );
	TEST_REAL_EQUAL( stats.sum(), stats2.sum() );
	TEST_REAL_EQUAL( stats.mean(), stats2.mean() / magic2 );
	TEST_REAL_EQUAL( stats.variance(), stats2.variance() / magic2 / magic2 );


}
RESULT

CHECK((template< typename IteratorType1, typename IteratorType2 > static RealType meanSquareError( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
{
	std::list<DoubleReal> numbers1(20, 1.5);
	std::list<DoubleReal> numbers2(20, 1.3);
	DoubleReal result = 0;

	PRECISION(0.000001);
	result = BasicStatistics<DoubleReal>::meanSquareError(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_EQUAL(result, 0.04);
}
RESULT

CHECK((template< typename IteratorType1, typename IteratorType2 > static RealType classificationRate( IteratorType1 begin_a, const IteratorType1 end_a, IteratorType2 begin_b, const IteratorType2 end_b )))
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

	result = BasicStatistics<DoubleReal>::classificationRate(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_EQUAL(result, 0.75);
}
RESULT

CHECK((template< typename IteratorType1, typename IteratorType2 > static RealType pearsonCorrelationCoefficient( const IteratorType1 begin_a, const IteratorType1 end_a, const IteratorType2 begin_b, const IteratorType2 end_b )))
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
	result = BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(numbers1.begin(), numbers1.end(), numbers2.begin(), numbers2.end());
	TEST_REAL_EQUAL(result, 0.897811);

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

	TEST_REAL_EQUAL(BasicStatistics<Real>::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),1);

	v2.clear();
	v2.push_back(-1);
	v2.push_back(-2);
	v2.push_back(-3);
	v2.push_back(-4);
	v2.push_back(-5);


	TEST_REAL_EQUAL(BasicStatistics<Real>::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),-1);


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

	TEST_REAL_EQUAL(BasicStatistics<Real>::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);

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

	TEST_REAL_EQUAL(BasicStatistics<Real>::pearsonCorrelationCoefficient(v1.begin(), v1.end(), v2.begin(), v2.end()),0);
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
