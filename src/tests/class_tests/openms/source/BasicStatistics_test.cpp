// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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

	double dvector_zero_data[] =
	{
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	};

	size_t num_numbers = sizeof (dvector_data) / sizeof (*dvector_data);
	size_t num_zeros = 12;

} // namespace OpenMS

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( BasicStatistics, "$Id$" );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION((BasicStatistics()))
{

	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	TOLERANCE_ABSOLUTE(0.1);
	STATUS( stats );
	TEST_REAL_SIMILAR( stats.sum(), 15228.2 );
	TEST_REAL_SIMILAR( stats.mean(), 96.4639 );
	TEST_REAL_SIMILAR( stats.variance(), 3276.51 );

	float fvector_coord[195];

	for ( int i = 0; i < 195; fvector_coord[i] = 1000.f - i, ++i ) ;

	BasicStatistics < double > stats2;
	stats2.update( &*dvector_data, dvector_data + num_numbers, &*fvector_coord );

	STATUS( stats );
	TEST_REAL_SIMILAR( stats2.sum(), stats.sum() );
	TEST_REAL_SIMILAR( stats2.mean(), 1000. - stats.mean() );
	TEST_REAL_SIMILAR( stats2.variance(), 3276.51 );

}
END_SECTION
//-----------------------------------------------------------
START_SECTION((BasicStatistics(BasicStatistics const &arg)))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	TOLERANCE_ABSOLUTE(0.1);

	BasicStatistics<> const & stats_cref = stats;

	BasicStatistics<> stats_copy (stats_cref);

	TEST_REAL_SIMILAR( stats_copy.sum(), stats.sum() );
	TEST_REAL_SIMILAR( stats_copy.mean(), stats.mean() );
	TEST_REAL_SIMILAR( stats_copy.variance(), stats.variance() );
}
END_SECTION
//-----------------------------------------------------------
START_SECTION((BasicStatistics& operator=(BasicStatistics const &arg)))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TEST_EQUAL(num_numbers,195);
	TOLERANCE_ABSOLUTE(0.1);

	BasicStatistics<> const & stats_cref = stats;

	BasicStatistics<> stats_copy;
	stats_copy = stats_cref;

	TEST_REAL_SIMILAR( stats_copy.sum(), stats.sum() );
	TEST_REAL_SIMILAR( stats_copy.mean(), stats.mean() );
	TEST_REAL_SIMILAR( stats_copy.variance(), stats.variance() );

}
END_SECTION
//-----------------------------------------------------------
START_SECTION((void clear()))
{
	BasicStatistics <> stats;
	stats.update( &*dvector_data, dvector_data + num_numbers );

	TOLERANCE_ABSOLUTE(0.1);
	TEST_REAL_SIMILAR( stats.sum(), 15228.2 );
	TEST_REAL_SIMILAR( stats.mean(), 96.4639 );
	TEST_REAL_SIMILAR( stats.variance(), 3276.51 );
	stats.clear();
	TEST_REAL_SIMILAR( stats.sum(), 0. );
	TEST_REAL_SIMILAR( stats.mean(), 0. );
	TEST_REAL_SIMILAR( stats.variance(), 0. );
}
END_SECTION
//-----------------------------------------------------------
START_SECTION((template <typename ProbabilityIterator> void update(ProbabilityIterator probability_begin, ProbabilityIterator const probability_end)))
{
	BasicStatistics<> stats;
	TEST_REAL_SIMILAR( stats.sum(), 0. );
	TEST_REAL_SIMILAR( stats.mean(), 0. );
	TEST_REAL_SIMILAR( stats.variance(), 0. );
	stats.update( &*dvector_data, dvector_data + num_numbers );
	TOLERANCE_ABSOLUTE(0.1);
	TEST_REAL_SIMILAR( stats.sum(), 15228.2 );
	TEST_REAL_SIMILAR( stats.mean(), 96.4639 );
	TEST_REAL_SIMILAR( stats.variance(), 3276.51 );

	stats.update( &*dvector_zero_data, dvector_zero_data + num_zeros );
	TEST_REAL_SIMILAR( stats.sum(), 0 );
	TEST_REAL_SIMILAR( stats.mean(), 0 );
	TEST_REAL_SIMILAR( stats.variance(), 0 );
}
END_SECTION
//-----------------------------------------------------------
START_SECTION((template <typename ProbabilityIterator, typename CoordinateIterator> void update(ProbabilityIterator const probability_begin, ProbabilityIterator const probability_end, CoordinateIterator const coordinate_begin)))
{
	BasicStatistics<> stats;
	TEST_REAL_SIMILAR( stats.sum(), 0. );
	TEST_REAL_SIMILAR( stats.mean(), 0. );
	TEST_REAL_SIMILAR( stats.variance(), 0. );

	float fvector_coord[195];
	for ( int i = 0; i < 195; fvector_coord[i] = 1000.f - i, ++i ) ;

	stats.update( &*dvector_data, dvector_data + num_numbers, &*fvector_coord );

	TOLERANCE_ABSOLUTE(0.1);
	TEST_REAL_SIMILAR( stats.sum(), 15228.2 );
	TEST_REAL_SIMILAR( stats.mean(), 1000.-96.4639 );
	TEST_REAL_SIMILAR( stats.variance(), 3276.51 );
}
END_SECTION;
//-----------------------------------------------------------

BasicStatistics<double> bid;

START_SECTION((RealType mean() const))
{
	TEST_EQUAL(bid.mean(),0.);
	// continued below
}
END_SECTION

START_SECTION((void setMean(RealType const &mean)))
{
	TEST_EQUAL(bid.mean(),0.);
	bid.setMean(17.);
	TEST_EQUAL(bid.mean(),17.);
}
END_SECTION

START_SECTION((RealType variance() const))
{
	TEST_EQUAL(bid.variance(),0.);
	// continued below
}
END_SECTION

START_SECTION((void setVariance(RealType const &variance)))
{
	TEST_EQUAL(bid.variance(),0.);
	bid.setVariance(18.);
	TEST_EQUAL(bid.variance(),18.);
}
END_SECTION

START_SECTION((RealType sum() const))
{
	TEST_EQUAL(bid.sum(),0.);
	// continued below
}
END_SECTION

START_SECTION((void setSum(RealType const &sum)))
{
	TEST_EQUAL(bid.sum(),0.);
	bid.setSum(19.);
	TEST_EQUAL(bid.sum(),19.);
}
END_SECTION

//-----------------------------------------------------------

START_SECTION((static RealType sqrt2pi()))
{
	TEST_EQUAL(BasicStatistics<>::sqrt2pi(),2.50662827463100050240);
}
END_SECTION

START_SECTION((RealType normalDensity_sqrt2pi(RealType coordinate) const))
{
	bid.clear();
	bid.setMean(10.);
	bid.setVariance(3.);
	TOLERANCE_ABSOLUTE(.0001);
	TEST_REAL_SIMILAR(bid.normalDensity_sqrt2pi(10.),1.);
	TEST_REAL_SIMILAR(bid.normalDensity_sqrt2pi(7.),.22313016014842982893);
	TEST_REAL_SIMILAR(bid.normalDensity_sqrt2pi(9.),.84648172489061407405);
	TEST_REAL_SIMILAR(bid.normalDensity_sqrt2pi(11.),.84648172489061407405);
}
END_SECTION

START_SECTION((RealType normalDensity(RealType const coordinate) const))
{
	bid.clear();
	bid.setMean(10.);
	bid.setVariance(3.);
	TOLERANCE_ABSOLUTE(.0001);
	TEST_REAL_SIMILAR(bid.normalDensity(10.),1./2.50662827463100050240);
	TEST_REAL_SIMILAR(bid.normalDensity(7.),.22313016014842982893/2.50662827463100050240);
	TEST_REAL_SIMILAR(bid.normalDensity(9.),.84648172489061407405/2.50662827463100050240);
	TEST_REAL_SIMILAR(bid.normalDensity(11.),.84648172489061407405/2.50662827463100050240);
}
END_SECTION

//-----------------------------------------------------------
START_SECTION((void normalApproximation(probability_container &probability, typename probability_container::size_type const size)))
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
	TOLERANCE_ABSOLUTE(0.1);
	STATUS( stats );
	TEST_REAL_SIMILAR( stats.sum(), 6. );
	TEST_REAL_SIMILAR( stats.mean(), 13. / 6. );
	TEST_REAL_SIMILAR( stats.variance(), 17. / 36. );

	TOLERANCE_ABSOLUTE(0.000001);
	TEST_REAL_SIMILAR( stats.sqrt2pi(), 2.50662827463100050240 );
	TOLERANCE_ABSOLUTE(0.1);

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
		TEST_REAL_SIMILAR( probs[i], good_probs[i] );
	}

	// testing START_SECTION((void normalApproximation(probability_container &probability)))
	std::vector < double > probs2(6);
	stats.normalApproximation ( probs2 );

	for ( UInt i = 0; i < probs2.size(); ++i )
	{
		// STATUS( std::showpoint << "i:" << i << "  probs[i]:" << probs[i] << '\n');
		TEST_REAL_SIMILAR( probs2[i], good_probs[i] );
	}

}
END_SECTION
//-----------------------------------------------------------
START_SECTION((void normalApproximation(probability_container &probability)))
{
	// already tested in START_SECTION((void normalApproximation(probability_container &probability, typename probability_container::size_type const size)))
	NOT_TESTABLE;
}
END_SECTION
//-----------------------------------------------------------
START_SECTION((void normalApproximation(probability_container &probability, coordinate_container const &coordinate)))
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

	TOLERANCE_ABSOLUTE(0.1);
	TEST_REAL_SIMILAR( stats.sum(), stats2.sum() );
	TEST_REAL_SIMILAR( stats.mean(), stats2.mean() );
	TEST_REAL_SIMILAR( stats.variance(), stats2.variance() );

	std::vector < double > pos2;
	for ( double i = 0; i < fit.size(); pos2.push_back(i), i += 1./magic2 ) ;
	std::vector < double > fit2;
	stats.normalApproximation ( fit2, pos2 );
	stats2.update ( fit2.begin(), fit2.end() );
	STATUS( stats  );
	STATUS( stats2 );
	TEST_REAL_SIMILAR( stats.sum(), stats2.sum() );
	TEST_REAL_SIMILAR( stats.mean(), stats2.mean() / magic2 );
	TEST_REAL_SIMILAR( stats.variance(), stats2.variance() / magic2 / magic2 );


}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
