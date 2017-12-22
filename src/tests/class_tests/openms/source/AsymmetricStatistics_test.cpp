// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

///////////////////////////

// This one is going to be tested.
//
#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>

///////////////////////////

// More headers

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>


///////////////////////////

namespace OpenMS
{

  // test data
  double dvector_data[] =
    { 142.99623, 71.69667, 140.45532, 78.81924, 57.99051, 19.66125, 29.71268,
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
using namespace std;

/////////////////////////////////////////////////////////////

START_TEST(AsymmetricStatistics, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AsymmetricStatistics<double>* ptr = nullptr;
AsymmetricStatistics<double>* nullPointer = nullptr;
START_SECTION(AsymmetricStatistics())
{
	ptr = new AsymmetricStatistics<double>();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AsymmetricStatistics())
{
	delete ptr;
}
END_SECTION

START_SECTION((RealType variance1() const))
{
  // summy subtest
	TEST_EQUAL(0, 0)
}
END_SECTION

START_SECTION((RealType variance2() const))
{
   // summy subtest
	TEST_EQUAL(0, 0)
}
END_SECTION

START_SECTION((template <typename ProbabilityIterator, typename CoordinateIterator> void update(ProbabilityIterator const probability_begin, ProbabilityIterator const probability_end, CoordinateIterator const coordinate_begin)))
{

	// set the beginning of coordinates
	float fvector_coord[90];
	for ( int i = 0; i < 90; fvector_coord[i] = 1000.f - i, ++i ) ;

	// set basic statistics
	BasicStatistics < double > stats2;
	stats2.update( &*dvector_data, dvector_data + num_numbers, &*fvector_coord );

	TEST_EQUAL(num_numbers,90);
	TOLERANCE_ABSOLUTE(0.1);
	STATUS( stats2 );

	TEST_REAL_SIMILAR( stats2.sum(), 7096.78 );
	TEST_REAL_SIMILAR( stats2.mean(), 954.86 );
	TEST_REAL_SIMILAR( stats2.variance(), 638.663 );

	AsymmetricStatistics < double > asy;

	// test default values for variance1 and variance2
	TEST_REAL_SIMILAR( asy.variance1(), 0 );
	TEST_REAL_SIMILAR( asy.variance2(), 0 );

	// compute variance1 and variance2
	asy.update(&*dvector_data, dvector_data + num_numbers, &*fvector_coord);

	// test basic statistics
	TEST_REAL_SIMILAR( asy.sum(), 7096.78 );
	TEST_REAL_SIMILAR( asy.mean(), 954.86 );
	TEST_REAL_SIMILAR( asy.variance(), 638.663 );

	// test advanced statistics, computed in method update
	//
 	// Note: Marcel had some other numbers here,
	// but Clemens changed the algorithm since then.
	// Not clear what's right here, but anyway we could detect way-off errors.
	TEST_REAL_SIMILAR( asy.variance1(), 612.229 );
	TEST_REAL_SIMILAR( asy.variance2() , 665.783 );

}
END_SECTION


// The following test might explain and check a bit more thoroughly how the asy stats are computed.
START_SECTION([EXTRA](template <typename ProbabilityIterator, typename CoordinateIterator> void update(ProbabilityIterator const probability_begin, ProbabilityIterator const probability_end, CoordinateIterator const coordinate_begin)))
{
	AsymmetricStatistics < double > asy;

	double vector_coord[] = { 0, 1, 2, 3, 4,    5, 6, 7, 8, 9 };

	{
		double vector_data[]  = { 0, 0, 0, 2, 997,  0, 1, 0, 0, 0 };
		TEST_EQUAL( sizeof (vector_data) / sizeof (*vector_data),  sizeof (vector_coord) / sizeof (*vector_coord));
		UInt num_numbers = sizeof (vector_data) / sizeof (*vector_data);
		asy.update(&*vector_data, vector_data + num_numbers, &*vector_coord);

		// test basic statistics
		TEST_REAL_SIMILAR( asy.sum(), 1000 );

		TOLERANCE_ABSOLUTE(1E-10);
		TEST_REAL_SIMILAR( asy.mean(), 4 );
		TEST_REAL_SIMILAR( asy.variance(), .006 );

		// test advanced statistics, computed in method update
		TEST_REAL_SIMILAR( asy.variance1(), 2.*2./1001. );
		TEST_REAL_SIMILAR( asy.variance2(), 2.*4./999. );
	}

	{
		double vector_data[]  = { 0, 0, 0, 5, 994,  0, 0, 0, 0, 1 };
		TEST_EQUAL( sizeof (vector_data) / sizeof (*vector_data),  sizeof (vector_coord) / sizeof (*vector_coord));
		UInt num_numbers = sizeof (vector_data) / sizeof (*vector_data);
		asy.update(&*vector_data, vector_data + num_numbers, &*vector_coord);

		// test basic statistics
		TEST_REAL_SIMILAR( asy.sum(), 1000 );

		TOLERANCE_ABSOLUTE(1E-10);
		TEST_REAL_SIMILAR( asy.mean(), 4 );
		TEST_REAL_SIMILAR( asy.variance(), .030 );

		// test advanced statistics, computed in method update
		TEST_REAL_SIMILAR( asy.variance1(), (5.*1.)/(994./2.+5.) );
		TEST_REAL_SIMILAR( asy.variance2(), (25.*1.)/(994./2.+1.) );
	}

}
END_SECTION

//-----------------------------------------------------------

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
