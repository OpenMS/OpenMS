// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: George Rosenberger, Hannes Roest $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(MRMRTNormalizer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// no constructor / destructor of static class

// MRMRTNormalizer() 
// ~MRMRTNormalizer() 
//


START_SECTION((static int outlier_candidate(std::vector<double> & x, std::vector<double> & y)))
{
}
END_SECTION

START_SECTION((static std::vector<std::pair<double, double> > rm_outliers(std::vector<std::pair<double, double> > & pairs, double rsq_limit, double coverage_limit)))
{
  static const double arr1[] = { 1,  3,4,2,10,12,4 };
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = { 100,3,4,2,11,100,7 };
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  std::vector<std::pair<double, double> > input;
  for (Size i = 0; i < data1.size(); i++)
  {
    input.push_back(std::make_pair(data1[i], data2[i]));
  }

  std::vector<std::pair<double, double> > output = MRMRTNormalizer::rm_outliers(input, 0.9, 0.5);
  TEST_EQUAL( output.size() , input.size() - 2 );

  TEST_EQUAL( output[0].first,  input[1].first );
  TEST_EQUAL( output[0].second, input[1].second );

  TEST_EQUAL( output[4].first,  input[6].first );
  TEST_EQUAL( output[4].second, input[6].second );

}
END_SECTION

START_SECTION((static bool chauvenet_probability(std::vector<double> & residuals, int pos)))
{

  static const double arr1[] = { 1,2,3,4,2,10,11,75,5,8,3,5,6,9,130 };
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 0), 0.61831553);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 1), 0.6387955);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 2), 0.65955473);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 3), 0.68057951);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 4), 0.6387955);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 5), 0.81146293);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 6), 0.8339146);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 7), 0.10161557);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 8), 0.70185552);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 9), 0.76703896);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 10), 0.65955473);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 11), 0.70185552);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 12), 0.72336784);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 13), 0.78916526);
  TEST_REAL_SIMILAR( MRMRTNormalizer::chauvenet_probability(data1, 14), 0.00126358);

}
END_SECTION

START_SECTION((static bool chauvenet(std::vector<double> & residuals, int pos)))
{

  static const double arr1[] = { 1,2,3,4,2,10,11,75,5,8,3,5,6,9,130 };
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 0), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 1), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 2), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 3), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 4), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 5), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 6), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 7), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 8), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 9), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 10), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 11), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 12), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 13), false);
  TEST_EQUAL( MRMRTNormalizer::chauvenet(data1, 14), true);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
