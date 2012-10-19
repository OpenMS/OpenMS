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
// $Maintainer: George Rosenberger $
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
  static const double arrx1[] = { 1.1, 2.0,3.3,3.9,4.9,6.2  };
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = { 0.9, 1.9,3.0,3.7,5.2,6.1  };
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );

  int c1 = MRMRTNormalizer::outlier_candidate(x1,y1);
  TEST_EQUAL(c1,4);

  static const double arrx2[] = { 1,2,3,4,5,6  };
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = { 1,2,3,4,5,6};
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );

  int c2 = MRMRTNormalizer::outlier_candidate(x2,y2);
  TEST_EQUAL(c2,0);

}
END_SECTION

START_SECTION((static std::vector<std::pair<double, double> > rm_outliers(std::vector<std::pair<double, double> > & pairs, double rsq_limit, double coverage_limit)))
{
  static const double arrx1[] = { 1.1,2.0,3.3,3.9,4.9,6.2 };
  std::vector<double> x1 (arrx1, arrx1 + sizeof(arrx1) / sizeof(arrx1[0]) );
  static const double arry1[] = { 0.9,1.9,3.0,3.7,5.2,6.1 };
  std::vector<double> y1 (arry1, arry1 + sizeof(arry1) / sizeof(arry1[0]) );

  std::vector<std::pair<double, double> > input1;
  for (Size i = 0; i < x1.size(); i++)
  {
    input1.push_back(std::make_pair(x1[i], y1[i]));
  }

  std::vector<std::pair<double, double> > output1 = MRMRTNormalizer::rm_outliers(input1, 0.9, 0.5);
  TEST_EQUAL( output1.size() , input1.size() );

  static const double arrx2[] = { 1.1,2.0,3.3,3.9,4.9,6.2 };
  std::vector<double> x2 (arrx2, arrx2 + sizeof(arrx2) / sizeof(arrx2[0]) );
  static const double arry2[] = { 0.9,1.9,7.0,3.7,5.2,6.1 };
  std::vector<double> y2 (arry2, arry2 + sizeof(arry2) / sizeof(arry2[0]) );

  std::vector<std::pair<double, double> > input2;
  for (Size i = 0; i < x2.size(); i++)
  { 
    input2.push_back(std::make_pair(x2[i], y2[i]));
  }
  
  std::vector<std::pair<double, double> > output2 = MRMRTNormalizer::rm_outliers(input2, 0.9, 0.5);
  TEST_EQUAL( output2.size() , input2.size() - 1 );

  TEST_EQUAL( output2[0].first,  input2[0].first );
  TEST_EQUAL( output2[1].second, input2[1].second );

  TEST_EQUAL( output2[2].first,  input2[3].first );
  TEST_EQUAL( output2[3].second, input2[4].second );

  static const double arrx3[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,1,21,22,23,24,25,26,27,28,29,30 };
  std::vector<double> x3 (arrx3, arrx3 + sizeof(arrx3) / sizeof(arrx3[0]) );
  static const double arry3[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,22,23,24,25,26,27,28,29,30 };
  std::vector<double> y3 (arry3, arry3 + sizeof(arry3) / sizeof(arry3[0]) );

  std::vector<std::pair<double, double> > input3;
  for (Size i = 0; i < x3.size(); i++)
  { 
    input3.push_back(std::make_pair(x3[i], y3[i]));
  }

  std::vector<std::pair<double, double> > output3 = MRMRTNormalizer::rm_outliers(input3, 0.9, 0.2);
  TEST_EQUAL( output3.size() , input3.size() - 2 );

  TEST_EQUAL( output3[18].first,  input3[18].first );
  TEST_EQUAL( output3[19].second, input3[21].second );
}
END_SECTION

START_SECTION(( static double chauvenet_probability(std::vector< double > &residuals, int pos) ))
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
