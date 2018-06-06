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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/EuclideanSimilarity.h>
#include <cmath>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EuclideanSimilarity, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EuclideanSimilarity* ptr = nullptr;
EuclideanSimilarity* nullPointer = nullptr;
START_SECTION(EuclideanSimilarity())
{
	ptr = new EuclideanSimilarity();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~EuclideanSimilarity())
{
	delete ptr;
}
END_SECTION

START_SECTION((EuclideanSimilarity(const EuclideanSimilarity &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((EuclideanSimilarity& operator=(const EuclideanSimilarity &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((float operator()(const std::pair< float, float > &a, const std::pair< float, float > &b) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2.0f,2.0f),make_pair(4.f,4.f)), 1-sqrt(8.0));
			TEST_REAL_SIMILAR(es(make_pair(9.f,0.1f),make_pair(2.8f,2.f)), 1-sqrt(42.05));
			TEST_REAL_SIMILAR(es(make_pair(12.f,0.0f),make_pair(2.f,0.0f)), 1-sqrt(100.0));
			es.setScale(sqrt(233.28f));
}
END_SECTION

START_SECTION((float operator()(const std::pair< float, float > &c) const ))
{
			EuclideanSimilarity es;
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(9.0f,0.1f)), 1-0);
			TEST_REAL_SIMILAR(es(make_pair(2.8f,2.0f)), 1-0);
}
END_SECTION

START_SECTION((void setScale(float x)))
{
			EuclideanSimilarity es;
			es.setScale(10);
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(es(make_pair(2.0f,2.0f),make_pair(4.f,4.f)), 1-(sqrt(8.0)/10));
			TEST_REAL_SIMILAR(es(make_pair(9.0f,0.1f),make_pair(2.8f,2.f)), 1-(sqrt(42.05)/10));
			TEST_REAL_SIMILAR(es(make_pair(12.0f,0.0f),make_pair(2.f,0.0f)), 1-(sqrt(100.0)/10));
			es.setScale(sqrt(233.28f));
			TEST_REAL_SIMILAR(es(make_pair(0.1f,0.1f),make_pair(10.9f,10.9f)), 1-1);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



