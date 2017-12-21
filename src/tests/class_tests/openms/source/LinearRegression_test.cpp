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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

///////////////////////////

START_TEST(LinearRegression<Iterator>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

LinearRegression* ptr;
LinearRegression* nullPointer = nullptr;
START_SECTION((LinearRegression()))
  ptr = new LinearRegression;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~LinearRegression()))
  delete ptr;
END_SECTION

// Create a test data set
vector<double> x_axis(10);
vector<double> y_axis(10);
vector<double> y_axis0(10);
vector<double> weight(10);
for (int i=0; i < 10; ++i)
{
  x_axis[i]=i;
  y_axis[i]=2*i+4;
  y_axis0[i]=2*i; // no intercept
  weight[i]=1+i;
}

LinearRegression lin_reg, lin_reg2;

START_SECTION((template < typename Iterator > void computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, bool compute_goodness = true)))
  lin_reg.computeRegression(0.95,x_axis.begin(),x_axis.end(),y_axis.begin());
  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),4.0)
  TEST_REAL_SIMILAR(lin_reg.getChiSquared(),0.0)

  lin_reg2.computeRegression(0.95,x_axis.begin(),x_axis.end(),y_axis0.begin());
  TEST_REAL_SIMILAR(lin_reg2.getSlope(), 2.0)
  TEST_REAL_SIMILAR(lin_reg2.getIntercept(), 0.0)
  TEST_REAL_SIMILAR(lin_reg2.getChiSquared(), 0.0)
END_SECTION

START_SECTION((template < typename Iterator > void computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, bool compute_goodness = true)))
  lin_reg.computeRegressionWeighted(0.95,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin(), false);
  TEST_REAL_SIMILAR(lin_reg.getSlope(), 2.0)
  TEST_REAL_SIMILAR(lin_reg.getIntercept(), 4.0)
  lin_reg2.computeRegressionWeighted(0.95,x_axis.begin(),x_axis.end(),y_axis0.begin(),weight.begin(), false);
  TEST_REAL_SIMILAR(lin_reg2.getSlope(), 2.0)
  TEST_REAL_SIMILAR(lin_reg2.getIntercept(), 0.0)


  lin_reg.computeRegressionWeighted(0.95,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin(), true); // to get meta stats (tested below)
END_SECTION

START_SECTION((double getChiSquared() const))
  TEST_REAL_SIMILAR(lin_reg.getChiSquared(),0)
END_SECTION

START_SECTION((double getIntercept() const))
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),4.0)
END_SECTION

START_SECTION((double getLower() const))
  TEST_REAL_SIMILAR(lin_reg.getLower(),-2.0)
END_SECTION

START_SECTION((double getUpper() const))
  TEST_REAL_SIMILAR(lin_reg.getUpper(),-2.0)
END_SECTION

START_SECTION((double getSlope() const))
  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
END_SECTION

START_SECTION((double getStandDevRes() const))
  TEST_REAL_SIMILAR(lin_reg.getStandDevRes(),0.0)
END_SECTION

START_SECTION((double getStandErrSlope() const))
  TEST_REAL_SIMILAR(lin_reg.getStandErrSlope(),0.0)
END_SECTION

START_SECTION((double getRSquared() const))
  TEST_REAL_SIMILAR(lin_reg.getRSquared(),1.0)
END_SECTION

START_SECTION((double getTValue() const))
  TEST_REAL_SIMILAR(lin_reg.getTValue(),2.306)
END_SECTION

START_SECTION((double getXIntercept() const))
  TEST_REAL_SIMILAR(lin_reg.getXIntercept(),-2.0)
END_SECTION

START_SECTION((double getRSD() const))
  TEST_REAL_SIMILAR(lin_reg.getRSD(),0.0)
END_SECTION

START_SECTION((double getMeanRes() const))
  TEST_REAL_SIMILAR(lin_reg.getMeanRes(),0.0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
