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
// $Authors: Christian Ehrlich, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace Math;

START_TEST(QuadraticRegression, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QuadraticRegression* ptr = nullptr;
QuadraticRegression* null_ptr = nullptr;
START_SECTION(QuadraticRegression())
{
	ptr = new QuadraticRegression();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuadraticRegression())
{
	delete ptr;
}
END_SECTION

  // Create a test data set
vector<double> x_axis(10);
vector<double> y_axis(10);
vector<double> y_axis0(10);
vector<double> weight(10);
for (int i=0; i < 10; ++i)
{
  x_axis[i] = i;
  y_axis[i] = 5.5*i*i + 2*i + 4;
  y_axis0[i] = 5.5*i*i + 2*i; // no intercept
  weight[i]=1+i;
}

QuadraticRegression q_reg, q_reg2;

START_SECTION((template < typename Iterator > void computeRegression(Iterator x_begin, Iterator x_end, Iterator y_begin)))
{
  q_reg.computeRegression(x_axis.begin(), x_axis.end(), y_axis.begin());
  TEST_REAL_SIMILAR(q_reg.getA(), 4.0)
  TEST_REAL_SIMILAR(q_reg.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg.getChiSquared(), 0.0)

  q_reg2.computeRegression(x_axis.begin(), x_axis.end(), y_axis0.begin());
  TEST_REAL_SIMILAR(q_reg2.getA(), 0.0)
  TEST_REAL_SIMILAR(q_reg2.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg2.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg2.getChiSquared(), 0.0)
}
END_SECTION

START_SECTION((template < typename Iterator > void computeRegressionWeighted(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
{
  q_reg.computeRegressionWeighted(x_axis.begin(), x_axis.end(), y_axis.begin(), weight.begin());
  TEST_REAL_SIMILAR(q_reg.getA(), 4.0)
  TEST_REAL_SIMILAR(q_reg.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg.getChiSquared(), 0.0)

  q_reg2.computeRegressionWeighted(x_axis.begin(), x_axis.end(), y_axis0.begin(), weight.begin());
  TEST_REAL_SIMILAR(q_reg2.getA(), 0.0)
  TEST_REAL_SIMILAR(q_reg2.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg2.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg2.getChiSquared(), 0.0)
}
END_SECTION

START_SECTION((double eval(double x) const ))
{
  double x = 100.0;
  TEST_REAL_SIMILAR(q_reg.eval(x), x*x*5.5 + x*2 + 4)
}
END_SECTION

START_SECTION(static double eval(double A, double B, double C, double x))
{
  double x = 100.0;
  TEST_REAL_SIMILAR(QuadraticRegression::eval(4.0, 2.0, 5.5, x), x*x*5.5 + x*2 + 4)
}
END_SECTION

START_SECTION((double getA() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getB() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getC() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getChiSquared() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



