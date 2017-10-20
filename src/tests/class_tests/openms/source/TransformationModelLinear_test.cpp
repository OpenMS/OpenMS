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
// $Authors: Hendrik Weisser, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>

///////////////////////////

START_TEST(TransformationModelLinear, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModelLinear* ptr = 0;
TransformationModelLinear* nullPointer = 0;

TransformationModel::DataPoints data, empty;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 2.0));
data.push_back(make_pair(1.0, 4.0));

START_SECTION((TransformationModelLinear(const DataPoints &, const Param &)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelLinear lm(empty, Param())); // need data
  ptr = new TransformationModelLinear(data, Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TransformationModelLinear()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual double evaluate(double value) const))
{
  ptr = new TransformationModelLinear(data, Param());

  TEST_REAL_SIMILAR(ptr->evaluate(-0.5), 0.0);
  TEST_REAL_SIMILAR(ptr->evaluate(0.0), 1.0);
  TEST_REAL_SIMILAR(ptr->evaluate(0.5), 2.0);
  TEST_REAL_SIMILAR(ptr->evaluate(1.0), 3.0);
  TEST_REAL_SIMILAR(ptr->evaluate(1.5), 4.0);

  delete ptr;
  data.push_back(make_pair(2.0, 2.0));
}
END_SECTION

START_SECTION((void getParameters(Param & params) const))
{
  Param p_in;
  p_in.setValue("symmetric_regression", "true");
  TransformationModelLinear lm(data, p_in);
  Param p_out = p_in;
  p_out.setValue("slope", 0.5);
  p_out.setValue("intercept", 1.75);

  TEST_EQUAL(lm.getParameters(), p_out);
  p_in.clear();
  p_in.setValue("slope", 12.3);
  p_in.setValue("intercept", -45.6);
  TransformationModelLinear lm2(empty, p_in);
  TEST_EQUAL(lm2.getParameters(), p_in);
}
END_SECTION

START_SECTION(([EXTRA] void getParameters(double&, double&)))
{
  Param param;
  param.setValue("slope", 12.3);
  param.setValue("intercept", -45.6);
  TransformationModelLinear lm(empty, param);
  double slope, intercept;
  lm.getParameters(slope, intercept);
  TEST_REAL_SIMILAR(param.getValue("slope"), slope);
  TEST_REAL_SIMILAR(param.getValue("intercept"), intercept);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
