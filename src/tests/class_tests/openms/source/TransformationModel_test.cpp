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
// $Maintainer: $
// $Authors: Hendrik Weisser, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
// #include <OpenMS/test_config.h>

///////////////////////////

// #include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include "/home/user/code/OpenMS/include/TransformationModel.h"

///////////////////////////

START_TEST(TransformationModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModel* ptr = 0;
TransformationModel* nullPointer = 0;

TransformationModel::DataPoints data, empty;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 2.0));
data.push_back(make_pair(1.0, 4.0));

START_SECTION((TransformationModel()))
{
  ptr = new TransformationModel();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((TransformationModel(const DataPoints &, const Param &)))
{
  ptr = new TransformationModel(TransformationModel::DataPoints(), Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TransformationModel()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual double evaluate(double value) const))
{
  // null model (identity):
  TransformationModel tm;
  TEST_REAL_SIMILAR(tm.evaluate(-3.14159), -3.14159);
  TEST_REAL_SIMILAR(tm.evaluate(0.0), 0.0);
  TEST_REAL_SIMILAR(tm.evaluate(12345678.9), 12345678.9);
}
END_SECTION

START_SECTION((void getParameters(Param & params) const))
{
  TransformationModel tm;
  Param p = tm.getParameters();
  TEST_EQUAL(p.empty(), true)
}
END_SECTION

START_SECTION(([EXTRA] static void getDefaultParameters(Param & params)))
{
  Param param;
  param.setValue("some-value", 12.3);
  TransformationModel::getDefaultParameters(param);
  TEST_EQUAL(param.empty(), true)
}
END_SECTION

START_SECTION((bool checkValidWeight(const string& weight, const vector<string>& valid_weights) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "ln(x)";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "1/y";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidYWeights()), true);
  test = "1/x2";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "none";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), false);
  test = "x2";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), false);
}
END_SECTION

START_SECTION((double weightDatum(double& datum, const string& weight) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  test = "none";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  test = "ln(x)";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), std::log(10e-5));
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), std::log(2.0));
  test = "1/x";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 1/10e-5);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/x2";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 1/std::pow(10e-5,2));
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(std::pow(2.0,2)));
  test = "ln(y)";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), std::log(10e-8));
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), std::log(2.0));
  test = "1/y";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 1/10e-8);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/y2";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 1/std::pow(10e-8,2));
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(std::pow(2.0,2)));
}
END_SECTION

START_SECTION((virtual void weightData(DataPoints& data, const Param& params)))
{
  TransformationModel::DataPoints data1;
  TransformationModel::DataPoints test1;
  Param param;
  TransformationModel::getDefaultParameters(param);
  TransformationModel dw(data, param);

  param.setValue("x_weight", "ln(x)");
  param.setValue("y_weight", "");
  test1.clear();
  test1.push_back(make_pair(std::log(10e-5), 1.0));
  test1.push_back(make_pair(std::log(1.0), 2.0));
  test1.push_back(make_pair(std::log(2.0), 4.0));  
  data1.clear();
  data1.push_back(make_pair(0.0, 1.0));
  data1.push_back(make_pair(1.0, 2.0));
  data1.push_back(make_pair(2.0, 4.0));
  dw.weightData(data1,param);
  for (size_t i = 0; i < data1.size(); ++i)
  {
    TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
    TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
  }

  param.setValue("x_weight", "");
  param.setValue("y_weight", "ln(y)");
  test1.clear();
  test1.push_back(make_pair(0.0, std::log(1.0)));
  test1.push_back(make_pair(1.0, std::log(2.0)));
  test1.push_back(make_pair(2.0, std::log(4.0)));  
  data1.clear();
  data1.push_back(make_pair(0.0, 1.0));
  data1.push_back(make_pair(1.0, 2.0));
  data1.push_back(make_pair(2.0, 4.0));
  dw.weightData(data1,param);
  for (size_t i = 0; i < data1.size(); ++i)
  {
    TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
    TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
  }
}
END_SECTION

START_SECTION((double unWeightDatum(double& datum, const string& weight) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "";
  TEST_REAL_SIMILAR(dw.unWeightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 2.0);
  test = "none";
  TEST_REAL_SIMILAR(dw.unWeightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 2.0);
  test = "ln(x)";
  TEST_REAL_SIMILAR(dw.unWeightDatum(std::log(11.0e5),test), 10e5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::exp(2.0));
  test = "1/x";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::abs(9.0e-5),test), 10e-5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/x2";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::pow(9.0e-5,2),test), 10e-5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::sqrt(1/std::abs(2.0)));
  test = "ln(y)";
  TEST_REAL_SIMILAR(dw.unWeightDatum(std::log(11.0e8),test), 10e8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::abs(std::exp(2.0)));
  test = "1/y";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::abs(9.0e-8),test), 10e-8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/y2";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::pow(9.0e-8,2),test), 10e-8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::sqrt(1/std::abs(2.0)));
}
END_SECTION

START_SECTION((virtual void unWeightData(DataPoints& data, const Param& params)))
{

  TransformationModel::DataPoints data1;
  TransformationModel::DataPoints test1;
  Param param;
  TransformationModel::getDefaultParameters(param);
  TransformationModel dw(data, param);

  param.setValue("x_weight", "ln(x)");
  param.setValue("y_weight", "");
  test1.clear();
  test1.push_back(make_pair(std::exp(0.0), 1.0));
  test1.push_back(make_pair(std::exp(1.0), 2.0));
  test1.push_back(make_pair(std::exp(2.0), 4.0));  
  data1.clear();
  data1.push_back(make_pair(0.0, 1.0));
  data1.push_back(make_pair(1.0, 2.0));
  data1.push_back(make_pair(2.0, 4.0));
  dw.unWeightData(data1,param);
  for (size_t i = 0; i < data1.size(); ++i)
  {
    TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
    TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
  }

  param.setValue("x_weight", "");
  param.setValue("y_weight", "ln(y)");
  test1.clear();
  test1.push_back(make_pair(0.0, std::exp(1.0)));
  test1.push_back(make_pair(1.0, std::exp(2.0)));
  test1.push_back(make_pair(2.0, std::exp(4.0)));  
  data1.clear();
  data1.push_back(make_pair(0.0, 1.0));
  data1.push_back(make_pair(1.0, 2.0));
  data1.push_back(make_pair(2.0, 4.0));
  dw.unWeightData(data1,param);
  for (size_t i = 0; i < data1.size(); ++i)
  {
    TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
    TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
