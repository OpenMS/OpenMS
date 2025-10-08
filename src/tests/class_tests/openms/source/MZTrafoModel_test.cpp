// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
///////////////////////////

#include <OpenMS/MATH/MathFunctions.h>

using namespace OpenMS;
using namespace std;

using namespace OpenMS;
using namespace std;

START_TEST(MZTrafoModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MZTrafoModel* ptr = nullptr;
MZTrafoModel* null_ptr = nullptr;
START_SECTION(MZTrafoModel())
  ptr = new MZTrafoModel();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(~MZTrafoModel())
  delete ptr;
END_SECTION

CalibrationData cd;
for (Size i = 0; i < 10; ++i)
{
  cd.insertCalibrationPoint(100.100 + i, 200.200 + i, 128.5 + i, 200.0 + i, 1, 66);
  cd.insertCalibrationPoint(120.100 + i + 0.5, 400.200 + i, 128.5 + i, 200.0 + i, 1, 77);
}

START_SECTION(MZTrafoModel(bool ppm_model))
  NOT_TESTABLE // see predict()
END_SECTION

START_SECTION(static const std::string names_of_modeltype[])
END_SECTION

START_SECTION(static MODELTYPE nameToEnum(const std::string& name))
//"linear", "linear_weighted", "quadratic", "quadratic_weighted", "size_of_modeltype"
  TEST_EQUAL(MZTrafoModel::nameToEnum("linear"), MZTrafoModel::LINEAR)
  TEST_EQUAL(MZTrafoModel::nameToEnum("linear_weighted"), MZTrafoModel::LINEAR_WEIGHTED)
  TEST_EQUAL(MZTrafoModel::nameToEnum("quadratic"), MZTrafoModel::QUADRATIC)
  TEST_EQUAL(MZTrafoModel::nameToEnum("quadratic_weighted"), MZTrafoModel::QUADRATIC_WEIGHTED)
  TEST_EQUAL(MZTrafoModel::nameToEnum("size_of_modeltype"), MZTrafoModel::SIZE_OF_MODELTYPE)
  TEST_EQUAL(MZTrafoModel::nameToEnum("something_different_______"), MZTrafoModel::SIZE_OF_MODELTYPE)
END_SECTION

START_SECTION(static const std::string& enumToName(MODELTYPE mt))
  TEST_EQUAL(MZTrafoModel::enumToName(MZTrafoModel::LINEAR), "linear")
  TEST_EQUAL(MZTrafoModel::enumToName(MZTrafoModel::LINEAR_WEIGHTED), "linear_weighted")
  TEST_EQUAL(MZTrafoModel::enumToName(MZTrafoModel::QUADRATIC), "quadratic")
  TEST_EQUAL(MZTrafoModel::enumToName(MZTrafoModel::QUADRATIC_WEIGHTED), "quadratic_weighted")
  TEST_EQUAL(MZTrafoModel::enumToName(MZTrafoModel::SIZE_OF_MODELTYPE), "size_of_modeltype")
END_SECTION 

START_SECTION(static void setRANSACParams(const Math::RANSACParam& p))
  Math::RANSACParam p(10, 1000, 2.0, 25, false);
  MZTrafoModel::setRANSACParams(p);
END_SECTION


START_SECTION(static void setCoefficientLimits(double offset, double scale, double power))
  MZTrafoModel m;
  m.setCoefficientLimits(30, 4, 2);
  m.setCoefficients(25, 3, 1);
  TEST_EQUAL(MZTrafoModel::isValidModel(m), true)
  TEST_EQUAL(m.isTrained(), true)
  m.setCoefficients(-25, -3, -1);
  TEST_EQUAL(MZTrafoModel::isValidModel(m), true)
  TEST_EQUAL(m.isTrained(), true)

  m.setCoefficients(33, 3, 1);
  TEST_EQUAL(MZTrafoModel::isValidModel(m), false)
  TEST_EQUAL(m.isTrained(), true)
  m.setCoefficients(25, 5, 1);
  TEST_EQUAL(MZTrafoModel::isValidModel(m), false)
  TEST_EQUAL(m.isTrained(), true)
  m.setCoefficients(25, 3, 3);
  TEST_EQUAL(MZTrafoModel::isValidModel(m), false)
  TEST_EQUAL(m.isTrained(), true)
END_SECTION

START_SECTION(static bool isValidModel(const MZTrafoModel& trafo))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool isTrained() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(double getRT() const)
  NOT_TESTABLE // tested below
END_SECTION

START_SECTION(double predict(double mz) const)
  MZTrafoModel m(true);
  m.setCoefficients(25, 0, 0);

  double mz_theo = 100.0;
  double mz_obs = mz_theo + Math::ppmToMass(25.0, mz_theo);
  TEST_REAL_SIMILAR(m.predict(mz_obs), mz_theo);

  MZTrafoModel m2(false);
  m2.setCoefficients(0.25, 0, 0);

  mz_theo = 100.0;
  mz_obs = mz_theo + 0.25;
  TEST_REAL_SIMILAR(m2.predict(mz_obs), mz_theo);
END_SECTION


START_SECTION(static Size findNearest(const std::vector<MZTrafoModel>& tms, double rt))
  std::vector<MZTrafoModel> tms;
  MZTrafoModel m;
  // unsorted RT
  m.train(cd, MZTrafoModel::LINEAR, false, 100.0, 104.0); // RT = 102
  tms.push_back(m);
  m.train(cd, MZTrafoModel::LINEAR, false, 110.0, 114.0); // RT = 112
  tms.push_back(m);
  m.train(cd, MZTrafoModel::LINEAR, false, 106.0, 108.0); // RT = 107
  tms.push_back(m);
  m.train(cd, MZTrafoModel::LINEAR, false, 126.0, 128.0); // RT = 127
  tms.push_back(m);
  std::sort(tms.begin(), tms.end(), MZTrafoModel::RTLess());
  TEST_REAL_SIMILAR(tms[0].getRT(), 102.0)
  TEST_REAL_SIMILAR(tms[1].getRT(), 107.0)
  TEST_REAL_SIMILAR(tms[2].getRT(), 112.0)
  TEST_REAL_SIMILAR(tms[3].getRT(), 127.0)
  
  TEST_EQUAL(MZTrafoModel::findNearest(tms, 0.0), 0)
  TEST_EQUAL(MZTrafoModel::findNearest(tms, 100.0), 0)
  TEST_EQUAL(MZTrafoModel::findNearest(tms, 105.0), 1)
  TEST_EQUAL(MZTrafoModel::findNearest(tms, 140.0), 3)
END_SECTION


START_SECTION(bool train(const CalibrationData& cd, MODELTYPE md, bool use_RANSAC, double rt_left = -std::numeric_limits<double>::max(), double rt_right = std::numeric_limits<double>::max()))

  MZTrafoModel m;
  m.train(cd, MZTrafoModel::LINEAR, false);
  std::cout << m.toString() << std::endl;
  TEST_REAL_SIMILAR(m.getRT(), 0.0);

END_SECTION

START_SECTION(bool train(std::vector<double> error_mz, std::vector<double> theo_mz, std::vector<double> weights, MODELTYPE md, bool use_RANSAC))

  MZTrafoModel m;
  std::vector<double> error_mz = ListUtils::create<double>("10,11,9,10,9,11");
  std::vector<double> theo_mz = ListUtils::create<double>("100,200,300,400,500,600");
  std::vector<double> weights;
  Math::RANSACParam p(3, 1000, 4.0, 1, false);
  MZTrafoModel::setRANSACParams(p);
  m.train(error_mz, theo_mz, weights, MZTrafoModel::LINEAR, true);
  std::cout << m.toString() << std::endl;
  TEST_REAL_SIMILAR(m.predict(300.0 + Math::ppmToMass(10.0, 300.0)), 300.0);

  double a,b,c;
  m.getCoefficients(a,b,c);
  TEST_REAL_SIMILAR(a, 10.0)
  TEST_REAL_SIMILAR(b, 0.0)
  TEST_REAL_SIMILAR(c, 0.0)

  MZTrafoModel m2;
  m2.setCoefficients(m);
  m2.getCoefficients(a,b,c);
  TEST_REAL_SIMILAR(a, 10.0)
  TEST_REAL_SIMILAR(b, 0.0)
  TEST_REAL_SIMILAR(c, 0.0)
  
  m2.setCoefficients(1.0, 2.0, 3.0);
  m2.getCoefficients(a,b,c);
  TEST_REAL_SIMILAR(a, 1.0)
  TEST_REAL_SIMILAR(b, 2.0)
  TEST_REAL_SIMILAR(c, 3.0)               
END_SECTION

START_SECTION(void getCoefficients(double& intercept, double& slope, double& power))
  MZTrafoModel m;
  double a,b,c;
  TEST_EXCEPTION(Exception::Precondition, m.getCoefficients(a,b,c))
  // more tests see above
END_SECTION

START_SECTION(void setCoefficients(const MZTrafoModel& rhs))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setCoefficients(double intercept, double slope, double power))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(String toString() const)
  NOT_TESTABLE // tested above
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
