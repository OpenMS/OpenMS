// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/CsvFile.h>
///////////////////////////
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GumbelDistributionFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GumbelDistributionFitter* ptr = nullptr;
GumbelDistributionFitter* nullPointer = nullptr;
START_SECTION(GumbelDistributionFitter())
	ptr = new GumbelDistributionFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~GumbelDistributionFitter()))
  delete ptr;
	NOT_TESTABLE
END_SECTION

START_SECTION((GumbelDistributionFitResult fit(std::vector<DPosition<2> >& points)))
	DPosition<2> pos;
  vector<DPosition<2> > points;

	pos.setX(-2.7); pos.setY(0.017); points.push_back(pos);
	pos.setX(-2.5); pos.setY(0.025); points.push_back(pos);
	pos.setX(-2); pos.setY(0.052); points.push_back(pos);
	pos.setX(-1); pos.setY(0.127); points.push_back(pos);
	pos.setX(-0.7); pos.setY(0.147); points.push_back(pos);
	pos.setX(-0.01); pos.setY(0.178); points.push_back(pos);
	pos.setX(0); pos.setY(0.178); points.push_back(pos);
	pos.setX(0.2); pos.setY(0.182); points.push_back(pos);
	pos.setX(0.5); pos.setY(0.184); points.push_back(pos);
	pos.setX(1); pos.setY(0.179); points.push_back(pos);
	pos.setX(1.3); pos.setY(0.171); points.push_back(pos);
	pos.setX(1.9); pos.setY(0.151); points.push_back(pos);
	pos.setX(2.5); pos.setY(0.127); points.push_back(pos);
	pos.setX(2.6); pos.setY(0.123); points.push_back(pos);
	pos.setX(2.7); pos.setY(0.119); points.push_back(pos);
	pos.setX(2.8); pos.setY(0.115); points.push_back(pos);
	pos.setX(2.9); pos.setY(0.111); points.push_back(pos);
	pos.setX(3); pos.setY(0.108); points.push_back(pos);
	pos.setX(3.5); pos.setY(0.089); points.push_back(pos);
	pos.setX(3.9); pos.setY(0.076); points.push_back(pos);
	pos.setX(4.01); pos.setY(0.073); points.push_back(pos);
	pos.setX(4.22); pos.setY(0.067); points.push_back(pos);
	pos.setX(4.7); pos.setY(0.054); points.push_back(pos);
	pos.setX(4.9); pos.setY(0.05); points.push_back(pos);
	pos.setX(5); pos.setY(0.047); points.push_back(pos);
	pos.setX(6); pos.setY(0.03); points.push_back(pos);
	pos.setX(7); pos.setY(0.017); points.push_back(pos);
	pos.setX(7.5); pos.setY(0.015); points.push_back(pos);
	pos.setX(7.9); pos.setY(0.012); points.push_back(pos);
	pos.setX(8.03); pos.setY(0.011); points.push_back(pos);
	//a= 0.5, b = 2
	

	ptr = new GumbelDistributionFitter;
	GumbelDistributionFitter::GumbelDistributionFitResult init_param;
	init_param.a = 1.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result = ptr->fit(points);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result.a, 0.5)
	TEST_REAL_SIMILAR(result.b, 2.0)	
	
	vector<DPosition<2> > points2;
	pos.setX(0); pos.setY(0.18); points2.push_back(pos);
	pos.setX(0.2); pos.setY(0.24); points2.push_back(pos);
	pos.setX(0.5); pos.setY(0.32); points2.push_back(pos);
	pos.setX(1); pos.setY(0.37); points2.push_back(pos);
	pos.setX(1.3); pos.setY(0.35); points2.push_back(pos);
	pos.setX(1.9); pos.setY(0.27); points2.push_back(pos);
	pos.setX(2.5); pos.setY(0.18); points2.push_back(pos);
	pos.setX(2.6); pos.setY(0.16); points2.push_back(pos);
	pos.setX(3); pos.setY(0.12); points2.push_back(pos);
	pos.setX(5); pos.setY(0.02); points2.push_back(pos);
	//a = 1, b = 1
	
	init_param.a = 3.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result2 = ptr->fit(points2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result2.a, 1.0)
	TEST_REAL_SIMILAR(result2.b, 1.0)	
	delete ptr;
END_SECTION

START_SECTION((void setInitialParameters(const GumbelDistributionFitResult& result)))
{
  GumbelDistributionFitter f1;
  GumbelDistributionFitter::GumbelDistributionFitResult result;
  f1.setInitialParameters(result);
	
	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((GumbelDistributionFitter(const GumbelDistributionFitter& rhs)))
  NOT_TESTABLE
END_SECTION

START_SECTION((GumbelDistributionFitter& operator = (const GumbelDistributionFitter& rhs)))
  NOT_TESTABLE
END_SECTION

GumbelDistributionFitter::GumbelDistributionFitResult* p = nullptr;

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult()))
  p =  new GumbelDistributionFitter::GumbelDistributionFitResult;
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_REAL_SIMILAR(p->a, 1.0)
  TEST_REAL_SIMILAR(p->b, 2.0)
END_SECTION

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult(const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))
  p-> a = 5.0;
  p->b = 4.0;
  GumbelDistributionFitter::GumbelDistributionFitResult obj(*p);
  TEST_REAL_SIMILAR(obj.a, 5.0)
  TEST_REAL_SIMILAR(obj.b, 4.0)
END_SECTION

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult& operator = (const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))
  p-> a = 3.0;
  p->b = 2.2;
  GumbelDistributionFitter::GumbelDistributionFitResult obj = *p;
  TEST_REAL_SIMILAR(obj.a, 3.0)
  TEST_REAL_SIMILAR(obj.b, 2.2)
  delete p;
END_SECTION

START_SECTION(MLE)
  vector<double> rand_score_vector;

  CsvFile gumbeldata (OPENMS_GET_TEST_DATA_PATH("Gumbel_1D.csv"));
  StringList gumbeldata_strings;
  gumbeldata.getRow(0, gumbeldata_strings);

  // Load mixture of Gumbel and Gaussian (1D) from provided csv
  for (StringList::const_iterator it = gumbeldata_strings.begin(); it != gumbeldata_strings.end(); ++it)
  {
    if(!it->empty())
    {
      rand_score_vector.push_back(it->toDouble());
    }
  }
  vector<double> w (rand_score_vector.size(),1.0);

  TEST_EQUAL(rand_score_vector.size(),1200)

  GumbelMaxLikelihoodFitter gmlf({4.0,2.0});

  auto res = gmlf.fitWeighted(rand_score_vector, w);
  TEST_REAL_SIMILAR(res.a, 2)
  TEST_REAL_SIMILAR(res.b, 0.6)
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



