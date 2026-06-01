// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/CLUSTERING/EuclideanSimilarity.h>
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



