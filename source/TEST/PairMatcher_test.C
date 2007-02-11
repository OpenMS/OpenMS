// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PairMatcher, "$Id PairMatcher_test.C 139 2006-07-14 10:08:39Z jjoachim $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

enum DimensionId
{
	RT = DimensionDescription < LCMS_Tag >::RT,
	MZ = DimensionDescription < LCMS_Tag >::MZ
};

typedef PairMatcher::FeatureMapType Features;
typedef PairMatcher::PairVectorType Pairs;

Features features;
features.resize(4);
features[0].getPosition()[MZ] = 0.0f;
features[0].getPosition()[RT] = 0.1f;
features[1].getPosition()[MZ] = 4.0f;
features[1].getPosition()[RT] = 0.0f;
features[2].getPosition()[MZ] = 5.0f;
features[2].getPosition()[RT] = 0.9f;
features[3].getPosition()[MZ] = 4.0f;
features[3].getPosition()[RT] = 0.2f;

features[0].getIntensity() = 1.0f;
features[1].getIntensity() = 2.0f;
features[2].getIntensity() = 3.0f;
features[3].getIntensity() = 4.0f;

features[0].setCharge(1);
features[1].setCharge(1);
features[2].setCharge(1);
features[3].setCharge(1);

features[0].setOverallQuality(1.0f);
features[1].setOverallQuality(1.0f);
features[2].setOverallQuality(1.0f);
features[3].setOverallQuality(1.0f);

PairMatcher* ptr = 0;
CHECK((PairMatcher(FeatureMapType& features)))
	ptr = new PairMatcher(features);
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(PairMatcher::getProductName(), "PairMatcher")
	TEST_EQUAL(PairMatcher(features).getName(), "PairMatcher")
RESULT

CHECK((~PairMatcher()))
	delete ptr;
RESULT

CHECK((PairMatcher& operator = (const PairMatcher& source)))
	Param p;	
	p.setValue("rt_stdev_low",0.1);
	p.setValue("rt_stdev_high",0.2);
	p.setValue("mz_stdev",0.3);
	p.setValue("mz_pair_dist",5.0);
	p.setValue("rt_pair_dist",0.4);
  
  PairMatcher pm1(features);
	pm1.setParameters(p);
	
	Features empty_features;
  PairMatcher pm2(empty_features);
  pm2 = pm1;

  TEST_EQUAL(pm1==pm2,true);
RESULT

CHECK((PairMatcher(const PairMatcher& source)))
  PairMatcher pm1(features);

  PairMatcher pm2(pm1);

  TEST_EQUAL(pm1==pm2,true);

RESULT


CHECK((const PairVectorType& run()))
	PairMatcher pm(features);
	const Pairs& pairs = pm.run();
	TEST_EQUAL(pairs.size(),2)
	ABORT_IF(pairs.size()!=2)

	TEST_REAL_EQUAL(pairs[0].getFirst().getPosition()[MZ],0.0f);
	TEST_REAL_EQUAL(pairs[0].getFirst().getPosition()[RT],0.1f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getPosition()[MZ],4.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getPosition()[RT],0.2f);
	TEST_REAL_EQUAL(pairs[1].getFirst().getPosition()[MZ],0.0f);
	TEST_REAL_EQUAL(pairs[1].getFirst().getPosition()[RT],0.1f);
	TEST_REAL_EQUAL(pairs[1].getSecond().getPosition()[MZ],4.0f);
	TEST_REAL_EQUAL(pairs[1].getSecond().getPosition()[RT],0.0f);
RESULT


CHECK((const PairVectorType& getBestPairs()))
	PairMatcher pm(features);
	pm.run();
	const Pairs& pairs = pm.getBestPairs();
	TEST_EQUAL(pairs.size(),1);
	ABORT_IF(pairs.size()!=1)
	TEST_REAL_EQUAL(pairs[0].getFirst().getPosition()[MZ],0.0f);
	TEST_REAL_EQUAL(pairs[0].getFirst().getPosition()[RT],0.1f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getPosition()[MZ],4.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getPosition()[RT],0.0f);
RESULT

CHECK((static void printInfo(std::ostream& out, const PairVectorType& pairs)))
	PairMatcher pm(features);
	pm.run();
	const Pairs& pairs = pm.getBestPairs();
	stringstream s;
	PairMatcher::printInfo(s,pairs);
	TEST_EQUAL(s.str(), "Found the following 1 pairs:\nQuality\tFirst[RT]\tFirst[MZ]\tFirst[Int]\tFirst[Corr]\tSecond[RT]\tSecond[MZ]\tSecond[Int]\tSecond[Corr]\tRatio\tCharge\tDiff[RT]\tDiff[MZ]\n0.36\t0.10\t0.00\t1.00\t1.00\t0.00\t4.00\t2.00\t1.00\t0.50\t1\t0.10\t-4.00\n");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



