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

PairMatcher::FeatureMapType features;
features.resize(1);
features[0].setMZ(0.0f);
features[0].setRT(0.1f);

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
	p.setValue("rt_pair_dist",0.4);
	p.setValue("rt_stdev_low",0.1);
	p.setValue("rt_stdev_high",0.2);

	p.setValue("mz_pair_dist",5.0);
	p.setValue("mz_stdev",0.3);

  
  PairMatcher pm1(features);
	pm1.setParameters(p);
	
	PairMatcher::FeatureMapType empty_features;
  PairMatcher pm2(empty_features);
  pm2 = pm1;

  TEST_EQUAL(pm1==pm2,true);
RESULT

CHECK((PairMatcher(const PairMatcher& source)))
  PairMatcher pm1(features);

  PairMatcher pm2(pm1);

  TEST_EQUAL(pm1==pm2,true);
RESULT


features.resize(10);

//start
features[0].setRT(1.0f);
features[0].setMZ(1.0f);
features[0].setCharge(1);
features[0].setOverallQuality(1);
features[0].setIntensity(4.0);


//best
features[1].setRT(1.5f);
features[1].setMZ(5.0f);
features[1].setCharge(1);
features[1].setOverallQuality(1);
features[1].setIntensity(2.0);

//inside (down, up, left, right)
features[2].setRT(1.0f);
features[2].setMZ(5.0f);
features[2].setCharge(1);
features[2].setOverallQuality(1);

features[3].setRT(3.0f);
features[3].setMZ(5.0f);
features[3].setCharge(1);
features[3].setOverallQuality(1);

features[4].setRT(1.5f);
features[4].setMZ(4.8f);
features[4].setCharge(1);
features[4].setOverallQuality(1);

features[5].setRT(1.5f);
features[5].setMZ(5.2f);
features[5].setCharge(1);
features[5].setOverallQuality(1);

//outside (down, up, left, right)
features[6].setRT(0.0f);
features[6].setMZ(5.0f);
features[6].setCharge(1);
features[6].setOverallQuality(1);

features[7].setRT(4.0f);
features[7].setMZ(5.0f);
features[7].setCharge(1);
features[7].setOverallQuality(1);

features[8].setRT(1.5f);
features[8].setMZ(4.0f);
features[8].setCharge(1);
features[8].setOverallQuality(1);

features[9].setRT(1.5f);
features[9].setMZ(6.0f);
features[9].setCharge(1);
features[9].setOverallQuality(1);

PairMatcher pm(features);
Param p;
p.setValue("rt_pair_dist",0.4);
p.setValue("rt_stdev_low",0.5);
p.setValue("rt_stdev_high",1);
p.setValue("mz_pair_dist",4.0);
p.setValue("mz_stdev",0.3);
pm.setParameters(p);

CHECK((const PairVectorType& run()))	
	const PairMatcher::PairVectorType& pairs = pm.run();
	
	TEST_EQUAL(pairs.size(),5)
	ABORT_IF(pairs.size()!=5)
	
	PRECISION(0.01)
	
	TEST_REAL_EQUAL(pairs[0].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[0].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getMZ(),5.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[0].getQuality(),0.4237f);
	
	TEST_REAL_EQUAL(pairs[1].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[1].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[1].getSecond().getMZ(),4.8f);
	TEST_REAL_EQUAL(pairs[1].getSecond().getRT(),1.5f);
	TEST_REAL_EQUAL(pairs[1].getQuality(),0.4647f);
		
	TEST_REAL_EQUAL(pairs[2].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[2].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[2].getSecond().getMZ(),5.0f);
	TEST_REAL_EQUAL(pairs[2].getSecond().getRT(),1.5f);
	TEST_REAL_EQUAL(pairs[2].getQuality(),0.9203f);
		
	TEST_REAL_EQUAL(pairs[3].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[3].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[3].getSecond().getMZ(),5.2f);
	TEST_REAL_EQUAL(pairs[3].getSecond().getRT(),1.5f);
	TEST_REAL_EQUAL(pairs[3].getQuality(),0.4647f);
		
	TEST_REAL_EQUAL(pairs[4].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[4].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[4].getSecond().getMZ(),5.0f);
	TEST_REAL_EQUAL(pairs[4].getSecond().getRT(),3.0f);
	TEST_REAL_EQUAL(pairs[4].getQuality(),0.1095f);
RESULT

CHECK((const PairVectorType& getBestPairs()))
	const PairMatcher::PairVectorType& pairs = pm.getBestPairs();
	TEST_EQUAL(pairs.size(),1);
	ABORT_IF(pairs.size()!=1)
	TEST_REAL_EQUAL(pairs[0].getFirst().getMZ(),1.0f);
	TEST_REAL_EQUAL(pairs[0].getFirst().getRT(),1.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getMZ(),5.0f);
	TEST_REAL_EQUAL(pairs[0].getSecond().getRT(),1.5f);
	TEST_REAL_EQUAL(pairs[0].getQuality(),0.9203f);
RESULT

CHECK((static void printInfo(std::ostream& out, const PairVectorType& pairs)))
	const PairMatcher::PairVectorType& pairs = pm.getBestPairs();
	stringstream s;
	PairMatcher::printInfo(s,pairs);
	TEST_EQUAL(s.str(), "Found the following 1 pairs:\nQuality\tFirst[RT]\tFirst[MZ]\tFirst[Int]\tFirst[Corr]\tSecond[RT]\tSecond[MZ]\tSecond[Int]\tSecond[Corr]\tRatio\tCharge\tDiff[RT]\tDiff[MZ]\n0.92\t1.00\t1.00\t4.00\t1.00\t1.50\t5.00\t2.00\t1.00\t2.00\t1\t0.50\t4.00\n");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



