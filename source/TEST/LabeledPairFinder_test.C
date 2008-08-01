// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/KERNEL/FeatureMap.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LabeledPairFinder, "$Id LabeledPairFinder_test.C 139 2006-07-14 10:08:39Z jjoachim $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LabeledPairFinder* ptr = 0;
CHECK((LabeledPairFinder()))
	ptr = new LabeledPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~LabeledPairFinder()))
	delete ptr;
RESULT

FeatureMap<> features;
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

CHECK((virtual void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map)))
	LabeledPairFinder pm;
	Param p;
	p.setValue("rt_pair_dist",0.4);
	p.setValue("rt_dev_low",1.0);
	p.setValue("rt_dev_high",2.0);
	p.setValue("mz_pair_dist",4.0);
	p.setValue("mz_dev",0.6);
	pm.setParameters(p);

	ConsensusMap output;
	TEST_EXCEPTION(Exception::IllegalArgument,pm.run(vector<ConsensusMap>(),output));
	vector<ConsensusMap> input(1);
	ConsensusMap::convert(5,features,input[0]);
	output.getFileDescriptions()[5].label = "light";
	output.getFileDescriptions()[5].filename = "filename";
	output.getFileDescriptions()[8] = output.getFileDescriptions()[5];
	output.getFileDescriptions()[8].label = "heavy";
	
	pm.run(input,output);

	TEST_EQUAL(output.size(),1);
	ABORT_IF(output.size()!=1)
	TEST_REAL_EQUAL(output[0].begin()->getMZ(),1.0f);
	TEST_REAL_EQUAL(output[0].begin()->getRT(),1.0f);
	TEST_REAL_EQUAL(output[0].rbegin()->getMZ(),5.0f);
	TEST_REAL_EQUAL(output[0].rbegin()->getRT(),1.5f);
	TEST_REAL_EQUAL(output[0].getQuality(),0.959346);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



