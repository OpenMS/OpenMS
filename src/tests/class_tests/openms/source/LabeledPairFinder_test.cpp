// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
///////////////////////////

#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(LabeledPairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LabeledPairFinder* ptr = nullptr;
LabeledPairFinder* nullPointer = nullptr;
START_SECTION((LabeledPairFinder()))
	ptr = new LabeledPairFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~LabeledPairFinder()))
	delete ptr;
END_SECTION

FeatureMap features;
features.resize(10);
//start
features[0].setRT(1.0f);
features[0].setMZ(1.0f);
features[0].setCharge(1);
features[0].setOverallQuality(1);
features[0].setIntensity(4.0f);
//best
features[1].setRT(1.5f);
features[1].setMZ(5.0f);
features[1].setCharge(1);
features[1].setOverallQuality(1);
features[1].setIntensity(2.0f);
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

START_SECTION((virtual void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map)))
	LabeledPairFinder pm;
	Param p;
	p.setValue("rt_estimate","false");
	p.setValue("rt_pair_dist",0.4);
	p.setValue("rt_dev_low",1.0);
	p.setValue("rt_dev_high",2.0);
	p.setValue("mz_pair_dists",ListUtils::create<double>(4.0));
	p.setValue("mz_dev",0.6);
	pm.setParameters(p);

	ConsensusMap output;
	TEST_EXCEPTION(Exception::IllegalArgument,pm.run(vector<ConsensusMap>(),output));
	vector<ConsensusMap> input(1);
	MapConversion::convert(5,features,input[0]);
	output.getColumnHeaders()[5].label = "light";
	output.getColumnHeaders()[5].filename = "filename";
	output.getColumnHeaders()[8] = output.getColumnHeaders()[5];
	output.getColumnHeaders()[8].label = "heavy";

	pm.run(input,output);

	TEST_EQUAL(output.size(),1);
	ABORT_IF(output.size()!=1)
	TEST_REAL_SIMILAR(output[0].begin()->getMZ(),1.0f);
	TEST_REAL_SIMILAR(output[0].begin()->getRT(),1.0f);
	TEST_REAL_SIMILAR(output[0].rbegin()->getMZ(),5.0f);
	TEST_REAL_SIMILAR(output[0].rbegin()->getRT(),1.5f);
	TEST_REAL_SIMILAR(output[0].getQuality(),0.959346);
	TEST_EQUAL(output[0].getCharge(),1);

	//test automated RT parameter estimation
	LabeledPairFinder pm2;
	Param p2;
	p2.setValue("rt_estimate","true");
	p2.setValue("mz_pair_dists", ListUtils::create<double>(4.0));
	p2.setValue("mz_dev",0.2);
	pm2.setParameters(p2);

	FeatureMap features2;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("LabeledPairFinder.featureXML"),features2);

	ConsensusMap output2;
	vector<ConsensusMap> input2(1);
	MapConversion::convert(5,features2,input2[0]);
	output2.getColumnHeaders()[5].label = "light";
	output2.getColumnHeaders()[5].filename = "filename";
	output2.getColumnHeaders()[8] = output.getColumnHeaders()[5];
	output2.getColumnHeaders()[8].label = "heavy";
	pm2.run(input2,output2);
	TEST_EQUAL(output2.size(),250);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



