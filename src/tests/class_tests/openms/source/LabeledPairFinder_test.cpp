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

START_SECTION((static BaseGroupFinder* create()))
	BaseGroupFinder* base_ptr = nullptr;
  BaseGroupFinder* base_nullPointer = nullptr;
	base_ptr = LabeledPairFinder::create();
  TEST_NOT_EQUAL(base_ptr, base_nullPointer)
END_SECTION

START_SECTION((static const String getProductName()))
  LabeledPairFinder spf;

  TEST_STRING_EQUAL(spf.getProductName(),"labeled_pair_finder")
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
	output.getFileDescriptions()[5].label = "light";
	output.getFileDescriptions()[5].filename = "filename";
	output.getFileDescriptions()[8] = output.getFileDescriptions()[5];
	output.getFileDescriptions()[8].label = "heavy";

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
	output2.getFileDescriptions()[5].label = "light";
	output2.getFileDescriptions()[5].filename = "filename";
	output2.getFileDescriptions()[8] = output.getFileDescriptions()[5];
	output2.getFileDescriptions()[8].label = "heavy";
	pm2.run(input2,output2);
	TEST_EQUAL(output2.size(),250);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



