// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmLabeled, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmLabeled* ptr = nullptr;
FeatureGroupingAlgorithmLabeled* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmLabeled()))
	ptr = new FeatureGroupingAlgorithmLabeled();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmLabeled()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = nullptr;
  FeatureGroupingAlgorithm* base_NullPointer = nullptr;
	ptr2 = FeatureGroupingAlgorithmLabeled::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmLabeled::getProductName(),"labeled")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap > &maps, ConsensusMap &out)))
	TOLERANCE_ABSOLUTE(0.001)

	FeatureGroupingAlgorithmLabeled fga;
	std::vector< FeatureMap > in;
	ConsensusMap out;

	//test exception (no input)
	TEST_EXCEPTION(Exception::IllegalArgument, fga.group(in,out));

	//real test
	in.resize(1);
	in[0].resize(10);
	//start
	in[0][0].setRT(1.0f);
	in[0][0].setMZ(1.0f);
	in[0][0].setCharge(1);
	in[0][0].setOverallQuality(1);
	in[0][0].setIntensity(4.0f);
	//best
	in[0][1].setRT(1.5f);
	in[0][1].setMZ(5.0f);
	in[0][1].setCharge(1);
	in[0][1].setOverallQuality(1);
	in[0][1].setIntensity(2.0f);
	//inside (down, up, left, right)
	in[0][2].setRT(1.0f);
	in[0][2].setMZ(5.0f);
	in[0][2].setCharge(1);
	in[0][2].setOverallQuality(1);

	in[0][3].setRT(3.0f);
	in[0][3].setMZ(5.0f);
	in[0][3].setCharge(1);
	in[0][3].setOverallQuality(1);

	in[0][4].setRT(1.5f);
	in[0][4].setMZ(4.8f);
	in[0][4].setCharge(1);
	in[0][4].setOverallQuality(1);

	in[0][5].setRT(1.5f);
	in[0][5].setMZ(5.2f);
	in[0][5].setCharge(1);
	in[0][5].setOverallQuality(1);

	//outside (down, up, left, right)
	in[0][6].setRT(0.0f);
	in[0][6].setMZ(5.0f);
	in[0][6].setCharge(1);
	in[0][6].setOverallQuality(1);

	in[0][7].setRT(4.0f);
	in[0][7].setMZ(5.0f);
	in[0][7].setCharge(1);
	in[0][7].setOverallQuality(1);

	in[0][8].setRT(1.5f);
	in[0][8].setMZ(4.0f);
	in[0][8].setCharge(1);
	in[0][8].setOverallQuality(1);

	in[0][9].setRT(1.5f);
	in[0][9].setMZ(6.0f);
	in[0][9].setCharge(1);
	in[0][9].setOverallQuality(1);

	Param p;
	p.setValue("rt_estimate","false");
	p.setValue("rt_pair_dist",0.4);
	p.setValue("rt_dev_low",1.0);
	p.setValue("rt_dev_high",2.0);
	p.setValue("mz_pair_dists",ListUtils::create<double>(4.0));
	p.setValue("mz_dev",0.6);
	fga.setParameters(p);


	//test exception (no file name set in out)
	TEST_EXCEPTION(Exception::IllegalArgument, fga.group(in,out));

	out.getColumnHeaders()[5].label = "light";
	out.getColumnHeaders()[5].filename = "filename";
	out.getColumnHeaders()[8] = out.getColumnHeaders()[5];
	out.getColumnHeaders()[8].label = "heavy";
	fga.group(in,out);

	TEST_EQUAL(out.size(),1)
	TEST_REAL_SIMILAR(out[0].getQuality(),0.959346);
	TEST_EQUAL(out[0].size(),2)
	ConsensusFeature::HandleSetType::const_iterator it = out[0].begin();
	TEST_REAL_SIMILAR(it->getMZ(),1.0f);
	TEST_REAL_SIMILAR(it->getRT(),1.0f);
	++it;
	TEST_REAL_SIMILAR(it->getMZ(),5.0f);
	TEST_REAL_SIMILAR(it->getRT(),1.5f);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



