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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include <OpenMS/CONCEPT/Factory.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class FGA
	 : public FeatureGroupingAlgorithm
	{
		public:
			void group(const vector< FeatureMap >&, ConsensusMap& map) override
			{
			  map.getFileDescriptions()[0].filename = "bla";
				map.getFileDescriptions()[0].size = 5;
			}
	};
}

START_TEST(FeatureGroupingAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FGA* ptr = nullptr;
FGA* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithm()))
	ptr = new FGA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const vector< FeatureMap > &maps, ConsensusMap &out)=0))
	FGA fga;
	vector< FeatureMap > in;
	ConsensusMap map;
	fga.group(in,map);
	TEST_EQUAL(map.getFileDescriptions()[0].filename, "bla")
END_SECTION

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[0],FeatureGroupingAlgorithmLabeled::getProductName());
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[1],FeatureGroupingAlgorithmUnlabeled::getProductName());
  TEST_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts().size(), 4)
}
END_SECTION

START_SECTION((void transferSubelements(const vector<ConsensusMap>& maps, ConsensusMap& out) const))
{
	vector<ConsensusMap> maps(2);
	maps[0].getFileDescriptions()[0].filename = "file1";
	maps[0].getFileDescriptions()[0].size = 1;
	maps[0].getFileDescriptions()[1].filename = "file2";
	maps[0].getFileDescriptions()[1].size = 1;
	maps[1].getFileDescriptions()[0].filename = "file3";
	maps[1].getFileDescriptions()[0].size = 1;
	maps[1].getFileDescriptions()[1].filename = "file4";
	maps[1].getFileDescriptions()[1].size = 1;

	Feature feat1, feat2, feat3, feat4;

  FeatureHandle handle1(0, feat1), handle2(1, feat2), handle3(0, feat3),
		handle4(1, feat4);

	maps[0].resize(1);
	maps[0][0].insert(handle1);
	maps[0][0].insert(handle2);
	maps[0][0].setUniqueId(1);
	maps[1].resize(1);
	maps[1][0].insert(handle3);
	maps[1][0].insert(handle4);
	maps[1][0].setUniqueId(2);

	ConsensusMap out;
	FeatureHandle handle5(0, static_cast<BaseFeature>(maps[0][0]));
	FeatureHandle handle6(1, static_cast<BaseFeature>(maps[1][0]));
	out.resize(1);
	out[0].insert(handle5);
	out[0].insert(handle6);

	// need an instance of FeatureGroupingAlgorithm:
	String algo_name = Factory<FeatureGroupingAlgorithm>::registeredProducts()[0];
	FeatureGroupingAlgorithm* algo = Factory<FeatureGroupingAlgorithm>::create(
		algo_name);

	algo->transferSubelements(maps, out);

	TEST_EQUAL(out.getFileDescriptions().size(), 4);
	TEST_EQUAL(out.getFileDescriptions()[0].filename, "file1");
	TEST_EQUAL(out.getFileDescriptions()[3].filename, "file4");
	TEST_EQUAL(out.size(), 1);
	TEST_EQUAL(out[0].size(), 4);

	ConsensusFeature::HandleSetType group = out[0].getFeatures();
	ConsensusFeature::HandleSetType::const_iterator it = group.begin();
	handle3.setMapIndex(2);
	handle4.setMapIndex(3);
	TEST_EQUAL(*it++ == handle1, true);
	TEST_EQUAL(*it++ == handle2, true);
	TEST_EQUAL(*it++ == handle3, true);
	TEST_EQUAL(*it++ == handle4, true);
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
