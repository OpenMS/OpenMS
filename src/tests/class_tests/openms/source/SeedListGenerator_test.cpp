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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SeedListGenerator, "$Id$")

/////////////////////////////////////////////////////////////

SeedListGenerator* slg_ptr = nullptr;
SeedListGenerator* slg_nullPointer = nullptr;

START_SECTION((SeedListGenerator()))
{
	slg_ptr = new SeedListGenerator();
  TEST_NOT_EQUAL(slg_ptr, slg_nullPointer);
}
END_SECTION


START_SECTION(([EXTRA] ~SeedListGenerator()))
{
	delete slg_ptr;
}
END_SECTION


START_SECTION((void generateSeedList(const PeakMap& experiment, SeedList& seeds)))
{
	PeakMap experiment;
	String path = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.mzML");
	MzMLFile().load(path, experiment);
	SeedListGenerator::SeedList seeds;
	SeedListGenerator().generateSeedList(experiment, seeds);
	TEST_EQUAL(seeds.size(), 9);
	TEST_EQUAL(seeds[0], DPosition<2>(0.5927, 538.605));
	TEST_EQUAL(seeds[1], DPosition<2>(0.5927, 637.885));
	TEST_EQUAL(seeds[2], DPosition<2>(0.5927, 678.384));
	// ...
	TEST_EQUAL(seeds[8], DPosition<2>(3.7572, 512.784));
}
END_SECTION


START_SECTION((void generateSeedList(vector<PeptideIdentification>& peptides, SeedList& seeds, bool use_peptide_mass = false)))
{
	vector<PeptideIdentification> peptides(3);
	peptides[0].setRT(1.1);
	peptides[0].setMZ(111.111);
	peptides[1].setRT(2.2);
	peptides[1].setMZ(222.222);
	peptides[2].setRT(3.3);
	peptides[2].setMZ(333.333);
	SeedListGenerator::SeedList seeds;
	SeedListGenerator().generateSeedList(peptides, seeds);
	TEST_EQUAL(seeds.size(), 3);
	TEST_EQUAL(seeds[0], DPosition<2>(1.1, 111.111));
	TEST_EQUAL(seeds[1], DPosition<2>(2.2, 222.222));
	TEST_EQUAL(seeds[2], DPosition<2>(3.3, 333.333));
	PeptideHit hit;
	hit.setSequence(AASequence::fromString("TEST"));
	hit.setCharge(2);
	peptides[0].insertHit(hit);
	peptides.resize(1);
	SeedListGenerator().generateSeedList(peptides, seeds, true);
	TEST_EQUAL(seeds.size(), 1);
	TEST_REAL_SIMILAR(seeds[0][1], 219.09755);
}
END_SECTION


START_SECTION((void generateSeedLists(const ConsensusMap& consensus, Map<UInt64, SeedList>& seed_lists)))
{
	ConsensusMap consensus;
	String path = OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML");
	ConsensusXMLFile().load(path, consensus);
	Map<UInt64, SeedListGenerator::SeedList> seed_lists;
	SeedListGenerator().generateSeedLists(consensus, seed_lists);
	TEST_EQUAL(seed_lists.size(), 2);
	TEST_EQUAL(seed_lists[0].size(), 0);
	TEST_EQUAL(seed_lists[1].size(), 2);
	TEST_EQUAL(seed_lists[1][0], DPosition<2>(1273.27, 904.47));
	TEST_EQUAL(seed_lists[1][1], DPosition<2>(1184.46, 953.368));
}
END_SECTION


START_SECTION((void convertSeedList(const SeedList& seeds, FeatureMap& features)))
{
	SeedListGenerator::SeedList seeds(3);
	seeds[0] = DPosition<2>(1.1, 111.111);
	seeds[1] = DPosition<2>(2.2, 222.222);
	seeds[2] = DPosition<2>(3.3, 333.333);
	FeatureMap features;
	SeedListGenerator().convertSeedList(seeds, features);
	TEST_EQUAL(features.size(), 3);
	TEST_EQUAL(features[0].getRT(), 1.1);
	TEST_EQUAL(features[0].getMZ(), 111.111);
	TEST_EQUAL(features[1].getRT(), 2.2);
	TEST_EQUAL(features[1].getMZ(), 222.222);
	TEST_EQUAL(features[2].getRT(), 3.3);
	TEST_EQUAL(features[2].getMZ(), 333.333);
}
END_SECTION


START_SECTION((void convertSeedList(const FeatureMap& features, SeedList& seeds)))
{
	FeatureMap features;
	features.resize(3);
	features[0].setRT(1.1);
	features[0].setMZ(111.111);
	features[1].setRT(2.2);
	features[1].setMZ(222.222);
	features[2].setRT(3.3);
	features[2].setMZ(333.333);
	SeedListGenerator::SeedList seeds;
	SeedListGenerator().convertSeedList(features, seeds);
	TEST_EQUAL(seeds.size(), 3);
	TEST_EQUAL(seeds[0], DPosition<2>(1.1, 111.111));
	TEST_EQUAL(seeds[1], DPosition<2>(2.2, 222.222));
	TEST_EQUAL(seeds[2], DPosition<2>(3.3, 333.333));
}
END_SECTION


END_TEST
