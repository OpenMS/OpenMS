// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <map>

///////////////////////////

#include <OpenMS/FEATUREFINDER/SeedListGenerator.h>
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
	std::map<UInt64, SeedListGenerator::SeedList> seed_lists;
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
