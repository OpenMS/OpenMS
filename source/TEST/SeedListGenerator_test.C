// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SeedListGenerator, "$Id$")

/////////////////////////////////////////////////////////////

SeedListGenerator* slg_ptr = 0;

START_SECTION((SeedListGenerator()))
{
	slg_ptr = new SeedListGenerator();
  TEST_NOT_EQUAL(slg_ptr, 0);
}
END_SECTION


START_SECTION(([EXTRA] ~SeedListGenerator()))
{
	delete slg_ptr;
}
END_SECTION


START_SECTION((void generateSeedList(const MSExperiment<>& experiment, SeedList& seeds)))
{
	MSExperiment<> experiment;
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
	peptides[0].setMetaValue("RT", 1.1);
	peptides[0].setMetaValue("MZ", 111.111);
	peptides[1].setMetaValue("RT", 2.2);
	peptides[1].setMetaValue("MZ", 222.222);
	peptides[2].setMetaValue("RT", 3.3);
	peptides[2].setMetaValue("MZ", 333.333);
	SeedListGenerator::SeedList seeds;
	SeedListGenerator().generateSeedList(peptides, seeds);
	TEST_EQUAL(seeds.size(), 3);
	TEST_EQUAL(seeds[0], DPosition<2>(1.1, 111.111));
	TEST_EQUAL(seeds[1], DPosition<2>(2.2, 222.222));
	TEST_EQUAL(seeds[2], DPosition<2>(3.3, 333.333));
	PeptideHit hit;
	hit.setSequence("TEST");
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


START_SECTION((void convertSeedList(const SeedList& seeds, FeatureMap<>& features)))
{
	SeedListGenerator::SeedList seeds(3);
	seeds[0] = DPosition<2>(1.1, 111.111);
	seeds[1] = DPosition<2>(2.2, 222.222);
	seeds[2] = DPosition<2>(3.3, 333.333);
	FeatureMap<> features;
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


START_SECTION((void convertSeedList(const FeatureMap<>& features, SeedList& seeds)))
{
	FeatureMap<> features;
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
