// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmIdentification.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class FGA
	 : public FeatureGroupingAlgorithm
	{
		public:
			void group(const vector< FeatureMap<> >&, ConsensusMap& map)
			{
			  map.getFileDescriptions()[0].filename = "bla";
				map.getFileDescriptions()[0].size = 5;
			}
	};
}

START_TEST(FeatureGroupingAlgorithm, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FGA* ptr = 0;
FGA* nullPointer = 0;
START_SECTION((FeatureGroupingAlgorithm()))
	ptr = new FGA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const vector< FeatureMap<> > &maps, ConsensusMap &out)=0))
	FGA fga;
	vector< FeatureMap<> > in;
	ConsensusMap map;
	fga.group(in,map);
	TEST_EQUAL(map.getFileDescriptions()[0].filename, "bla")
END_SECTION

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[0],FeatureGroupingAlgorithmLabeled::getProductName());
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[1],FeatureGroupingAlgorithmUnlabeled::getProductName());
	TEST_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts().size(), 3)
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
