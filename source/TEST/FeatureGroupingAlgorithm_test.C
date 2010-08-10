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
			void group(const std::vector< FeatureMap<> >&, ConsensusMap& map)
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
START_SECTION((FeatureGroupingAlgorithm()))
	ptr = new FGA();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap<> > &maps, ConsensusMap &out)=0))
	FGA fga;
	std::vector< FeatureMap<> > in;
	ConsensusMap map; 
	fga.group(in,map);
	TEST_EQUAL(map.getFileDescriptions()[0].filename, "bla")
END_SECTION

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[0],FeatureGroupingAlgorithmLabeled::getProductName());
	TEST_STRING_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts()[1],FeatureGroupingAlgorithmUnlabeled::getProductName());
	TEST_EQUAL(Factory<FeatureGroupingAlgorithm>::registeredProducts().size(), 2)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
