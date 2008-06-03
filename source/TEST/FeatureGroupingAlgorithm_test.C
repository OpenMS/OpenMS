// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

///////////////////////////

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
			  map.setFileDescription(0);
			  map.getFileDescriptions()[0].filename = "bla";
				map.getFileDescriptions()[0].size = 5;
			}
	};
}

START_TEST(FeatureGroupingAlgorithm, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FGA* ptr = 0;
CHECK((FeatureFinderAlgorithm()))
	ptr = new FGA();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~FeatureFinderAlgorithm()))
	delete ptr;
RESULT

CHECK(virtual void group(const std::vector< FeatureMap<> >&, ConsensusMap&))
	FGA fga;
	std::vector< FeatureMap<> > in;
	ConsensusMap map; 
	fga.group(in,map);
	TEST_EQUAL(map.getFileDescriptions()[0].filename, "bla")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



