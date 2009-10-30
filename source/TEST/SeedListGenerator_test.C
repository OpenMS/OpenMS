// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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


START_SECTION((void getSeedList(const MSExperiment<>& experiment, SeedList& seeds)))
{
}
END_SECTION


START_SECTION((void getSeedList(const vector<PeptideIdentification>& peptides, SeedList& seeds)))
{
}
END_SECTION


START_SECTION((void getSeedLists(const ConsensusMap& consensus, Map<UInt64, SeedList>& seed_lists)))
{
}
END_SECTION


START_SECTION((void convert(const SeedList& seeds, FeatureMap<>& features)))
{
}
END_SECTION


START_SECTION((void convert(const FeatureMap<>& features, SeedList& seeds)))
{
}
END_SECTION


END_TEST
