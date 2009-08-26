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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition<2> PositionType;

START_TEST(StablePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StablePairFinder* ptr = 0;
START_SECTION((StablePairFinder()))
	ptr = new StablePairFinder();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~StablePairFinder()))
	delete ptr;
END_SECTION

START_SECTION((static BaseGroupFinder* create()))
	BaseGroupFinder* base_ptr = 0;
	base_ptr = StablePairFinder::create();
	TEST_NOT_EQUAL(base_ptr, 0)
END_SECTION

START_SECTION((static const String getProductName()))
	StablePairFinder spf;

  TEST_EQUAL(spf.getName() == "stable",true)
END_SECTION

START_SECTION((void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map)))
{

  std::vector<ConsensusMap> input(2);
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0f);
  feat2.setPosition(pos2);
  feat2.setIntensity(300.0f);
  feat3.setPosition(pos3);
  feat3.setIntensity(400.0f);
  ConsensusFeature cons1(0,0,feat1);
  ConsensusFeature cons2(0,1,feat2);
  ConsensusFeature cons3(0,2,feat3);
  input[0].push_back(cons1);
  input[0].push_back(cons2);
  input[0].push_back(cons3);

  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,0.04);
  PositionType pos5(204,300.04);
  PositionType pos6(404,500.04);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  feat5.setPosition(pos5);
  feat5.setIntensity(300.0f);
  feat6.setPosition(pos6);
  feat6.setIntensity(400.0f);
  ConsensusFeature cons4(1,0,feat4);
  ConsensusFeature cons5(1,1,feat5);
  ConsensusFeature cons6(1,2,feat6);
  input[1].push_back(cons4);
  input[1].push_back(cons5);
  input[1].push_back(cons6);

  StablePairFinder spf;
	Param param = spf.getDefaults();
//	param.setValue("similarity:max_pair_distance:RT",50.0);
//	param.setValue("similarity:max_pair_distance:MZ",5.0);
//	param.setValue("similarity:second_nearest_gap",  2.0);
	spf.setParameters(param);
	ConsensusMap result;
	spf.run(input,result);
	TEST_EQUAL(result.size(),3);
	ABORT_IF(result.size()!=3);

  ConsensusFeature::HandleSetType group1 = result[0].getFeatures();
  ConsensusFeature::HandleSetType group2 = result[1].getFeatures();
  ConsensusFeature::HandleSetType group3 = result[2].getFeatures();

  FeatureHandle ind1(0,0,feat1);
  FeatureHandle ind2(0,1,feat2);
  FeatureHandle ind3(0,2,feat3);
  FeatureHandle ind4(1,0,feat4);
  FeatureHandle ind5(1,1,feat5);
  FeatureHandle ind6(1,2,feat6);

  ConsensusFeature::HandleSetType::const_iterator it;
	it = group1.begin();
  STATUS(*it);
	STATUS(ind1);
	TEST_EQUAL(*(it) == ind1, true)
	++it;
  STATUS(*it);
	STATUS(ind4);
  TEST_EQUAL(*(it) == ind4, true)
	it = group2.begin();
  STATUS(*it);
	STATUS(ind2);
  TEST_EQUAL(*(it) == ind2, true)
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true)
  it = group3.begin();
  STATUS(*it);
	STATUS(ind3);
  TEST_EQUAL(*(it) == ind3, true)
	++it;
  STATUS(*it);
	STATUS(ind6);
  TEST_EQUAL(*(it) == ind6, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
