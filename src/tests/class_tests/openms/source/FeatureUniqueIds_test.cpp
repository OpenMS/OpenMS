// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/FeatureUniqueIds.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;

namespace
{
  template <typename MapType>
  bool hasValidUniqueIds(const MapType& map)
  {
    try
    {
      FeatureUniqueIds::check(map);
      return true;
    }
    catch (const Exception::InvalidValue&)
    {
      return false;
    }
  }

  Feature featureWithId(UInt64 id)
  {
    Feature feature;
    feature.setUniqueId(id);
    return feature;
  }
}

START_TEST(FeatureUniqueIds, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void check(const FeatureMap& map)))
{
  FeatureMap map;
  TEST_TRUE(hasValidUniqueIds(map)); // empty map is valid

  map.push_back(featureWithId(1));
  map.push_back(featureWithId(2));
  TEST_TRUE(hasValidUniqueIds(map));

  // invalid id on a top-level feature
  FeatureMap invalid = map;
  invalid.push_back(Feature());
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(invalid));

  // duplicate id on top-level features
  FeatureMap duplicate = map;
  duplicate.push_back(featureWithId(2));
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(duplicate));

  // invalid id on a subordinate
  FeatureMap invalid_sub = map;
  invalid_sub[0].getSubordinates().push_back(Feature());
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(invalid_sub));

  // subordinate sharing its id with a top-level feature
  FeatureMap duplicate_sub = map;
  duplicate_sub[0].getSubordinates().push_back(featureWithId(2));
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(duplicate_sub));

  FeatureMap valid_sub = map;
  valid_sub[0].getSubordinates().push_back(featureWithId(3));
  TEST_TRUE(hasValidUniqueIds(valid_sub));
}
END_SECTION

START_SECTION((static void check(const ConsensusMap& map)))
{
  ConsensusMap map;
  TEST_TRUE(hasValidUniqueIds(map));

  ConsensusFeature cf;
  cf.setUniqueId(1);
  map.push_back(cf);
  cf.setUniqueId(2);
  map.push_back(cf);
  TEST_TRUE(hasValidUniqueIds(map));

  ConsensusMap duplicate = map;
  duplicate.push_back(cf);
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(duplicate));

  ConsensusMap invalid = map;
  invalid.push_back(ConsensusFeature());
  TEST_EXCEPTION(Exception::InvalidValue, FeatureUniqueIds::check(invalid));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
