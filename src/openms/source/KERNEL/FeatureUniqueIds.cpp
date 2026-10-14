// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureUniqueIds.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <algorithm>
#include <string>
#include <vector>

namespace OpenMS
{
  namespace
  {
    void collectIds(const Feature& feature, std::vector<UInt64>& ids)
    {
      ids.push_back(feature.getUniqueId());
      for (const auto& subordinate : feature.getSubordinates())
      {
        collectIds(subordinate, ids);
      }
    }

    // sorting is about three times faster than hashing for maps with a million features
    void checkIds(std::vector<UInt64>& ids, const char* what)
    {
      std::sort(ids.begin(), ids.end());
      // the invalid id is 0, so invalid ids come first
      const auto first_valid = std::find_if(ids.begin(), ids.end(), [](UInt64 id) { return UniqueIdInterface::isValid(id); });
      const auto invalid = std::distance(ids.begin(), first_valid);
      if (invalid > 0)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::to_string(invalid) + " " + what + " without a valid unique id.",
                                      "0");
      }
      const auto duplicate = std::adjacent_find(ids.begin(), ids.end());
      if (duplicate != ids.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string("Unique id shared by several ") + what + ".",
                                      std::to_string(*duplicate));
      }
    }
  } // namespace

  void FeatureUniqueIds::check(const FeatureMap& map)
  {
    std::vector<UInt64> ids;
    ids.reserve(map.size());
    for (const auto& feature : map)
    {
      collectIds(feature, ids);
    }
    checkIds(ids, "features (including subordinates)");
  }

  void FeatureUniqueIds::check(const ConsensusMap& map)
  {
    std::vector<UInt64> ids;
    ids.reserve(map.size());
    for (const auto& feature : map)
    {
      ids.push_back(feature.getUniqueId());
    }
    checkIds(ids, "consensus features");
  }
} // namespace OpenMS
