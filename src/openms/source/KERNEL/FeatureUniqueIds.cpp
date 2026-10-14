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

#include <string>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    struct IdTally
    {
      std::unordered_set<UInt64> seen;
      Size invalid = 0;
      bool has_duplicate = false;
      UInt64 duplicate = 0;

      void add(UInt64 id)
      {
        if (!UniqueIdInterface::isValid(id))
        {
          ++invalid;
        }
        else if (!seen.insert(id).second && !has_duplicate)
        {
          has_duplicate = true;
          duplicate = id;
        }
      }

      void addFeature(const Feature& feature)
      {
        add(feature.getUniqueId());
        for (const auto& subordinate : feature.getSubordinates())
        {
          addFeature(subordinate);
        }
      }

      void throwIfViolated(const char* what) const
      {
        if (invalid > 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::to_string(invalid) + " " + what + " without a valid unique id.",
                                        "0");
        }
        if (has_duplicate)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Unique id shared by several ") + what + ".",
                                        std::to_string(duplicate));
        }
      }
    };
  } // namespace

  void FeatureUniqueIds::check(const FeatureMap& map)
  {
    IdTally tally;
    for (const auto& feature : map)
    {
      tally.addFeature(feature);
    }
    tally.throwIfViolated("features (including subordinates)");
  }

  void FeatureUniqueIds::check(const ConsensusMap& map)
  {
    IdTally tally;
    for (const auto& feature : map)
    {
      tally.add(feature.getUniqueId());
    }
    tally.throwIfViolated("consensus features");
  }
} // namespace OpenMS
