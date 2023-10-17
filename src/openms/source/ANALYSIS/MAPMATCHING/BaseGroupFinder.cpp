// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{

  BaseGroupFinder::BaseGroupFinder() :
    DefaultParamHandler("BaseGroupFinder")
  {
  }

  BaseGroupFinder::~BaseGroupFinder() = default;

  void BaseGroupFinder::registerChildren()
  {
    Factory<BaseGroupFinder>::registerProduct(
      SimplePairFinder::getProductName(), &SimplePairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      LabeledPairFinder::getProductName(), &LabeledPairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      StablePairFinder::getProductName(), &StablePairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      QTClusterFinder::getProductName(), &QTClusterFinder::create);
  }

  void BaseGroupFinder::checkIds_(const std::vector<ConsensusMap>& maps) const
  {
    std::set<Size> used_ids;
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap& map = maps[i];
      for (ConsensusMap::ColumnHeaders::const_iterator it = map.getColumnHeaders().begin(); it != map.getColumnHeaders().end(); ++it)
      {
        if (used_ids.find(it->first) != used_ids.end())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "file ids have to be unique");
        }
        else
        {
          used_ids.insert(it->first);
        }
      }
    }
  }

}
