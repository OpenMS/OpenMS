// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

using namespace std;

namespace OpenMS
{

  GridFeature::GridFeature(const BaseFeature& feature, Size map_index,
                           Size feature_index) :
    feature_(feature),
    map_index_(map_index),
    feature_index_(feature_index),
    annotations_()
  {
  }

  GridFeature::GridFeature(const BaseFeature& feature, Size map_index,
                           Size feature_index, const IdentificationData& data) :
    GridFeature(feature, map_index, feature_index)
  {
    for (const auto& linked : feature.getLinkedIdentifications(data))
    {
      // the top match (the first, if none has a score), like the first hit of a sorted peptide identification:
      const auto* best = linked.getBestMatch();
      if (!best && !linked.matches.empty())
      {
        best = linked.matches.front();
      }
      if (!best || best->encoding != IdentificationData::Encoding::AA_SEQUENCE)
      {
        continue; // no peptide match
      }
      annotations_.insert(AASequence::fromString(best->representation));
    }
  }

  GridFeature::~GridFeature() = default;

  const BaseFeature& GridFeature::getFeature() const
  {
    return feature_;
  }

  Size GridFeature::getMapIndex() const
  {
    return map_index_;
  }

  Size GridFeature::getFeatureIndex() const
  {
    return feature_index_;
  }

  Int GridFeature::getID() const
  {
    return (Int)feature_index_;
  }

  const set<AASequence>& GridFeature::getAnnotations() const
  {
    return annotations_;
  }

  double GridFeature::getRT() const
  {
    return feature_.getRT();
  }

  double GridFeature::getMZ() const
  {
    return feature_.getMZ();
  }

}
