// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

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
    const vector<PeptideIdentification>& peptides =
      feature.getPeptideIdentifications();
    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty())
      {
        continue; // shouldn't be the case
      }
      annotations_.insert(pep_it->getHits()[0].getSequence());
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
