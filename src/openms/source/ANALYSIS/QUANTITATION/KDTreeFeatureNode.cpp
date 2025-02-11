// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>

namespace OpenMS
{

KDTreeFeatureNode::KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx) :
  data_(data),
  idx_(idx)
{
}

KDTreeFeatureNode::KDTreeFeatureNode(const KDTreeFeatureNode& rhs) 
  
= default;

KDTreeFeatureNode& KDTreeFeatureNode::operator=(KDTreeFeatureNode const& rhs)
= default;

KDTreeFeatureNode::~KDTreeFeatureNode()
= default;

Size KDTreeFeatureNode::getIndex() const
{
  return idx_;
}

KDTreeFeatureNode::value_type KDTreeFeatureNode::operator[](Size i) const
{
  if (i == 0)
  {
    return data_->rt(idx_);
  }
  else if (i == 1)
  {
    return data_->mz(idx_);
  }
  else
  {
    const String& err_msg = "Indices other than 0 (RT) and 1 (m/z) are not allowed!";
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, err_msg);
  }
}

}
