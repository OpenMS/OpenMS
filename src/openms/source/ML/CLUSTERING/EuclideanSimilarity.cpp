// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/EuclideanSimilarity.h>

namespace OpenMS
{
  EuclideanSimilarity::EuclideanSimilarity() :
    scale_(1)
  {
  }

  EuclideanSimilarity::EuclideanSimilarity(const EuclideanSimilarity & source) = default;

  EuclideanSimilarity::~EuclideanSimilarity() = default;

  EuclideanSimilarity & EuclideanSimilarity::operator=(const EuclideanSimilarity & source)
  {
    if (this != &source)
    {
      scale_ = source.scale_;
    }
    return *this;
  }

  float EuclideanSimilarity::operator()(const std::pair<float, float> & c) const
  {
    return operator()(c, c);
  }

  // calculates euclidean distance between two points
  float EuclideanSimilarity::operator()(const std::pair<float, float> & a, const std::pair<float, float> & b) const
  {
    if (scale_ == 0)
    {
      //inapplicable scaling
      throw Exception::DivisionByZero(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    return 1 - (sqrtf((a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second)) / scale_);
  }

  void EuclideanSimilarity::setScale(float x)
  {
    scale_ = x;
  }

}
