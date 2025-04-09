// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

namespace OpenMS
{

  // TODO consider using (RT,MZ,IT) as a unique identifier ?
  bool MapAlignmentEvaluationAlgorithm::isSameHandle(const FeatureHandle& lhs, const FeatureHandle& rhs, const double& rt_dev, const double& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge)
  {
#if 1
    // use (RT,MZ,IT) as "unique" identifier
    if (fabs(lhs.getRT() - rhs.getRT()) > rt_dev)
      return false; // TODO MAGIC_ALERT

    if (fabs(lhs.getMZ() - rhs.getMZ()) > mz_dev)
      return false; // TODO MAGIC_ALERT

    if (fabs(lhs.getIntensity() - rhs.getIntensity()) > int_dev)
      return false; // TODO MAGIC_ALERT

    if (use_charge && (lhs.getCharge() != rhs.getCharge()))
      return false;

    return true;

#else
    // use (map index, element index) as unique identifier
    return lhs.getMapIndex() == rhs.getMapIndex() && lhs.getElementIndex() == rhs.getElementIndex();

#endif
  }

  MapAlignmentEvaluationAlgorithm::MapAlignmentEvaluationAlgorithm() = default;

  MapAlignmentEvaluationAlgorithm::~MapAlignmentEvaluationAlgorithm() = default;

}
