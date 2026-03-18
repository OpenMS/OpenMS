// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
      @brief Base class for all Caap evaluation algorithms

      These algorithms evaluates alignment results against a ground truth.
  */
  class OPENMS_DLLAPI MapAlignmentEvaluationAlgorithm
  {

protected:
    typedef ConsensusFeature::HandleSetType::const_iterator HandleIterator;

public:

    /// Default constructor
    MapAlignmentEvaluationAlgorithm();

    /// Destructor
    virtual ~MapAlignmentEvaluationAlgorithm();


    ///Applies the algorithm. The input consensus map is compared to the ground truth.
    virtual void evaluate(const ConsensusMap & conensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out) = 0;

    ///Decides if two features are the same, based on maximum allowed deviations for retention time, m/z and intensity.
    bool isSameHandle(const FeatureHandle & lhs, const FeatureHandle & rhs, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge);

private:
    ///Copy constructor is not implemented -> private
    MapAlignmentEvaluationAlgorithm(const MapAlignmentEvaluationAlgorithm &);
    ///Assignment operator is not implemented -> private
    MapAlignmentEvaluationAlgorithm & operator=(const MapAlignmentEvaluationAlgorithm &);

  };

} // namespace OpenMS

