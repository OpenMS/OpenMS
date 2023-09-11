// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

namespace OpenMS
{
  /**
      @brief Caap evaluation algorithm to obtain a precision value.

      It evaluates an input consensus map with respect to a ground truth.

      @ingroup MapAlignmentEvaluation
  */
  class OPENMS_DLLAPI MapAlignmentEvaluationAlgorithmPrecision :
    public MapAlignmentEvaluationAlgorithm
  {
public:
    /// Default constructor
    MapAlignmentEvaluationAlgorithmPrecision();

    /// Destructor
    ~MapAlignmentEvaluationAlgorithmPrecision() override;

    /**
        @brief Applies the algorithm
    */
    void evaluate(const ConsensusMap & consensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out) override;

    /// Creates a new instance of this class (for Factory)
    static MapAlignmentEvaluationAlgorithm * create()
    {
      return new MapAlignmentEvaluationAlgorithmPrecision();
    }

    /// Returns the product name (for the Factory)
    static String getProductName()
    {
      return "precision";
    }

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentEvaluationAlgorithmPrecision(const MapAlignmentEvaluationAlgorithmPrecision &);
    /// Assignment operator intentionally not implemented -> private
    MapAlignmentEvaluationAlgorithmPrecision & operator=(const MapAlignmentEvaluationAlgorithmPrecision &);

  };

} // namespace OpenMS

