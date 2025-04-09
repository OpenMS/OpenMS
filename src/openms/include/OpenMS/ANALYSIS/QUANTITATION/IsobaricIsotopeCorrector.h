// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>

// forward decl
namespace Eigen
{
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;
}

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  class IsobaricQuantifierStatistics;
  class ConsensusMap;
  class ConsensusFeature;

  /**
    @brief Performs isotope impurity correction on the intensities extracted from an isobaric labeling experiment.
  */
  class OPENMS_DLLAPI IsobaricIsotopeCorrector
  {
public:
    /**
     @brief Apply isotope correction to the given input map and store the corrected values in the output map.

     @param consensus_map_in The map containing the values that should be corrected.
     @param consensus_map_out The map where the corrected values should be stored.
     @param quant_method IsobaricQuantitationMethod (e.g., iTRAQ 4 plex)

     @throws Exception::FailedAPICall If the least-squares fit fails.
     @throws Exception::InvalidParameter If the given correction matrix is invalid.
     */
    static IsobaricQuantifierStatistics correctIsotopicImpurities(const ConsensusMap& consensus_map_in,
                                                                  ConsensusMap& consensus_map_out,
                                                                  const IsobaricQuantitationMethod* quant_method);

private:
    /**
     @brief Fills the input vector for the Eigen/NNLS step given the ConsensusFeature.
     */
    static void fillInputVector_(Eigen::VectorXd& b,
                                 Matrix<double>& m_b,
                                 const ConsensusFeature& cf,
                                 const ConsensusMap& cm);

    /**
     @brief
     */
    static void solveNNLS_(const Matrix<double>& correction_matrix,
                           const Matrix<double>& m_b, Matrix<double>& m_x);

    /**
     @brief
     */
    static void computeStats_(const Matrix<double>& m_x,
                              const Eigen::MatrixXd& x,
                              const float cf_intensity,
                              const IsobaricQuantitationMethod* quant_method,
                              IsobaricQuantifierStatistics& stats);

    /**
     @brief
     */
    static float updateOutpuMap_(const ConsensusMap& consensus_map_in,
                                 ConsensusMap& consensus_map_out,
                                 Size current_cf,
                                 const Matrix<double>& m_x);
  };
} // namespace

