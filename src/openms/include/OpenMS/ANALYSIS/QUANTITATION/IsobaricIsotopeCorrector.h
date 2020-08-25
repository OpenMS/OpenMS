// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <Eigen/Core>

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  struct IsobaricQuantifierStatistics;
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
     @param IsobaricQuantitationMethod (e.g., iTRAQ 4 plex)

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

