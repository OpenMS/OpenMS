// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ISOBARICISOTOPECORRECTOR_H
#define OPENMS_ANALYSIS_QUANTITATION_ISOBARICISOTOPECORRECTOR_H

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

namespace OpenMS
{
  /**
    @brief Performs isotope impurity correction on the intensities extracted from an isobaric labeling experiment.
  */
  class OPENMS_DLLAPI IsobaricIsotopeCorrector
  {
public:
    /**
      @brief Constructor given an IsobaricQuantitationMethod (e.g., iTRAQ 4 plex).

      @param quant_method The quantification method used for the data set to analyze.
     */
    IsobaricIsotopeCorrector(const IsobaricQuantitationMethod* const quant_method);

    /// Copy c'tor
    IsobaricIsotopeCorrector(const IsobaricIsotopeCorrector& other);

    /// Assignment operator
    IsobaricIsotopeCorrector& operator=(const IsobaricIsotopeCorrector& rhs);

    virtual ~IsobaricIsotopeCorrector();
    
    /**
     @brief Apply isotope correction to the given input map and store the corrected values in the output map.

     @param consensus_map_in The map containing the values that should be corrected.
     @param consensus_map_out The map where the corrected values should be stored.

     @throws Exception::FailedAPICall If the least-squares fit fails.
     @throws Exception::InvalidParameter If the given correction matrix is invalid.
     */
    IsobaricQuantifierStatistics correctIsotopicImpurities(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out);

private:
    /// The quantification method used for the dataset to be analyzed.
    const IsobaricQuantitationMethod* quant_method_;
    
    /// @brief GSL objects used for the isotope correction.
    /// @{
    gsl_matrix* gsl_m_;
    gsl_permutation* gsl_p_;
    gsl_vector* gsl_b_;
    gsl_vector* gsl_x_;
    
    /// Indicates wether memory was allocated for the gsl vector/matrix pointers.
    bool gsl_allocated_;
    
    /// Free all memory allocated by GSL objects.
    void freeGSLMemory_();
    /// @}
    
    /**
     @brief Checks if the given matrix is an identity matrix.
     
     @param channel_frequency The matrix to check.
     @return True if the given matrix is an identity matrix, false otherwise.
     */
    bool isIdentityMatrix_(const Matrix<double>& channel_frequency) const;
    
    /**
     @brief Checks if the given matrix is invertible.
     
     @param channel_frequency The matrix to test.
     @return True if the matrix is invertible, false otherwise.
     */
    bool isInvertible_(Matrix<double>& channel_frequency) const;
    
    /**
     @brief Fills the input vector for the gsl/NNLS step given the ConsensusFeature.
     */
    void fillInputVector_(gsl_vector* gsl_b, Matrix<double>& m_b, const ConsensusFeature& cf, const ConsensusMap& cm) const;
    
    /**
     @brief Solves the
     */
    void solveGSL_(const gsl_matrix* gsl_m, const gsl_permutation* gsl_p, const gsl_vector* gsl_b, gsl_vector* gsl_x) const;
    
    /**
     @brief
     */
    void solveNNLS_(const Matrix<double>& correction_matrix, const Matrix<double>& m_b, Matrix<double>& m_x) const;
    
    /**
     @brief
     */
    void computeStats_(const Matrix<double>& m_x, gsl_vector* gsl_x, const ConsensusFeature::IntensityType cf_intensity, IsobaricQuantifierStatistics& stats);
    
    /**
     @brief
     */
    ConsensusFeature::IntensityType updateOutpuMap_(const ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, ConsensusMap::size_type current_cf, const Matrix<double> & m_x) const;
  };
} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ISOBARICISOTOPECORRECTOR_H
