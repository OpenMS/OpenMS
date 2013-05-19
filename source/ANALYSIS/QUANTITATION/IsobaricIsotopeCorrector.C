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

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>

// NNLS isotope correction
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>

// LINALG methods for isotope correction
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// #define ISOBARIC_QUANT_DEBUG

namespace OpenMS
{

  IsobaricIsotopeCorrector::IsobaricIsotopeCorrector(const IsobaricQuantitationMethod* const quant_method) :
    quant_method_(quant_method),
    gsl_m_(0),
    gsl_p_(0),
    gsl_b_(0),
    gsl_x_(0),
    gsl_allocated_(false)
  {
  }

  IsobaricIsotopeCorrector::IsobaricIsotopeCorrector(const IsobaricIsotopeCorrector& other)
  {
    quant_method_ = other.quant_method_;
  }

  IsobaricIsotopeCorrector& IsobaricIsotopeCorrector::operator=(const IsobaricIsotopeCorrector& rhs)
  {
    if (this == &rhs)
      return *this;

    quant_method_ = rhs.quant_method_;

    return *this;
  }

  IsobaricIsotopeCorrector::~IsobaricIsotopeCorrector()
  {
    freeGSLMemory_();
  }

  void IsobaricIsotopeCorrector::freeGSLMemory_()
  {
    // ensure the memory is cleared
    if (gsl_allocated_)
    {
      gsl_matrix_free(gsl_m_);
      gsl_permutation_free(gsl_p_);
      gsl_vector_free(gsl_b_);
      gsl_vector_free(gsl_x_);
      gsl_allocated_ = false;
    }
  }

  IsobaricQuantifierStatistics IsobaricIsotopeCorrector::correctIsotopicImpurities(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out)
  {
    OPENMS_PRECONDITION(consensus_map_in.size() == consensus_map_out.size(), "The in- and output map need to have the same size.")

    // the stats object to fill while correcting
    IsobaricQuantifierStatistics stats;
    stats.number_ms2_total = consensus_map_out.size();
    stats.channel_count = quant_method_->getNumberOfChannels();

    Matrix<double> correction_matrix = quant_method_->getIsotopeCorrectionMatrix();

    if (isIdentityMatrix_(correction_matrix))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricIsotopeCorrector: The given isotope correction matrix is an identity matrix leading to no correction. "
                                                                                 "Please provide a valid isotope_correction matrix as it was provided with the sample kit!");
    }

    // convert to GSL matrix and setup required gsl datastructures
    gsl_m_ = correction_matrix.toGslMatrix();
    gsl_p_ = gsl_permutation_alloc(quant_method_->getNumberOfChannels());
    gsl_b_ = gsl_vector_alloc(quant_method_->getNumberOfChannels());
    gsl_x_ = gsl_vector_alloc(quant_method_->getNumberOfChannels());
    gsl_allocated_ = true;

    if (!isInvertible_(correction_matrix))
    {
      // clean up before we leave
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricIsotopeCorrector: The given isotope correction matrix is not invertible!");
    }

    // data structures for NNLS
    Matrix<double> m_b(quant_method_->getNumberOfChannels(), 1);
    Matrix<double> m_x(quant_method_->getNumberOfChannels(), 1);

    // correct all consensus elements
    for (ConsensusMap::size_type i = 0; i < consensus_map_out.size(); ++i)
    {
#ifdef ISOBARIC_QUANT_DEBUG
      std::cout << "\nMAP element  #### " << i << " #### \n" << std::endl;
#endif
      // delete only the consensus handles from the output map
      consensus_map_out[i].clear();

      // fill b vector
      fillInputVector_(gsl_b_, m_b, consensus_map_in[i], consensus_map_in);

      // solve using gsl and NNLS for stability and QC reasons
      solveGSL_(gsl_m_, gsl_p_, gsl_b_, gsl_x_);
      solveNNLS_(correction_matrix, m_b, m_x);

      // update the ouput consensus map with the corrected intensities
      ConsensusFeature::IntensityType cf_intensity = updateOutpuMap_(consensus_map_in, consensus_map_out, i, m_x);

      // check consistency between GSL and NNLS results
      computeStats_(m_x, gsl_x_, cf_intensity, stats);
    }

    // free all memory allocated by GSL objects
    freeGSLMemory_();

    return stats;
  }

  bool IsobaricIsotopeCorrector::isIdentityMatrix_(const Matrix<double>& channel_frequency) const
  {
    bool is_identity = true;

    for (Matrix<double>::SizeType i = 0; i < channel_frequency.rows(); ++i)
    {
      for (Matrix<double>::SizeType j = 0; j < channel_frequency.rows(); ++j)
      {
        // check if the entries are those of a identity matrix;
        // i==j -> m(i,j) == 1.0 && i!=j -> m(i,j) == 0.0
        if ((i == j && channel_frequency(i, j) != 1.0) || channel_frequency(i, j) != 0.0)
        {
          is_identity = false;
          break;
        }
      }
      // leave outer loop if we have reached the abortion cirteria
      if (!is_identity) break;
    }
    return is_identity;
  }

  bool IsobaricIsotopeCorrector::isInvertible_(Matrix<double>& channel_frequency) const
  {
    // lets see if the matrix is invertible
    int* gsl_sign = new int(0);
    int  gsl_status = gsl_linalg_LU_decomp(gsl_m_, gsl_p_, gsl_sign);
    return gsl_status == 0;
  }

  void IsobaricIsotopeCorrector::fillInputVector_(gsl_vector* gsl_b, Matrix<double>& m_b, const ConsensusFeature& cf, const ConsensusMap& cm) const
  {
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = cf.getFeatures().begin();
         it_elements != cf.getFeatures().end();
         ++it_elements)
    {
      //find channel_id of current element
      Int index = Int(cm.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_id"));
#ifdef ISOBARIC_QUANT_DEBUG
      std::cout << "  map_index " << it_elements->getMapIndex() << "-> id " << index << " with intensity " << it_elements->getIntensity() << "\n" << std::endl;
#endif
      // this is deprecated, but serves as quality measurement
      gsl_vector_set(gsl_b, index, it_elements->getIntensity());
      m_b(index, 0) = it_elements->getIntensity();
    }
  }

  void IsobaricIsotopeCorrector::solveGSL_(const gsl_matrix* gsl_m, const gsl_permutation* gsl_p, const gsl_vector* gsl_b, gsl_vector* gsl_x) const
  {
    int gsl_status = gsl_linalg_LU_solve(gsl_m, gsl_p, gsl_b, gsl_x);
    if (gsl_status != 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricIsotopeCorrector: Invalid entry in Param 'isotope_correction_values'; Cannot multiply!");
    }
  }

  void IsobaricIsotopeCorrector::solveNNLS_(const Matrix<double>& correction_matrix, const Matrix<double>& m_b, Matrix<double>& m_x) const
  {
    Int status = NonNegativeLeastSquaresSolver::solve(correction_matrix, m_b, m_x);
    if (status != NonNegativeLeastSquaresSolver::SOLVED)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricIsotopeCorrector: Failed to find least-squares fit!");
    }
  }

  void IsobaricIsotopeCorrector::computeStats_(const Matrix<double>& m_x, gsl_vector* gsl_x, const ConsensusFeature::IntensityType cf_intensity, IsobaricQuantifierStatistics& stats)
  {
    Size s_negative(0);
    Size s_different_count(0); // happens when naive solution is negative in other channels
    DoubleReal s_different_intensity(0);

    // ISOTOPE CORRECTION: compare solutions of Matrix inversion vs. NNLS
    for (Size index = 0; index < quant_method_->getNumberOfChannels(); ++index)
    {
      if (gsl_vector_get(gsl_x, index) < 0.0)
      {
        ++s_negative;
      }
      else if (std::fabs(m_x(index, 0) - gsl_vector_get(gsl_x, index)) > 0.000001)
      {
        ++s_different_count;
        s_different_intensity += std::fabs(m_x(index, 0) - gsl_vector_get(gsl_x, index));
      }
    }

    if (s_negative == 0 && s_different_count > 0) // solutions are inconsistent, despite being positive! This should not happen!
    {
      throw Exception::Postcondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricIsotopeCorrector: Isotope correction values of alternative method differ!");
    }

    // update global stats
    stats.iso_number_reporter_negative += s_negative;
    stats.iso_number_reporter_different += s_different_count;
    stats.iso_solution_different_intensity += s_different_intensity;

    if (s_negative > 0)
    {
      ++stats.iso_number_ms2_negative;
      stats.iso_total_intensity_negative += cf_intensity;
    }
  }

  ConsensusFeature::IntensityType IsobaricIsotopeCorrector::updateOutpuMap_(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out, ConsensusMap::size_type current_cf, const Matrix<double>& m_x) const
  {
    ConsensusFeature::IntensityType cf_intensity(0);
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[current_cf].begin();
         it_elements != consensus_map_in[current_cf].end();
         ++it_elements)
    {
      FeatureHandle handle = *it_elements;
      //find channel_id of current element
      Int index = Int(consensus_map_out.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_id"));
      handle.setIntensity(ConsensusFeature::IntensityType(m_x(index, 0)));

      consensus_map_out[current_cf].insert(handle);
      cf_intensity += handle.getIntensity(); // sum up all channels for CF

#ifdef ISOBARIC_QUANT_DEBUG
      std::cout <<  it_elements->getIntensity() << " -> " << handle.getIntensity() << std::endl;
#endif
    }
    consensus_map_out[current_cf].setIntensity(cf_intensity); // set overall intensity of CF (sum of all channels)

    return cf_intensity;
  }

} // namespace
