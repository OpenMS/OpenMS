// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

#include <OpenMS/DATASTRUCTURES/Utils/MatrixUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

// NNLS isotope correction
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>

#include <Eigen/LU>

// #define ISOBARIC_QUANT_DEBUG

namespace OpenMS
{

  IsobaricQuantifierStatistics
  IsobaricIsotopeCorrector::correctIsotopicImpurities(
    const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out,
    const IsobaricQuantitationMethod* quant_method)
  {
    OPENMS_PRECONDITION(consensus_map_in.size() == consensus_map_out.size(),
                        "The in- and output map need to have the same size.")

    // the stats object to fill while correcting
    IsobaricQuantifierStatistics stats;
    stats.number_ms2_total = consensus_map_out.size();
    stats.channel_count = quant_method->getNumberOfChannels();

    Matrix<double> correction_matrix = quant_method->getIsotopeCorrectionMatrix();

    if (matrixIsIdentityMatrix(correction_matrix))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "IsobaricIsotopeCorrector: The given isotope correction matrix is an identity matrix leading to no correction. "
                                        "Please provide a valid isotope_correction matrix as it was provided with the sample kit!");
    }

    // convert to Eigen matrix
    EigenMatrixXdPtr m(convertOpenMSMatrix2EigenMatrixXd(correction_matrix));
    Eigen::FullPivLU<Eigen::MatrixXd> ludecomp(*m);
    Eigen::VectorXd b;
    b.resize(quant_method->getNumberOfChannels());
    b.setZero();
    std::vector<double> x(quant_method->getNumberOfChannels(), 0);

    if (!ludecomp.isInvertible())
    {
      // clean up before we leave
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: The given isotope correction matrix is not invertible!");
    }

    // data structures for NNLS
    Matrix<double> m_b(quant_method->getNumberOfChannels(), 1);
    Matrix<double> m_x(quant_method->getNumberOfChannels(), 1);

    // correct all consensus elements
    for (ConsensusMap::size_type i = 0; i < consensus_map_out.size(); ++i)
    {
#ifdef ISOBARIC_QUANT_DEBUG
      std::cout << "\nMAP element  #### " << i << " #### \n" << std::endl;
#endif
      // delete only the consensus handles from the output map
      consensus_map_out[i].clear();

      // fill b vector
      fillInputVector_(b, m_b, consensus_map_in[i], consensus_map_in);

      //solve
      Eigen::MatrixXd e_mx = ludecomp.solve(b);
      if (!((*m) * e_mx).isApprox(b))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Cannot multiply!");
      }
      solveNNLS_(correction_matrix, m_b, m_x);

      // update the output consensus map with the corrected intensities
      float cf_intensity = updateOutpuMap_(consensus_map_in, consensus_map_out, i, m_x);

      // check consistency
      computeStats_(m_x, e_mx, cf_intensity, quant_method, stats);    
    }

    return stats;
  }

  void
  IsobaricIsotopeCorrector::fillInputVector_(Eigen::VectorXd& b,
                                             Matrix<double>& m_b, const ConsensusFeature& cf, const ConsensusMap& cm)
  {
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = cf.getFeatures().begin();
         it_elements != cf.getFeatures().end();
         ++it_elements)
    {
      //find channel_id of current element
      Int index = Int(cm.getColumnHeaders().find(it_elements->getMapIndex())->second.getMetaValue("channel_id"));
#ifdef ISOBARIC_QUANT_DEBUG
      std::cout << "  map_index " << it_elements->getMapIndex() << "-> id " << index << " with intensity " << it_elements->getIntensity() << "\n" << std::endl;
#endif
      // this is deprecated, but serves as quality measurement
      b(index) = it_elements->getIntensity();
      m_b(index, 0) = it_elements->getIntensity();
    }
  }

  void
  IsobaricIsotopeCorrector::solveNNLS_(const Matrix<double>& correction_matrix,
                                       const Matrix<double>& m_b, Matrix<double>& m_x)
  {
    Int status = NonNegativeLeastSquaresSolver::solve(correction_matrix, m_b, m_x);
    if (status != NonNegativeLeastSquaresSolver::SOLVED)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Failed to find least-squares fit!");
    }
  }

  void
  IsobaricIsotopeCorrector::computeStats_(const Matrix<double>& m_x,
                                          const Eigen::MatrixXd& x, const float cf_intensity,
                                          const IsobaricQuantitationMethod* quant_method, IsobaricQuantifierStatistics& stats)
  {
    Size s_negative(0);
    Size s_different_count(0); // happens when naive solution is negative in other channels
    double s_different_intensity(0);

    // ISOTOPE CORRECTION: compare solutions of Matrix inversion vs. NNLS
    for (Size index = 0; index < quant_method->getNumberOfChannels(); ++index)
    {
      if (x(index) < 0.0)
      {
        ++s_negative;
      }
      else if ((((std::fabs(m_x(index, 0) - x(index)))/m_x(index, 0))*100) > 1)
      {
        ++s_different_count;
        s_different_intensity += std::fabs(m_x(index, 0) - x(index));
      }
    }

    if (s_negative == 0 && s_different_count > 0) //some solutions are inconsistent, despite being positive
    {
      OPENMS_LOG_WARN << "IsobaricIsotopeCorrector: Isotope correction values of alternative method differ!" << std::endl;
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

  float
  IsobaricIsotopeCorrector::updateOutpuMap_(
    const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out,
    ConsensusMap::size_type current_cf, const Matrix<double>& m_x)
  {
    float cf_intensity(0);
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[current_cf].begin();
         it_elements != consensus_map_in[current_cf].end();
         ++it_elements)
    {
      FeatureHandle handle = *it_elements;
      //find channel_id of current element
      Int index = Int(consensus_map_out.getColumnHeaders()[it_elements->getMapIndex()].getMetaValue("channel_id"));
      handle.setIntensity(float(m_x(index, 0)));

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
