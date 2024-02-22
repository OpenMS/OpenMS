// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/LogStream.h>

// NNLS isotope correction
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <memory>

// #define ISOBARIC_QUANT_DEBUG

namespace OpenMS
{

  void
  IsobaricIsotopeCorrector::correctIsotopicImpurities(std::vector<double> & intensities,
                                                     const IsobaricQuantitationMethod* quant_method)
  {
    std::vector<double> res(quant_method->getNumberOfChannels());
    // we need to copy anyway because NNLS will modify the input Matrix
    MutableEigenMatrixXdPtr m(convertOpenMSMatrix2MutableEigenMatrixXd(quant_method->getIsotopeCorrectionMatrix()));
    solveNNLS_(m, intensities, res);
    intensities = res;
  }

  /* That is the safe version with dozen of copies
  std::vector<double>
  IsobaricIsotopeCorrector::correctIsotopicImpurities(std::vector<double> & intensities,
                                                     const IsobaricQuantitationMethod* quant_method)
  {
    Size c = quant_method->getNumberOfChannels();
    Matrix<double> m(quant_method->getIsotopeCorrectionMatrix());
    std::vector<double> res(c, 0.);

    // data structures for NNLS
    Matrix<double> m_b(c, 1.);
    for (Size i = 0; i < c; ++i)
    {
      m_b(i,0) = intensities[i];
    }

    Matrix<double> m_x(c, 1.);
    
    solveNNLS_(m, m_b, m_x);
    for (Size i = 0; i < c; ++i)
    {
      res[i] = m_x(i,0);
    }
    return res;
  }*/

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

    if (correction_matrix.getEigenMatrix().isIdentity(0.0))
    {
      OPENMS_LOG_DEBUG << "Correction matrix is the identity matrix." << std::endl;
      OPENMS_LOG_DEBUG << correction_matrix << std::endl;

      // workaround: TMT11plex has a special case where the correction matrix is the identity matrix
      if (quant_method->getMethodName() != "tmt11plex")
      {        
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "IsobaricIsotopeCorrector: The given isotope correction matrix is an identity matrix leading to no correction. "
                                          "Please provide a valid isotope_correction matrix as it was provided with the sample kit!");
      }
    }
    
    Eigen::FullPivLU<Eigen::MatrixXd> ludecomp(correction_matrix.getEigenMatrix());
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
      if (!(correction_matrix.getEigenMatrix() * e_mx).isApprox(b))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Cannot multiply!");
      }
      solveNNLS_(correction_matrix, m_b, m_x);

      // update the output consensus map with the corrected intensities
      float cf_intensity = updateOutputMap_(consensus_map_in, consensus_map_out, i, m_x);

      // check consistency
      computeStats_(m_x, e_mx, cf_intensity, quant_method, stats);

      /* Try this when time permits
      
      std::vector<double> b = getIntensities_(quant_method, consensus_map_in[i], consensus_map_in);

      //solve
      const auto& b_eigen = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
      Eigen::MatrixXd e_x = ludecomp.solve(b_eigen);
      if (!((*m) * e_x).isApprox(b_eigen))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Cannot multiply!");
      }
      solveNNLS_(m, b, x);

      // update the output consensus map with the corrected intensities
      float cf_intensity = updateOutputMap_(consensus_map_in, consensus_map_out, i, x);

      // check consistency
      computeStats_(x, e_x, cf_intensity, quant_method, stats);
      */
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

  std::vector<double>
  IsobaricIsotopeCorrector::getIntensities_(const IsobaricQuantitationMethod* quant_method, const ConsensusFeature& cf, const ConsensusMap& cm)
  {
    int first_map_index = cf.getFeatures().begin()->getMapIndex();
    Int map_index_offset = Int(cm.getColumnHeaders().find(first_map_index)->second.getMetaValue("channel_id")) - first_map_index;
    std::vector<double> res;
    res.resize(quant_method->getNumberOfChannels());
    Int index = 0;
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = cf.getFeatures().begin();
         it_elements != cf.getFeatures().end();
         ++it_elements)
    {
      //find channel_id of current element
      index = Int(it_elements->getMapIndex() - map_index_offset);
      res[index] = it_elements->getIntensity();
    }
    return res;
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
  IsobaricIsotopeCorrector::solveNNLS_(std::shared_ptr<Eigen::MatrixXd> & correction_matrix, std::vector<double> & b, std::vector<double> & x)
  {
    Int status = NonNegativeLeastSquaresSolver::solve(correction_matrix, b, x);
    if (status != NonNegativeLeastSquaresSolver::SOLVED)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Failed to find least-squares fit!");
    }
  }

  void
  IsobaricIsotopeCorrector::solveNNLS_(std::shared_ptr<const Eigen::MatrixXd> & correction_matrix, std::vector<double> & b, std::vector<double> & x)
  {
    Eigen::MatrixXd copy = *correction_matrix;
    auto copy_ptr = std::make_shared<Eigen::MatrixXd>(copy);
    Int status = NonNegativeLeastSquaresSolver::solve(copy_ptr, b, x);
    if (status != NonNegativeLeastSquaresSolver::SOLVED)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsobaricIsotopeCorrector: Failed to find least-squares fit!");
    }
  }

  void
  IsobaricIsotopeCorrector::computeStats_(const std::vector<double>& m_x,
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
      else if ((((std::fabs(m_x[index] - x(index)))/m_x[index])*100) > 1)
      {
        ++s_different_count;
        s_different_intensity += std::fabs(m_x[index] - x(index));
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
      else if ((((std::fabs(m_x(index,0) - x(index)))/m_x(index,0))*100) > 1)
      {
        ++s_different_count;
        s_different_intensity += std::fabs(m_x(index,0) - x(index));
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
  IsobaricIsotopeCorrector::updateOutputMap_(
    const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out,
    ConsensusMap::size_type current_cf, const std::vector<double>& m_x)
  {
    float cf_intensity(0);
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[current_cf].begin();
         it_elements != consensus_map_in[current_cf].end();
         ++it_elements)
    {
      FeatureHandle handle = *it_elements;
      //find channel_id of current element
      Int index = Int(consensus_map_out.getColumnHeaders()[it_elements->getMapIndex()].getMetaValue("channel_id"));
      handle.setIntensity(float(m_x[index]));

      consensus_map_out[current_cf].insert(handle);
      cf_intensity += handle.getIntensity(); // sum up all channels for CF

#ifdef ISOBARIC_QUANT_DEBUG
      std::cout <<  it_elements->getIntensity() << " -> " << handle.getIntensity() << std::endl;
#endif
    }
    consensus_map_out[current_cf].setIntensity(cf_intensity); // set overall intensity of CF (sum of all channels)

    return cf_intensity;
  }

    float
  IsobaricIsotopeCorrector::updateOutputMap_(
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
      handle.setIntensity(float(m_x(index,0)));

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
