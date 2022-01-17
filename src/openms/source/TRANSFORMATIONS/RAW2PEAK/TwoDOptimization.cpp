// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>

#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <boost/math/special_functions/acosh.hpp> // TODO replace with std::acosh

using Eigen::VectorXd;

namespace OpenMS
{

  TwoDOptimization::TwoDOptimization() :
    DefaultParamHandler("TwoDOptimization")
  {
    // 2D optimization parameters
    defaults_.setValue("penalties:position", 0.0, "If the position changes more than 0.2Da during the fitting it can be penalized");
    defaults_.setValue("penalties:height", 1.0, "penalty term for the fitting of the intensity:" \
                                                "If it gets negative during the fitting it can be penalized.");
    defaults_.setValue("penalties:left_width", 0.0, "penalty term for the fitting of the left width:" \
                                                    "If the left width gets too broad or negative during the fitting it can be penalized.");
    defaults_.setValue("penalties:right_width", 0.0, "penalty term for the fitting of the right width:" \
                                                     "If the right width gets too broad or negative during the fitting it can be penalized.");
    defaults_.setValue("2d:tolerance_mz", 2.2, "mz tolerance for cluster construction", {"advanced"});
    defaults_.setValue("2d:max_peak_distance", 1.2, "maximal peak distance in mz in a cluster", {"advanced"});
    defaults_.setValue("iterations", 10, "maximal number of iterations for the fitting step");


    defaultsToParam_();
    updateMembers_();
  }

  TwoDOptimization::TwoDOptimization(const TwoDOptimization & opt) :
    DefaultParamHandler(opt)
  {
    updateMembers_();
  }

  TwoDOptimization & TwoDOptimization::operator=(const TwoDOptimization & opt)
  {
    if (&opt == this)
    {
      return *this;
    }
    DefaultParamHandler::operator=(opt);
    updateMembers_();

    return *this;
  }

  void TwoDOptimization::findMatchingPeaks_(std::multimap<double, IsotopeCluster>::iterator & it, PeakMap & ms_exp)
  {
    IsotopeCluster::ChargedIndexSet::const_iterator iter = it->second.peaks.begin();
    for (; iter != it->second.peaks.end(); ++iter)
    {

      double mz = (ms_exp[iter->first][iter->second]).getMZ();
      mz *= 10;
      matching_peaks_[(Int)(mz + 0.5)].push_back(PeakIndex(iter->first, iter->second));
    }


#ifdef DEBUG_2D
    std::map<Int, PeakIndex>::iterator it2 = matching_peaks_.begin();
    for (; it2 != matching_peaks_.end(); ++it2)
    {
      std::cout << it2->first << " has " << it2->second.size() << " elements:" << std::endl;
      for (Size i = 0; i < it2->second.size(); ++i)
      {
        std::cout << it2->second[i]->getPeak(ms_exp).getMZ() << "\t";
      }
      std::cout << std::endl;
    }
#endif

  }

  // Finds the neighbour of the peak denoted by @p current_mz in the previous scan
  std::vector<double>::iterator TwoDOptimization::searchInScan_(std::vector<double>::iterator scan_begin,
                                                                    std::vector<double>::iterator scan_end,
                                                                    double current_mz)
  {

    // perform binary search to find the neighbour in rt dimension
    //  lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
    std::vector<double>::iterator insert_iter = lower_bound(scan_begin, scan_end, current_mz);

    // the peak found by lower_bound does not have to be the closest one, therefore we have
    // to check both neighbours
    if (insert_iter == scan_end)   // we are at the and have only one choice
    {
      return --insert_iter;
    }
    else
    {
      // if the found peak is at the beginning of the spectrum,
      // there is not much we can do.
      if (insert_iter == scan_begin)
      {
        return insert_iter;
      }
      else           // see if the next smaller one fits better
      {
        double delta_mz = fabs(*insert_iter - current_mz);
        --insert_iter;

        if (fabs(*insert_iter - current_mz) < delta_mz)
        {
          return insert_iter;                       // peak to the left is closer (in m/z dimension)
        }
        else
        {
          return ++insert_iter;                          // peak to the right is closer
        }
      }
    }

  } // end of searchInScan_

  void TwoDOptimization::updateMembers_()
  {
    penalties_.height = (double)param_.getValue("penalties:height");
    penalties_.pos = (double)param_.getValue("penalties:position");
    penalties_.lWidth = (double)param_.getValue("penalties:left_width");
    penalties_.rWidth = (double)param_.getValue("penalties:right_width");
    max_peak_distance_ = (double)param_.getValue("2d:max_peak_distance");
    tolerance_mz_ = (double)param_.getValue("2d:tolerance_mz");
    max_iteration_ = (UInt)param_.getValue("iterations");

  }


  int TwoDOptimization::TwoDOptFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
  {

    // x contains the parameters to be optimized.
    // In our case these are the weighted average mz-position and left and right width for all
    // matching peaks in all considered scans. Additionally the intensities of all peaks are stored.
    // Params might contain any additional parameters. We handle these using class members
    // instead.
    // The vector f is supposed to contain the result when we return from this function.
    double computed_signal, experimental_signal, step, last_position;
    double p_height, p_position, p_width;
    Int count = 0;
    Int counter_posf = 0;
    const std::vector<std::pair<SignedSize, SignedSize> > & signal2D = m_data->signal2D;
    std::multimap<double, IsotopeCluster>::iterator iso_map_iter = m_data->iso_map_iter;
    Size total_nr_peaks = m_data->total_nr_peaks;
    const std::map<Int, std::vector<PeakIndex> > & matching_peaks = m_data->matching_peaks;
    const PeakMap & picked_peaks = m_data->picked_peaks;
    PeakMap::ConstIterator raw_data_first = m_data->raw_data_first;
    const OptimizationFunctions::PenaltyFactorsIntensity & penalties = m_data->penalties;

    Size num_scans = signal2D.size() / 2;
    IsotopeCluster::ChargedIndexSet::iterator peak_iter = iso_map_iter->second.peaks.begin();
    fvec.setZero();

    //iterate over all scans
    for (Size current_scan = 0; current_scan < num_scans; ++current_scan)
    {
      Size curr_scan_idx = current_scan + iso_map_iter->second.peaks.begin()->first;
      double current_position = ((raw_data_first
                           + signal2D[2 * current_scan].first)->begin()
                          + signal2D[2 * current_scan].second)->getMZ();
      //iterate over all points of the signal
      for (Int current_point = 1;
           current_point +  signal2D[2 * current_scan].second
           <= signal2D[2 * current_scan + 1].second;
           ++current_point)
      {
        last_position = current_position;

        computed_signal   = 0.;
        current_position  = ((raw_data_first
                              + signal2D[2 * current_scan].first)->begin()
                             + signal2D[2 * current_scan].second + current_point)->getMZ();
        experimental_signal = ((raw_data_first
                                + signal2D[2 * current_scan].first)->begin()
                               + signal2D[2 * current_scan].second + current_point)->getIntensity();
        step = current_position - last_position;
#ifdef DEBUG_2D
        std::cout << "experimental signal rt " << (raw_data_first
                                                   + signal2D[2 * current_scan].first)->getRT()
                  << "\tmz " << ((raw_data_first
                        + signal2D[2 * current_scan].first)->begin()
                       + signal2D[2 * current_scan].second + current_point)->getMZ()
        << "\tint " << ((raw_data_first
                         + signal2D[2 * current_scan].first)->begin()
                        + signal2D[2 * current_scan].second + current_point)->getIntensity() << std::endl;
#endif

        // size_t current_peak = 0;
        peak_iter = iso_map_iter->second.peaks.begin();
        while (peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first != curr_scan_idx)
        {
          ++peak_iter;
        }
        //iterate over all peaks of the current scan
        while (peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first == curr_scan_idx)
        {
          Int peak_idx = distance(iso_map_iter->second.peaks.begin(), peak_iter);
          double mz_in_hash = ((picked_peaks[peak_iter->first]).begin() + peak_iter->second)->getMZ() * 10;
          std::map<Int, std::vector<PeakIndex> >::const_iterator  m_spec_iter = matching_peaks.begin();
          Int map_idx = 0;
          while (m_spec_iter->first != (Int)(mz_in_hash + 0.5))
          {
            ++map_idx;
            ++m_spec_iter;
          }
          //Store the current parameters for this peak
          p_position    = x(total_nr_peaks + 3 * map_idx);
          p_height      = x(peak_idx);
          p_width       = (current_position <= p_position) ?
                          x(total_nr_peaks + 3 * map_idx + 1) :
                          x(total_nr_peaks + 3 * map_idx + 2);
          ++count;



          //is it a Lorentz or a Sech - Peak?
          if ((PeakShape::Type)(Int) Math::round((picked_peaks[peak_iter->first]).getFloatDataArrays()[5][peak_iter->second]) == PeakShape::LORENTZ_PEAK)
          {
#ifdef DEBUG_2D
            std::cout << "p_height " << p_height << "\tp_position " << p_position << "\tcurrent_position "
                      << current_position << std::endl;
#endif
            computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
          }
          // if it's a Sech - Peak
          else
          {
#ifdef DEBUG_2D
            std::cout << "p_height " << p_height << "\tp_position " << p_position << "\tcurrent_position "
                      << current_position << std::endl;
#endif
            computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
          }
          // ++current_peak;
          ++peak_iter;

        }                        // end while
#ifdef DEBUG_2D
        std::cout << "computed vs experimental signal: " << computed_signal << "\t"
                  << experimental_signal << std::endl;
#endif
        fvec(counter_posf) = step * (computed_signal - experimental_signal);
        ++counter_posf;
      }

    }

    // penalties : especially negative heights have to be penalised
    double penalty = 0.;


    //iterate over all peaks again to compute the penalties
    // first look at all positions and width parameters
    UInt peak = 0, current_peak = 0;
    std::map<Int, std::vector<PeakIndex> >::const_iterator map_iter = matching_peaks.begin();
    for (; map_iter != matching_peaks.end(); ++map_iter)
    {
      std::vector<PeakIndex>::const_iterator vec_iter = map_iter->second.begin();
      double old_position = 0, old_width_l = 0, old_width_r = 0;
      double weight = 0;
      for (; vec_iter != map_iter->second.end(); ++vec_iter)
      {
        // intensity at peak position is stored in the meta data array, not as intensity of the peak
        double old_height = picked_peaks[vec_iter->spectrum].getFloatDataArrays()[1][vec_iter->peak];
        weight += old_height;
        old_position += (vec_iter)->getPeak(picked_peaks).getMZ() * old_height;
        old_width_l += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[3][vec_iter->peak] * old_height;
        old_width_r += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[4][vec_iter->peak] * old_height;

        double p_height = x(peak);
        ++peak;

        if (p_height < 1)
        {
          penalty += 1000000 * penalties.height * pow(fabs(p_height - old_height), 2);
        }

      }
      old_position /= weight;
      old_width_l /= weight;
      old_width_r /= weight;

      double p_position   = x(total_nr_peaks + 3 * current_peak);
      double p_width_l    = x(total_nr_peaks + 3 * current_peak + 1);
      double p_width_r    = x(total_nr_peaks + 3 * current_peak + 2);
      if (p_width_l < 0)
      {
        penalty += 1e7 * penalties.lWidth * pow(fabs(p_width_l - old_width_l), 2);
      }
      else if (p_width_l < 1)
      {
        penalty += 1000 * penalties.lWidth * pow(fabs(p_width_l - old_width_l), 2);
      }
      if (p_width_r < 0)
      {
        penalty += 1e7 * penalties.rWidth * pow(fabs(p_width_r - old_width_r), 2);
      }
      else if (p_width_r < 1)
      {
        penalty += 1000 * penalties.rWidth * pow(fabs(p_width_r - old_width_r), 2);
      }
      if (p_position < 0)
      {
        penalty += 100 * penalties.pos * pow(p_position - old_position, 2);
      }
      if (fabs(old_width_r - p_width_r) > 1)
      {
        penalty += 1000 * penalties.rWidth * pow(old_width_r - p_width_r, 2);
      }
      if (fabs(old_width_l - p_width_l) > 1)
      {
        penalty += 1000 * penalties.lWidth * pow(old_width_l - p_width_l, 2);
      }
      if (fabs(old_position - p_position) > 0.2)
      {
        penalty += 1000 * penalties.pos * pow(p_position - old_position, 2);
      }

      ++current_peak;
    }

    fvec(fvec.size()-1) = penalty;

    return 0;
  }

  // compute Jacobian matrix for the different parameters
  int TwoDOptimization::TwoDOptFunctor::df(const Eigen::VectorXd &x, Eigen::MatrixXd &J)
  {
    // For the conventions on x and params c.f. the commentary in residual()
    //
    // The matrix J is supposed to contain the result when we return from this function.
    // Note: expects the Jacobian as follows:
    //          - each row corresponds to one data point
    //          - each column corresponds to one parameter

    double last_position, step;
    double p_height, p_position, p_width;
    double diff, denom_inv, ddl_left, ddl_right, ddx0, sinh_term;
    Int count = 0;
    Int counter_posf = 0;

    const std::vector<std::pair<SignedSize, SignedSize> > & signal2D = m_data->signal2D;
    std::multimap<double, IsotopeCluster>::iterator iso_map_iter = m_data->iso_map_iter;
    Size total_nr_peaks = m_data->total_nr_peaks;
    const std::map<Int, std::vector<PeakIndex> > & matching_peaks = m_data->matching_peaks;
    std::vector<double> ov_weight(matching_peaks.size(), 0);
    const PeakMap & picked_peaks = m_data->picked_peaks;
    PeakMap::ConstIterator raw_data_first = m_data->raw_data_first;
    const OptimizationFunctions::PenaltyFactorsIntensity & penalties = m_data->penalties;
//          std::vector<double> &positions=static_cast<TwoDOptimization::Data*> (params) ->positions;
//          std::vector<double> &signal=static_cast<TwoDOptimization::Data*> (params) ->signal;
    IsotopeCluster::ChargedIndexSet::iterator peak_iter = iso_map_iter->second.peaks.begin();
    Size num_scans = signal2D.size() / 2;
    //iterate over all scans
    for (Size current_scan = 0; current_scan < num_scans; ++current_scan)
    {
      Size curr_scan_idx = current_scan + iso_map_iter->second.peaks.begin()->first;
      double current_position = ((raw_data_first
                           + signal2D[2 * current_scan].first)->begin()
                          + signal2D[2 * current_scan].second)->getMZ();
      // iterate over all points of the signal
      for (Int current_point = 1;
           current_point +  signal2D[2 * current_scan].second
           <= signal2D[2 * current_scan + 1].second;
           ++current_point)
      {
        last_position = current_position;
        current_position  = ((raw_data_first
                              + signal2D[2 * current_scan].first)->begin()
                             + signal2D[2 * current_scan].second + current_point)->getMZ();
//       double experimental_signal = ((raw_data_first + signal2D[2 * current_scan].first)->begin()
//                               + signal2D[2 * current_scan].second + current_point)->getIntensity();

        step = current_position - last_position;

#ifdef DEBUG_2D
        std::cout << "experimental signal rt " << (raw_data_first
                                                   + signal2D[2 * current_scan].first)->getRT()
                  << "\tmz " << ((raw_data_first
                        + signal2D[2 * current_scan].first)->begin()
                       + signal2D[2 * current_scan].second + current_point)->getMZ()
        << "\tint " << ((raw_data_first
                         + signal2D[2 * current_scan].first)->begin()
                        + signal2D[2 * current_scan].second + current_point)->getIntensity() << std::endl;

#endif

        // size_t current_peak = 0;
        peak_iter = iso_map_iter->second.peaks.begin();
        while (peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first != curr_scan_idx)
        {
          ++peak_iter;
        }
        //iterate over all peaks of the current scan
        while (peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first == curr_scan_idx)
        {
          Int peak_idx = distance(iso_map_iter->second.peaks.begin(), peak_iter);
          double mz_in_hash = ((picked_peaks[peak_iter->first]).begin() + peak_iter->second)->getMZ() * 10;
          std::map<Int, std::vector<PeakIndex> >::const_iterator  m_spec_iter =  matching_peaks.begin();
          Int map_idx = 0;
          while (m_spec_iter->first != (Int)(mz_in_hash + 0.5))
          {
            ++map_idx;
            ++m_spec_iter;
          }
          // if the current peak is in the reference scan take all parameters from the vector x
#ifdef DEBUG_2D
          std::cout << "ref_scan : " << current_peak << "\t "
                    << x(3 * current_peak) << "\t" << x(3 * current_peak + 1)
                    << "\t" << x(3 * current_peak + 2) << std::endl;
#endif
          // Store the current parameters for this peak
          p_position    = x(total_nr_peaks + 3 * map_idx);
          p_height      = x(peak_idx);
          p_width       = (current_position <= p_position) ?
                          x(total_nr_peaks + 3 * map_idx + 1) :
                          x(total_nr_peaks + 3 * map_idx + 2);
          ++count;
          double weight = step * picked_peaks[peak_iter->first].getFloatDataArrays()[1][peak_iter->second];
          ov_weight[map_idx] += weight;
          double ddx0_old = J(counter_posf, total_nr_peaks + 3 * map_idx);
          double ddl_left_old = J(counter_posf, total_nr_peaks + 3 * map_idx + 1);
          double ddl_right_old = J(counter_posf, total_nr_peaks + 3 * map_idx + 2);
          //is it a Lorentz or a Sech - Peak?

          if ((PeakShape::Type)(Int) Math::round((picked_peaks[peak_iter->first]).getFloatDataArrays()[5][peak_iter->second]) == PeakShape::LORENTZ_PEAK)
          {
            diff      = current_position - p_position;
            // partial derivative with respect to the height,...
            denom_inv = 1. / (1. + pow(p_width * diff, 2));
            // left width,...
            ddl_left  = (current_position <= p_position)
                        ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) :
                          0;
            // right width ...
            ddl_right = (current_position  > p_position)
                        ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) :
                          0;

            // and position
            ddx0 = 2 * p_height * pow(p_width, 2) * diff * pow(denom_inv, 2);


            J(counter_posf, total_nr_peaks + 3 * map_idx) = ddx0 * weight  + ddx0_old;
            J(counter_posf, peak_idx) = step * denom_inv;
            J(counter_posf, total_nr_peaks + 3 * map_idx + 1) = ddl_left * weight + ddl_left_old;
            J(counter_posf, total_nr_peaks + 3 * map_idx + 2) = ddl_right * weight + ddl_right_old;
          }
          // if it's a Sech - Peak
          else
          {
            diff      = current_position - p_position;
            denom_inv = 1. / cosh(p_width * diff);
            // The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
            // and can assume that all derivatives vanish
            sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width * diff);


            ddl_left  = (current_position <= p_position)
                        ? -2 * p_height * sinh_term * diff * pow(denom_inv, 3) :
                          0;
            ddl_right = (current_position  > p_position)
                        ? -2 * p_height * sinh_term * diff * pow(denom_inv, 3) :
                          0;

            ddx0      = 2 * p_height * p_width * sinh_term * pow(denom_inv, 3);

            J(counter_posf, total_nr_peaks + 3 * map_idx) = ddx0 * weight + ddx0_old;
            J(counter_posf, peak_idx) = step * pow(denom_inv, 2);
            J(counter_posf, total_nr_peaks + 3 * map_idx  + 1) = ddl_left * weight + ddl_left_old;
            J(counter_posf, total_nr_peaks + 3 * map_idx  + 2) = ddl_right * weight + ddl_right_old;
          }
          // ++current_peak;
          ++peak_iter;

        }                        // end while

        ++counter_posf;
      }

    }

    for (Size cluster = 0; cluster < matching_peaks.size(); ++cluster)
    {
      for (int j = 0; j < J.rows() - 1; ++j)
      {
        J(j, total_nr_peaks + 3 * cluster)
              = J(j, total_nr_peaks + 3 * cluster) / ov_weight[cluster];
        J(j, total_nr_peaks + 3 * cluster + 1)
              = J(j, total_nr_peaks + 3 * cluster + 1) / ov_weight[cluster];
        J(j, total_nr_peaks + 3 * cluster + 2)
              = J(j, total_nr_peaks + 3 * cluster + 2) / ov_weight[cluster];
      }
    }

    //iterate over all peaks again to compute the penalties
    // first look at all positions and width parameters
    UInt peak = 0, current_peak = 0;
    std::map<Int, std::vector<PeakIndex> >::const_iterator map_iter = matching_peaks.begin();

    for (; map_iter != matching_peaks.end(); ++map_iter)
    {
      std::vector<PeakIndex>::const_iterator vec_iter   = map_iter->second.begin();
      double old_position = 0, old_width_l = 0, old_width_r = 0;
      double weight = 0;
      double penalty_h = 0, penalty_l = 0, penalty_r = 0, penalty_p = 0;
      for (; vec_iter != map_iter->second.end(); ++vec_iter)
      {
        double old_height = picked_peaks[vec_iter->spectrum].getFloatDataArrays()[1][vec_iter->peak];
        weight += old_height;
        old_position += (vec_iter)->getPeak(picked_peaks).getMZ() * old_height;
        old_width_l += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[3][vec_iter->peak] * old_height;
        old_width_r += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[4][vec_iter->peak] * old_height;

        double p_height     = x(peak);


        double penalty_height = 2. * penalties.height * fabs(p_height - old_height);
        if (p_height < 1)
        {
          penalty_h += 1000000 * penalty_height;
        }
        J(counter_posf, peak) = penalty_h;
        ++peak;
      }
      old_position /= weight;
      old_width_l /= weight;
      old_width_r /= weight;

      // std::cout << old_position << "vs. ";
      double p_position   = x(total_nr_peaks + 3 * current_peak);
      double p_width_l    = x(total_nr_peaks + 3 * current_peak + 1);
      double p_width_r    = x(total_nr_peaks + 3 * current_peak + 2);
      double penalty_lwidth = 2. * penalties.lWidth * fabs(p_width_l - old_width_l);
      double penalty_rwidth = 2. * penalties.rWidth * fabs(p_width_r - old_width_r);
      double penalty_pos    = 2. * penalties.pos * fabs(p_position - old_position);
      //std::cout << p_position<<std::endl;
#ifdef DEBUG_2D
      std::cout << "penalty_lwidth " << penalty_lwidth << "penalty_rwidth " << penalty_rwidth
                << "penalty_pos " << penalty_pos << std::endl;
#endif
      if (p_width_l < 0)
      {
        penalty_l += 1e7 * penalty_lwidth;
      }
      else if (p_width_l < 1)
      {
        penalty_l += 2000 * penalties.lWidth * (fabs(p_width_l - old_width_l));
      }
      if (p_width_r < 0)
      {
        penalty_r += 1e7 * penalty_rwidth;
      }
      else if (p_width_r < 1)
      {
        penalty_r += 2000 * penalties.rWidth * (fabs(p_width_r - old_width_r));
      }
      if (p_position < 0)
      {
        penalty_p += 200 * penalty_pos;
      }

      if (fabs(old_position - p_position) > 0.2)
      {
        penalty_p += 2000 * penalties.pos * fabs(p_position - old_position);
      }
      if (fabs(old_width_r - p_width_r) > 1)
      {
        penalty_r += 1000 * penalty_rwidth;
      }
      if (fabs(old_width_l - p_width_l) > 1)
      {
        penalty_l += 1000 * penalty_lwidth;
      }

      J(counter_posf, total_nr_peaks + 3 * current_peak + 1) = penalty_l;
      J(counter_posf, total_nr_peaks + 3 * current_peak + 2) = penalty_r;
      J(counter_posf, total_nr_peaks + 3 * current_peak) = penalty_p;

      ++current_peak;
    }
    return 0;
  }

  void TwoDOptimization::optimize(InputSpectrumIterator first, InputSpectrumIterator last, PeakMap& ms_exp, bool real2D)
  {
    //#define DEBUG_2D
    //check if the input maps have the same number of spectra
    if ((UInt)distance(first, last) != ms_exp.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error in Two2Optimization: Raw and peak map do not have the same number of spectra");
    }
    //do nothing if there are no scans
    if (ms_exp.empty())
    {
      return;
    }
    //check if required meta data arrays are present (for each scan)
    for (Size i = 0; i < ms_exp.size(); ++i)
    {
      //check if enough meta data arrays are present
      if (ms_exp[i].getFloatDataArrays().size() < 6)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error in Two2Optimization: Not enough meta data arrays present (1:area, 5:shape, 3:left width, 4:right width)");
      }
      bool area = ms_exp[i].getFloatDataArrays()[1].getName() == "maximumIntensity";
      bool wleft = ms_exp[i].getFloatDataArrays()[3].getName() == "leftWidth";
      bool wright = ms_exp[i].getFloatDataArrays()[4].getName() == "rightWidth";
      bool shape = ms_exp[i].getFloatDataArrays()[5].getName() == "peakShape";

      if (!area || !wleft || !wright || !shape)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error in Two2Optimization: One or several meta data arrays missing (1:intensity, 5:shape, 3:left width, 4:right width)");
      }
    }
    real_2D_ = real2D;
    typedef MSSpectrum SpectrumType;

    typename PeakMap::Iterator ms_exp_it = ms_exp.begin();
    typename PeakMap::Iterator ms_exp_it_end = ms_exp.end();
    if (ms_exp.empty())
    {
      return;
    }
    // stores the monoisotopic peaks of isotopic clusters
    std::vector<double> iso_last_scan;
    std::vector<double> iso_curr_scan;
    std::vector<std::multimap<double, IsotopeCluster>::iterator> clusters_last_scan;
    std::vector<std::multimap<double, IsotopeCluster>::iterator> clusters_curr_scan;
    std::multimap<double, IsotopeCluster>::iterator cluster_iter;
    double current_rt = ms_exp_it->getRT(), last_rt  = 0;

    // retrieve values for accepted peaks distances
    max_peak_distance_ = param_.getValue("2d:max_peak_distance");
    double tolerance_mz = param_.getValue("2d:tolerance_mz");

    UInt current_charge     = 0; // charge state of the current isotopic cluster
    double mz_in_hash   = 0; // used as reference to the current isotopic peak

    // sweep through scans
    for (UInt curr_scan = 0; ms_exp_it + curr_scan != ms_exp_it_end; ++curr_scan)
    {
      Size nr_peaks_in_scan = (ms_exp_it + curr_scan)->size();
      if (nr_peaks_in_scan == 0)
        continue;

      //last_rt = current_rt;
      current_rt = (ms_exp_it + curr_scan)->getRT();
      typename PeakMap::SpectrumType::Iterator peak_it  = (ms_exp_it + curr_scan)->begin();

      // copy cluster information of least scan
      iso_last_scan = iso_curr_scan;
      iso_curr_scan.clear();
      clusters_last_scan = clusters_curr_scan;
      clusters_curr_scan.clear();

#ifdef DEBUG_2D
      std::cout << "Next scan with rt: " << current_rt << std::endl;
      std::cout << "Next scan, rt = " << current_rt << " last_rt: " << last_rt << std::endl;
      std::cout << "---------------------------------------------------------------------------" << std::endl;
#endif
      MSSpectrum s;
      s.setRT(current_rt);
      // check if there were scans in between
      if (last_rt == 0 || // are we in the first scan
          ((lower_bound(first, last, s, typename SpectrumType::RTLess()) - 1)->getRT() == last_rt))
      {


        for (UInt curr_peak = 0; curr_peak < (ms_exp_it + curr_scan)->size() - 1; ++curr_peak)
        {

          // store the m/z of the current peak
          double curr_mz         = (peak_it + curr_peak)->getMZ();
          double dist2nextpeak = (peak_it + curr_peak + 1)->getMZ() - curr_mz;

          if (dist2nextpeak <= max_peak_distance_) // one single peak without neighbors isn't optimized
          {
#ifdef DEBUG_2D
            std::cout << "Isotopic pattern found ! " << std::endl;
            std::cout << "We are at: " << (peak_it + curr_peak)->getMZ()  << " " << curr_mz << std::endl;
#endif
            if (!iso_last_scan.empty()) // Did we find any isotopic cluster in the last scan?
            {
              std::sort(iso_last_scan.begin(), iso_last_scan.end());
              // there were some isotopic clusters in the last scan...
              std::vector<double>::iterator it =
                searchInScan_(iso_last_scan.begin(), iso_last_scan.end(), curr_mz);

              double delta_mz = fabs(*it - curr_mz);
              //std::cout << delta_mz << " "<< tolerance_mz << std::endl;
              if (delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
              {
                mz_in_hash = curr_mz; // update current hash key

                // create new isotopic cluster
// #ifdef DEBUG_2D
//                                                      std::cout << "Last peak cluster too far, creating new cluster at "<<curr_mz << std::endl;
// #endif
                IsotopeCluster new_cluster;
                new_cluster.peaks.charge  = current_charge;
                new_cluster.scans.push_back(curr_scan);
                cluster_iter = iso_map_.insert(std::pair<double, IsotopeCluster>(mz_in_hash, new_cluster));

              }
              else
              {
// //#ifdef DEBUG_2D
//                                                      std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
// //#endif
                cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(), it)];

                // check whether this scan is already contained
                if (find(cluster_iter->second.scans.begin(), cluster_iter->second.scans.end(), curr_scan)
                    == cluster_iter->second.scans.end())
                {
                  cluster_iter->second.scans.push_back(curr_scan);
                }

//                                                      //#ifdef DEBUG_2D
//                                                      std::cout << "Cluster with " << cluster_iter->second.peaks.size()
//                                                                          << " peaks retrieved." << std::endl;
//                                                      //#endif
              }

            }
            else                             // last scan did not contain any isotopic cluster
            {
//                                              //#ifdef DEBUG_2D
//                                              std::cout << "Last scan was empty => creating new cluster." << std::endl;
//                                              std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
//                                              //#endif

              mz_in_hash = curr_mz; // update current hash key

              // create new isotopic cluster
              IsotopeCluster new_cluster;
              new_cluster.peaks.charge  = current_charge;
              new_cluster.scans.push_back(curr_scan);
              cluster_iter = iso_map_.insert(std::pair<double, IsotopeCluster>(mz_in_hash, new_cluster));

            }

//                                      //#ifdef DEBUG_2D
//                                      std::cout << "Storing found peak in current isotopic cluster" << std::endl;
//                                      //#endif



            cluster_iter->second.peaks.insert(std::pair<UInt, UInt>(curr_scan, curr_peak));

            iso_curr_scan.push_back(mz_in_hash);
            clusters_curr_scan.push_back(cluster_iter);
            ++curr_peak;

            cluster_iter->second.peaks.insert(std::pair<UInt, UInt>(curr_scan, curr_peak));
            iso_curr_scan.push_back((peak_it + curr_peak)->getMZ());
            clusters_curr_scan.push_back(cluster_iter);

            // check distance to next peak
            if ((curr_peak + 1) >= nr_peaks_in_scan)
              break;
            dist2nextpeak = (peak_it + curr_peak + 1)->getMZ() -  (peak_it + curr_peak)->getMZ();


            // loop until end of isotopic pattern in this scan
            while (dist2nextpeak <= max_peak_distance_
                  &&  curr_peak < (nr_peaks_in_scan - 1))
            {
              cluster_iter->second.peaks.insert(std::pair<UInt, UInt>(curr_scan, curr_peak + 1)); // save peak in cluster
              iso_curr_scan.push_back((peak_it + curr_peak + 1)->getMZ());
              clusters_curr_scan.push_back(cluster_iter);
              // std::cout << "new entered: "<<(peak_it+curr_peak+1)->getMZ()<<" im while"<<std::endl;
              ++curr_peak;
              if (curr_peak >= nr_peaks_in_scan - 1)
                break;
              dist2nextpeak = (peak_it + curr_peak + 1)->getMZ() -  (peak_it + curr_peak)->getMZ(); // get distance to next peak


            } // end while(...)



          } // end of if (dist2nextpeak <= max_peak_distance_)
          else
          {
            if (!iso_last_scan.empty()) // Did we find any isotopic cluster in the last scan?
            {
              std::sort(iso_last_scan.begin(), iso_last_scan.end());
              // there were some isotopic clusters in the last scan...
              std::vector<double>::iterator it =
                searchInScan_(iso_last_scan.begin(), iso_last_scan.end(), curr_mz);

              double delta_mz = fabs(*it - curr_mz);
              // std::cout << delta_mz << " "<< tolerance_mz << std::endl;
              if (delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
              {
                mz_in_hash = curr_mz; // update current hash key

                // create new isotopic cluster
//                                                      //#ifdef DEBUG_2D
//                                                      std::cout << "Last peak cluster too far, creating new cluster at "<<curr_mz << std::endl;
//                                                      //#endif
                IsotopeCluster new_cluster;
                new_cluster.peaks.charge  = current_charge;
                new_cluster.scans.push_back(curr_scan);
                cluster_iter = iso_map_.insert(std::pair<double, IsotopeCluster>(mz_in_hash, new_cluster));

              }
              else
              {
//                                                      //#ifdef DEBUG_2D
//                                                      std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
//                                                      //#endif
                cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(), it)];

                // check whether this scan is already contained
                if (find(cluster_iter->second.scans.begin(), cluster_iter->second.scans.end(), curr_scan)
                    == cluster_iter->second.scans.end())
                {
                  cluster_iter->second.scans.push_back(curr_scan);
                }

//                                                      //#ifdef DEBUG_2D
//                                                      std::cout << "Cluster with " << cluster_iter->second.peaks.size()
//                                                                          << " peaks retrieved." << std::endl;
//                                                      //#endif
              }

            }
            else                             // last scan did not contain any isotopic cluster
            {
//                                              //#ifdef DEBUG_2D
//                                              std::cout << "Last scan was empty => creating new cluster." << std::endl;
//                                              std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
//                                              //#endif

              mz_in_hash = curr_mz; // update current hash key

              // create new isotopic cluster
              IsotopeCluster new_cluster;
              new_cluster.peaks.charge  = current_charge;
              new_cluster.scans.push_back(curr_scan);
              cluster_iter = iso_map_.insert(std::pair<double, IsotopeCluster>(mz_in_hash, new_cluster));

            }

//                                      //#ifdef DEBUG_2D
//                                      std::cout << "Storing found peak in current isotopic cluster" << std::endl;
//                                      //#endif



            cluster_iter->second.peaks.insert(std::pair<UInt, UInt>(curr_scan, curr_peak));

            iso_curr_scan.push_back(mz_in_hash);
            clusters_curr_scan.push_back(cluster_iter);


          }

          current_charge = 0; // reset charge
        } // end for (...)
      }
      last_rt = current_rt;
    }
    curr_region_ = iso_map_.begin();
#ifdef DEBUG_2D
    std::cout << iso_map_.size() << " isotopic clusters were found ! " << std::endl;
#endif

    if (real_2D_)
      optimizeRegions_(first, last, ms_exp);
    else
      optimizeRegionsScanwise_(first, last, ms_exp);
    //#undef DEBUG_2D
  }


  void TwoDOptimization::optimizeRegions_(InputSpectrumIterator& first,
                                          InputSpectrumIterator& last,
                                          PeakMap& ms_exp)
  {
    Int counter = 0;
    // go through the clusters
    for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
         it != iso_map_.end();
         ++it)
    {
#ifdef DEBUG_2D
      std::cout << "element: " << counter << std::endl;
      std::cout << "mz: " << it->first << std::endl << "rts: ";
//              for(Size i=0;i<it->second.scans.size();++i) std::cout << it->second.scans[i] << "\n";
      std::cout << std::endl << "peaks: ";
      IsotopeCluster::IndexSet::const_iterator iter = it->second.peaks.begin();
      for (; iter != it->second.peaks.end(); ++iter)
        std::cout << ms_exp[iter->first].getRT() << " " << (ms_exp[iter->first][iter->second]).getMZ() << std::endl;

//for(Size i=0;i<it->second.peaks.size();++i) std::cout << ms_exp[it->first].getRT() << " "<<(ms_exp[it->first][it->second]).getMZ()<<std::endl;
      std::cout << std::endl << std::endl;

#endif

      // prepare for optimization:
      // determine the matching peaks
      matching_peaks_.clear();
      findMatchingPeaks_(it, ms_exp);
      TwoDOptimization::Data twoD_data;
      twoD_data.penalties = penalties_;
      twoD_data.matching_peaks = matching_peaks_;
      // and the endpoints of each isotope pattern in the cluster
      getRegionEndpoints_(ms_exp, first, last, counter, 400, twoD_data);

      // peaks have to be stored globally
      twoD_data.iso_map_iter = it;

      twoD_data.picked_peaks = ms_exp;
      twoD_data.raw_data_first =  first;

      Size nr_diff_peaks = matching_peaks_.size();
      twoD_data.total_nr_peaks = it->second.peaks.size();

      Size nr_parameters = nr_diff_peaks * 3 + twoD_data.total_nr_peaks;

      // initialize and set parameters for optimization
      Eigen::VectorXd x_init (nr_parameters);
      x_init.setZero();

      std::map<Int, std::vector<PeakIndex> >::iterator m_peaks_it = twoD_data.matching_peaks.begin();
      Int peak_counter = 0;
      Int diff_peak_counter = 0;
      // go through the matching peaks
      for (; m_peaks_it != twoD_data.matching_peaks.end(); ++m_peaks_it)
      {
        double av_mz = 0, av_lw = 0, av_rw = 0, avr_height = 0, height;
        std::vector<PeakIndex>::iterator iter_iter = (m_peaks_it)->second.begin();
        for (; iter_iter != m_peaks_it->second.end(); ++iter_iter)
        {
          height = ms_exp[(iter_iter)->spectrum].getFloatDataArrays()[1][(iter_iter)->peak]; //(iter_iter)->getPeak(ms_exp).getIntensity();
          avr_height += height;
          av_mz += (iter_iter)->getPeak(ms_exp).getMZ() * height;
          av_lw += ms_exp[(iter_iter)->spectrum].getFloatDataArrays()[3][(iter_iter)->peak] * height; //left width
          av_rw +=    ms_exp[(iter_iter)->spectrum].getFloatDataArrays()[4][(iter_iter)->peak] * height; //right width
          x_init(peak_counter) = height;
          ++peak_counter;
        }
        x_init(twoD_data.total_nr_peaks + 3 * diff_peak_counter) = av_mz / avr_height;
        x_init(twoD_data.total_nr_peaks + 3 * diff_peak_counter + 1) = av_lw / avr_height;
        x_init(twoD_data.total_nr_peaks + 3 * diff_peak_counter + 2) = av_rw / avr_height;
        ++diff_peak_counter;
      }

#ifdef DEBUG_2D
      std::cout << "----------------------------\n\nstart_value: " << std::endl;
      for (Size k = 0; k < start_value->size; ++k)
      {
          std::cout << x_init(k) << std::endl;
      }
#endif
      Int num_positions = 0;
      for (Size i = 0; i < twoD_data.signal2D.size(); i += 2)
      {
        num_positions += (twoD_data.signal2D[i + 1].second - twoD_data.signal2D[i].second + 1);
#ifdef DEBUG_2D
        std::cout << twoD_data.signal2D[i + 1].second << " - " << twoD_data.signal2D[i].second << " +1 " << std::endl;
#endif

      }
#ifdef DEBUG_2D
      std::cout << "num_positions : " << num_positions << std::endl;
#endif
      
      TwoDOptFunctor functor (nr_parameters, std::max(num_positions + 1, (Int)(nr_parameters)), &twoD_data);
      Eigen::LevenbergMarquardt<TwoDOptFunctor> lmSolver (functor);
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-TwoDOptimization:", "Could not fit the data: Error " + String(status));
      }

      Int peak_idx = 0;
      std::map<Int, std::vector<PeakIndex> >::iterator itv
        = twoD_data.matching_peaks.begin();
      for (; itv != twoD_data.matching_peaks.end(); ++itv)
      {
        Int i = distance(twoD_data.matching_peaks.begin(), itv);
        for (Size j = 0; j < itv->second.size(); ++j)
        {

#ifdef DEBUG_2D
          std::cout << "pos: " << itv->second[j].getPeak(ms_exp).getMZ() << "\nint: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[1][itv->second[j].peak] //itv->second[j].getPeak(ms_exp).getIntensity()
                    << "\nlw: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[3][itv->second[j].peak]
                    << "\nrw: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[4][itv->second[j].peak] << "\n";

#endif
          double mz = x_init(twoD_data.total_nr_peaks + 3 * i);
          ms_exp[itv->second[j].spectrum][itv->second[j].peak].setMZ(mz);
          double height = x_init(peak_idx);
          ms_exp[itv->second[j].spectrum].getFloatDataArrays()[1][itv->second[j].peak] = height;
          double left_width = x_init(twoD_data.total_nr_peaks + 3 * i + 1);
          ms_exp[itv->second[j].spectrum].getFloatDataArrays()[3][itv->second[j].peak] = left_width;
          double right_width = x_init(twoD_data.total_nr_peaks + 3 * i + 2);

          ms_exp[itv->second[j].spectrum].getFloatDataArrays()[4][itv->second[j].peak] = right_width;
          // calculate area
          if ((PeakShape::Type)(Int)ms_exp[itv->second[j].spectrum].getFloatDataArrays()[5][itv->second[j].peak] == PeakShape::LORENTZ_PEAK)
          {
            double x_left_endpoint = mz - 1 / left_width* sqrt(height / 1 - 1);
            double x_right_endpoint = mz + 1 / right_width* sqrt(height / 1 - 1);
            double area_left = -height / left_width* atan(left_width * (x_left_endpoint - mz));
            double area_right = -height / right_width* atan(right_width * (mz - x_right_endpoint));
            ms_exp[itv->second[j].spectrum][itv->second[j].peak].setIntensity(area_left + area_right);
          }
          else         // it's a sech peak
          {
            double x_left_endpoint = mz - 1 / left_width* std::acosh(sqrt(height / 0.001));
            double x_right_endpoint = mz + 1 / right_width* std::acosh(sqrt(height / 0.001));
            double area_left = -height / left_width * (sinh(left_width * (mz - x_left_endpoint)) / cosh(left_width * (mz - x_left_endpoint)));
            double area_right = -height / right_width * (sinh(right_width * (mz - x_right_endpoint)) / cosh(right_width * (mz - x_right_endpoint)));
            ms_exp[itv->second[j].spectrum][itv->second[j].peak].setIntensity(area_left + area_right);
          }


#ifdef DEBUG_2D
          std::cout << "pos: " << itv->second[j].getPeak(ms_exp).getMZ() << "\nint: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[1][itv->second[j].peak] //itv->second[j].getPeak(ms_exp).getIntensity()
                    << "\nlw: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[3][itv->second[j].peak]
                    << "\nrw: " << itv->second[j].getSpectrum(ms_exp).getFloatDataArrays()[4][itv->second[j].peak] << "\n";
#endif

          ++peak_idx;


        }
      }

      ++counter;
    } // end for
    //#undef DEBUG_2D
  }

  void TwoDOptimization::optimizeRegionsScanwise_(InputSpectrumIterator& first,
                                                  InputSpectrumIterator& last,
                                                  PeakMap& ms_exp)
  {
    Int counter = 0;
    TwoDOptimization::Data d;
    d.picked_peaks = ms_exp;
    d.raw_data_first =  first;

    struct OpenMS::OptimizationFunctions::PenaltyFactors penalties;

    ParamValue pv = param_.getValue("penalties:position");
    if (pv.isEmpty() || pv.toString().empty())
      penalties.pos = 0.;
    else
      penalties.pos = (float)pv;

    pv = param_.getValue("penalties:left_width");
    if (pv.isEmpty() || pv.toString().empty())
      penalties.lWidth = 1.;
    else
      penalties.lWidth = (float)pv;

    pv = param_.getValue("penalties:right_width");
    if (pv.isEmpty() || pv.toString().empty())
      penalties.rWidth = 1.;
    else
      penalties.rWidth = (float)pv;
#ifdef DEBUG_2D
    std::cout << penalties.pos << " "
              << penalties.rWidth << " "
              << penalties.lWidth << std::endl;
#endif
//      PeakMap::const_iterator help = first;
//      // std::cout << "\n\n\n\n---------------------------------------------------------------";
//      while(help!=last)
//          {
//              // std::cout<<help->getRT()<<std::endl;
//              ++help;
//          }
    // std::cout << "---------------------------------------------------------------\n\n\n\n";

    UInt max_iteration;
    pv = param_.getValue("iterations");
    if (pv.isEmpty() || pv.toString().empty())
      max_iteration = 15;
    else
      max_iteration = (UInt)pv;

    std::vector<PeakShape> peak_shapes;


    // go through the clusters
    for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
         it != iso_map_.end();
         ++it)
    {
      d.iso_map_iter = it;
#ifdef DEBUG_2D
      std::cerr << "element: " << counter << std::endl;
      std::cerr << "mz: " << it->first << std::endl << "rts: ";
      for (Size i = 0; i < it->second.scans.size(); ++i)
        std::cerr << it->second.scans[i] << "\n";
      std::cerr << std::endl << "peaks: ";
      IsotopeCluster::IndexSet::const_iterator iter = it->second.peaks.begin();
      for (; iter != it->second.peaks.end(); ++iter)
        std::cerr << ms_exp[iter->first].getRT() << " " << (ms_exp[iter->first][iter->second]).getMZ() << std::endl;
      //for(Size i=0;i<it->second.peaks_.size();++i) std::cout << ms_exp[it->first].getRT() << " "<<(ms_exp[it->first][it->second]).getMZ()<<std::endl;
      std::cerr << std::endl << std::endl;

#endif
      // prepare for optimization:
      // determine the matching peaks
      // and the endpoints of each isotope pattern in the cluster

      getRegionEndpoints_(ms_exp, first, last, counter, 400, d);
      OptimizePick::Data data;


      Size idx = 0;
      for (Size i = 0; i < d.signal2D.size() / 2; ++i)
      {
        data.positions.clear();
        data.signal.clear();

        PeakMap::SpectrumType::const_iterator ms_it =
          (d.raw_data_first + d.signal2D[2 * i].first)->begin() + d.signal2D[2 * i].second;
        Int size = distance(ms_it, (d.raw_data_first + d.signal2D[2 * i].first)->begin() + d.signal2D[2 * i + 1].second);
        data.positions.reserve(size);
        data.signal.reserve(size);

        while (ms_it != (d.raw_data_first + d.signal2D[2 * i].first)->begin() + d.signal2D[2 * i + 1].second)
        {
          data.positions.push_back(ms_it->getMZ());
          data.signal.push_back(ms_it->getIntensity());
          ++ms_it;
        }


        IsotopeCluster::IndexPair pair;
        pair.first =  d.iso_map_iter->second.peaks.begin()->first + idx;

        IsotopeCluster::IndexSet::const_iterator set_iter = lower_bound(d.iso_map_iter->second.peaks.begin(),
                                                                        d.iso_map_iter->second.peaks.end(),
                                                                        pair, [](auto& left, auto& right){return left.first < right.first;});


        // find the last entry with this rt-value
        ++pair.first;
        IsotopeCluster::IndexSet::const_iterator set_iter2 = lower_bound(d.iso_map_iter->second.peaks.begin(),
                                                                         d.iso_map_iter->second.peaks.end(),
                                                                         pair, [](auto& left, auto& right){return left.first < right.first;});

        while (set_iter != set_iter2)
        {
          const Size peak_index = set_iter->second;
          const MSSpectrum& spec = ms_exp[set_iter->first];
          PeakShape shape(spec.getFloatDataArrays()[1][peak_index], //intensity
                          spec[peak_index].getMZ(),
                          spec.getFloatDataArrays()[3][peak_index], //left width
                          spec.getFloatDataArrays()[4][peak_index], //right width
                          spec[peak_index].getIntensity(), //area is stored in peak intensity
                          PeakShape::Type(Int(spec.getFloatDataArrays()[5][peak_index]))); //shape
          peak_shapes.push_back(shape);
          ++set_iter;
        }
#ifdef DEBUG_2D
        std::cout << "rt "
                  << (d.raw_data_first + d.signal2D[2 * i].first)->getRT()
                  << "\n";
#endif
        OptimizePick opt(penalties, max_iteration);
#ifdef DEBUG_2D
        std::cout << "vorher\n";

        for (Size p = 0; p < peak_shapes.size(); ++p)
        {
          std::cout << peak_shapes[p].mz_position << "\t" << peak_shapes[p].height
                    << "\t" << peak_shapes[p].left_width << "\t" << peak_shapes[p].right_width  << std::endl;
        }
#endif
        opt.optimize(peak_shapes, data);
#ifdef DEBUG_2D
        std::cout << "nachher\n";
        for (Size p = 0; p < peak_shapes.size(); ++p)
        {
          std::cout << peak_shapes[p].mz_position << "\t" << peak_shapes[p].height
                    << "\t" << peak_shapes[p].left_width << "\t" << peak_shapes[p].right_width  << std::endl;
        }
#endif
        std::sort(peak_shapes.begin(), peak_shapes.end(), PeakShape::PositionLess());
        pair.first =  d.iso_map_iter->second.peaks.begin()->first + idx;

        set_iter = lower_bound(d.iso_map_iter->second.peaks.begin(),
                               d.iso_map_iter->second.peaks.end(),
                               pair, [](auto& left, auto& right){return left.first < right.first;});
        Size p = 0;
        while (p < peak_shapes.size())
        {
          MSSpectrum& spec = ms_exp[set_iter->first];
          spec[set_iter->second].setMZ(peak_shapes[p].mz_position);
          spec.getFloatDataArrays()[3][set_iter->second] = peak_shapes[p].left_width;
          spec.getFloatDataArrays()[4][set_iter->second] = peak_shapes[p].right_width;
          spec.getFloatDataArrays()[1][set_iter->second] = peak_shapes[p].height; // maximum intensity
          // calculate area
          if (peak_shapes[p].type == PeakShape::LORENTZ_PEAK)
          {
            PeakShape& ps = peak_shapes[p];
            double x_left_endpoint = ps.mz_position - 1 / ps.left_width* sqrt(ps.height / 1 - 1);
            double x_right_endpoint = ps.mz_position + 1 / ps.right_width* sqrt(ps.height / 1 - 1);
            double area_left = -ps.height / ps.left_width* atan(ps.left_width * (x_left_endpoint - ps.mz_position));
            double area_right = -ps.height / ps.right_width* atan(ps.right_width * (ps.mz_position - x_right_endpoint));
            spec[set_iter->second].setIntensity(area_left + area_right); // area is stored as peak intensity
          }
          else        //It's a Sech - Peak
          {
            PeakShape& ps = peak_shapes[p];
            double x_left_endpoint = ps.mz_position - 1 / ps.left_width* std::acosh(sqrt(ps.height / 0.001));
            double x_right_endpoint = ps.mz_position + 1 / ps.right_width* std::acosh(sqrt(ps.height / 0.001));
            double area_left = ps.height / ps.left_width * (sinh(ps.left_width * (ps.mz_position - x_left_endpoint)) / cosh(ps.left_width * (ps.mz_position - x_left_endpoint)));
            double area_right = -ps.height / ps.right_width * (sinh(ps.right_width * (ps.mz_position - x_right_endpoint)) / cosh(ps.right_width * (ps.mz_position - x_right_endpoint)));
            spec[set_iter->second].setIntensity(area_left + area_right); // area is stored as peak intensity
          }
          ++set_iter;
          ++p;
        }
        ++idx;
        peak_shapes.clear();
      }

      ++counter;
    }
  }

  void TwoDOptimization::getRegionEndpoints_(PeakMap& exp,
                                             InputSpectrumIterator& first,
                                             InputSpectrumIterator& last,
                                             Size iso_map_idx,
                                             double noise_level,
                                             TwoDOptimization::Data& d)
  {
    d.signal2D.clear();
    typedef typename InputSpectrumIterator::value_type InputExperimentType;
    typedef typename InputExperimentType::value_type InputPeakType;
    typedef std::multimap<double, IsotopeCluster> MapType;

    double rt, first_peak_mz, last_peak_mz;

    typename PeakMap::SpectrumType spec;
    InputPeakType peak;

    MapType::iterator iso_map_iter = iso_map_.begin();
    for (Size i = 0; i < iso_map_idx; ++i)
      ++iso_map_iter;

#ifdef DEBUG2D
    std::cout << "rt begin: " << exp[iso_map_iter->second.scans[0]].getRT()
              << "\trt end: " << exp[iso_map_iter->second.scans[iso_map_iter->second.scans.size() - 1]].getRT()
              << " \t" << iso_map_iter->second.scans.size() << " scans"
              << std::endl;
#endif

    // get left and right endpoint for all scans in the current cluster
    for (Size i = 0; i < iso_map_iter->second.scans.size(); ++i)
    {
      typename PeakMap::iterator exp_it;

      // first the right scan through binary search
      rt = exp[iso_map_iter->second.scans[i]].getRT();
      spec.setRT(rt);
      InputSpectrumIterator iter = lower_bound(first, last, spec, MSSpectrum::RTLess());
      //  if(iter->getRT() != rt) --iter;
      exp_it = exp.RTBegin(rt);
#ifdef DEBUG2D
      std::cout << exp_it->getRT() << " vs " << iter->getRT() << std::endl;
#endif
      // now the right mz
      IsotopeCluster::IndexPair pair;
      pair.first =  iso_map_iter->second.peaks.begin()->first + i;
      // get iterator in peaks-set that points to the first peak in the current scan
      IsotopeCluster::IndexSet::const_iterator set_iter = lower_bound(iso_map_iter->second.peaks.begin(),
                                                                      iso_map_iter->second.peaks.end(),
                                                                      pair, [](auto& left, auto& right){return left.first < right.first;});

      // consider a bit more of the signal to the left
      first_peak_mz = (exp_it->begin() + set_iter->second)->getMZ() - 1;

      // find the last entry with this rt-value
      ++pair.first;
      IsotopeCluster::IndexSet::const_iterator set_iter2 = lower_bound(iso_map_iter->second.peaks.begin(),
                                                                       iso_map_iter->second.peaks.end(),
                                                                       pair, [](auto& left, auto& right){return left.first < right.first;});

      if (i == iso_map_iter->second.scans.size() - 1)
      {
        set_iter2 = iso_map_iter->second.peaks.end();
        --set_iter2;
      }
      else if (set_iter2 != iso_map_iter->second.peaks.begin())
        --set_iter2;

      last_peak_mz = (exp_it->begin() + set_iter2->second)->getMZ() + 1;

      //std::cout << rt<<": first peak mz "<<first_peak_mz << "\tlast peak mz "<<last_peak_mz <<std::endl;
      peak.setPosition(first_peak_mz);
      typename PeakMap::SpectrumType::const_iterator raw_data_iter
        = lower_bound(iter->begin(), iter->end(), peak, typename InputPeakType::PositionLess());
      if (raw_data_iter != iter->begin())
      {
        --raw_data_iter;
      }
      double intensity = raw_data_iter->getIntensity();
      // while the intensity is falling go to the left
      while (raw_data_iter != iter->begin() && (raw_data_iter - 1)->getIntensity() < intensity &&
             (raw_data_iter - 1)->getIntensity() > noise_level)
      {
        --raw_data_iter;
        intensity = raw_data_iter->getIntensity();
      }
      ++raw_data_iter;
      IsotopeCluster::IndexPair left, right;
      left.first = distance(first, iter);
      left.second = raw_data_iter - iter->begin();
#ifdef DEBUG2D
      std::cout << "left: " << iter->getRT() << "\t" << raw_data_iter->getMZ() << std::endl;
#endif
      // consider a bit more of the signal to the right
      peak.setPosition(last_peak_mz + 1);
      raw_data_iter
        = upper_bound(iter->begin(), iter->end(), peak, typename InputPeakType::PositionLess());
      if (raw_data_iter == iter->end())
        --raw_data_iter;
      intensity = raw_data_iter->getIntensity();
      // while the intensity is falling go to the right
      while (raw_data_iter + 1 != iter->end() && (raw_data_iter + 1)->getIntensity() < intensity)
      {
        ++raw_data_iter;
        intensity = raw_data_iter->getIntensity();
        if ((raw_data_iter + 1 != iter->end()) && (raw_data_iter + 1)->getIntensity() > noise_level)
          break;
      }
      right.first = left.first;
      right.second = raw_data_iter - iter->begin();
#ifdef DEBUG2D
      std::cout << "right: " << iter->getRT() << "\t" << raw_data_iter->getMZ() << std::endl;
#endif
      // region endpoints are stored in global vector
      d.signal2D.push_back(left);
      d.signal2D.push_back(right);
    }
#ifdef DEBUG2D
    std::cout << first_peak_mz << "\t" << last_peak_mz << std::endl;
#endif
  }


}

