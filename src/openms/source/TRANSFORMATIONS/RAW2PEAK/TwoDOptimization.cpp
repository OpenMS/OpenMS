// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <map>
#include <Eigen/Core>

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
    defaults_.setValue("2d:tolerance_mz", 2.2, "mz tolerance for cluster construction", ListUtils::create<String>("advanced"));
    defaults_.setValue("2d:max_peak_distance", 1.2, "maximal peak distance in mz in a cluster", ListUtils::create<String>("advanced"));
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
      return *this;

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
        std::cout << it2->second[i]->getPeak(ms_exp).getMZ() << "\t";
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
          ++peak_iter;
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
        penalty += 1000 * penalties.lWidth * pow(fabs(p_width_l - old_width_l), 2);
      if (p_width_r < 0)
      {
        penalty += 1e7 * penalties.rWidth * pow(fabs(p_width_r - old_width_r), 2);
      }
      else if (p_width_r < 1)
        penalty += 1000 * penalties.rWidth * pow(fabs(p_width_r - old_width_r), 2);
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
          ++peak_iter;
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
        penalty_l += 2000 * penalties.lWidth * (fabs(p_width_l - old_width_l));
      if (p_width_r < 0)
      {
        penalty_r += 1e7 * penalty_rwidth;
      }
      else if (p_width_r < 1)
        penalty_r += 2000 * penalties.rWidth * (fabs(p_width_r - old_width_r));
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
}
