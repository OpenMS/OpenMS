// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MassTrace.h>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
  MassTrace::MassTrace() :
    trace_peaks_(),
    centroid_mz_(),
    centroid_sd_(),
    centroid_rt_(),
    label_(),
    smoothed_intensities_(),
    fwhm_(0.0),
    //fwhm_num_scans_(),
    scan_time_(1.0),
    fwhm_start_idx_(0),
    fwhm_end_idx_(0)
  {
  }

  MassTrace::MassTrace(const std::list<PeakType>& tmp_lst, const DoubleReal& scan_time) :
    centroid_mz_(),
    centroid_sd_(),
    centroid_rt_(),
    label_(),
    smoothed_intensities_(),
    fwhm_(0.0),
    fwhm_start_idx_(0),
    fwhm_end_idx_(0)
  {
    trace_peaks_.clear();

    for (std::list<PeakType>::const_iterator l_it = tmp_lst.begin(); l_it != tmp_lst.end(); ++l_it)
    {
      trace_peaks_.push_back((*l_it));
    }

    scan_time_ = scan_time;
  }

  MassTrace::MassTrace(const std::vector<PeakType>& tmp_vec, const DoubleReal& scan_time) :
    centroid_mz_(),
    centroid_sd_(),
    centroid_rt_(),
    label_(),
    smoothed_intensities_(),
    fwhm_(0.0),
    fwhm_start_idx_(0),
    fwhm_end_idx_(0)
  {
    trace_peaks_ = tmp_vec;
    scan_time_ = scan_time;
  }

  MassTrace::~MassTrace()
  {
  }

  MassTrace::MassTrace(const MassTrace& mt) :
    trace_peaks_(mt.trace_peaks_),
    centroid_mz_(mt.centroid_mz_),
    centroid_rt_(mt.centroid_rt_),
    centroid_sd_(mt.centroid_sd_),
    label_(mt.label_),
    smoothed_intensities_(mt.smoothed_intensities_),
    fwhm_(mt.fwhm_),
    //fwhm_num_scans_(mt.fwhm_num_scans_),
    scan_time_(mt.scan_time_),
    fwhm_start_idx_(mt.fwhm_start_idx_),
    fwhm_end_idx_(mt.fwhm_end_idx_)
  {
  }

  MassTrace& MassTrace::operator=(const MassTrace& rhs)
  {
    if (this == &rhs)
      return *this;

    trace_peaks_ = rhs.trace_peaks_;
    centroid_mz_ = rhs.centroid_mz_;
    centroid_rt_ = rhs.centroid_rt_;
    centroid_sd_ = rhs.centroid_sd_;
    label_ = rhs.label_;
    smoothed_intensities_ = rhs.smoothed_intensities_;
    fwhm_ = rhs.fwhm_;
    // fwhm_num_scans_ = rhs.fwhm_num_scans_;
    scan_time_ = rhs.scan_time_;
    fwhm_start_idx_ = rhs.fwhm_start_idx_;
    fwhm_end_idx_ = rhs.fwhm_end_idx_;

    return *this;
  }

  PeakType& MassTrace::operator[](const Size& mt_idx)
  {
    return trace_peaks_[mt_idx];
  }

  const PeakType& MassTrace::operator[](const Size& mt_idx) const
  {
    return trace_peaks_[mt_idx];
  }

  DoubleReal MassTrace::computeSmoothedPeakArea()
  {
    // sum all non-negative (smoothed!) intensities in MassTrace
    DoubleReal peak_area(0.0);

    for (Size i = 0; i < smoothed_intensities_.size(); ++i)
    {
      if (smoothed_intensities_[i] > 0.0)
      {
        peak_area += smoothed_intensities_[i];
      }
    }

    peak_area *= scan_time_;

    return peak_area;
  }

  DoubleReal MassTrace::computePeakArea()
  {
    DoubleReal peak_area(0.0);

    if (trace_peaks_.empty())
      return peak_area;

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      peak_area += (*l_it).getIntensity();
    }

    peak_area *= scan_time_;

    return peak_area;
  }

  DoubleReal MassTrace::computePeakArea() const
  {
    DoubleReal peak_area(0.0);

    if (trace_peaks_.empty())
      return peak_area;

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      peak_area += (*l_it).getIntensity();
    }

    peak_area *= scan_time_;

    return peak_area;
  }

  Size MassTrace::findMaxByIntPeak(bool use_smoothed_ints = false) const
  {
    if (use_smoothed_ints && smoothed_intensities_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
    }

    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace appears to be empty! Aborting...", String(trace_peaks_.size()));
    }

    DoubleReal max_int;
    Size max_idx(0);

    if (use_smoothed_ints)
    {
      max_int = smoothed_intensities_[0];
    }
    else
    {
      max_int = trace_peaks_.begin()->getIntensity();
    }

    for (Size i = 0; i < trace_peaks_.size(); ++i)
    {
      DoubleReal act_int = use_smoothed_ints ? smoothed_intensities_[i] : trace_peaks_[i].getIntensity();

      if (act_int > max_int)
      {
        max_int = act_int;
        max_idx = i;
      }

    }

    return max_idx;
  }

  DoubleReal MassTrace::estimateFWHM(bool use_smoothed_ints = false)
  {
    Size max_idx(this->findMaxByIntPeak(use_smoothed_ints));

    std::vector<DoubleReal> tmp_ints;

    if (use_smoothed_ints)
    {
      tmp_ints = smoothed_intensities_;
    }
    else
    {
      for (Size vec_idx = 0; vec_idx < trace_peaks_.size(); ++vec_idx)
      {
        tmp_ints.push_back(trace_peaks_[vec_idx].getIntensity());
      }
    }

    DoubleReal half_max_int(tmp_ints[max_idx] / 2.0);

    Size left_border(max_idx), right_border(max_idx);

    while (left_border > 0 && tmp_ints[left_border] >= half_max_int)
    {
      --left_border;
    }

    while (right_border + 1 < tmp_ints.size() && tmp_ints[right_border] >= half_max_int)
    {
      ++right_border;
    }


    fwhm_start_idx_ = left_border;
    fwhm_end_idx_ = right_border;
    fwhm_ = std::fabs(trace_peaks_[right_border].getRT() - trace_peaks_[left_border].getRT());

    return fwhm_;
  }

  DoubleReal MassTrace::computeFwhmAreaSmooth()
  {
    if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FWHM beginning/ending indices not computed? Aborting...", String(fwhm_start_idx_) + String(" ") + String(fwhm_end_idx_));
    }

    DoubleReal t_area(0.0);

    for (Size i = fwhm_start_idx_; i < fwhm_end_idx_; ++i)
    {
      t_area += smoothed_intensities_[i];
    }

    return t_area * scan_time_;
  }

  DoubleReal MassTrace::computeFwhmArea()
  {
    DoubleReal t_area(0.0);

    for (Size i = fwhm_start_idx_; i < fwhm_end_idx_; ++i)
    {
      t_area += trace_peaks_[i].getIntensity();
    }

    return t_area * scan_time_;
  }

  DoubleReal MassTrace::computeFwhmAreaSmoothRobust()
  {
    if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FWHM beginning/ending indices not computed? Aborting...", String(fwhm_start_idx_) + String(" ") + String(fwhm_end_idx_));
    }

    DoubleReal t_area(0.0);

    for (Size i = fwhm_start_idx_; i < fwhm_end_idx_; ++i)
    {
      DoubleReal rt_diff(std::fabs(trace_peaks_[i + 1].getRT() - trace_peaks_[i].getRT()));

      if (rt_diff < 1.5 * scan_time_)
      {
        t_area += smoothed_intensities_[i] * rt_diff;
      }
      else
      {
        DoubleReal averaged_int((smoothed_intensities_[i] + smoothed_intensities_[i + 1]) / 2.0);
        DoubleReal averaged_rt_diff(rt_diff / 2.0);

        t_area += smoothed_intensities_[i] * averaged_rt_diff;
        t_area += averaged_int * averaged_rt_diff;
      }

    }

    return t_area;
  }

  DoubleReal MassTrace::computeFwhmAreaRobust()
  {
    if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FWHM beginning/ending indices not computed? Aborting...", String(fwhm_start_idx_) + String(" ") + String(fwhm_end_idx_));
    }

    DoubleReal t_area(0.0);

    for (Size i = fwhm_start_idx_; i < fwhm_end_idx_; ++i)
    {
      DoubleReal rt_diff(std::fabs(trace_peaks_[i + 1].getRT() - trace_peaks_[i].getRT()));

      if (rt_diff < 1.5 * scan_time_)
      {
        t_area += trace_peaks_[i].getIntensity() * rt_diff;
      }
      else
      {
        DoubleReal averaged_int((trace_peaks_[i].getIntensity() + trace_peaks_[i + 1].getIntensity()) / 2.0);
        DoubleReal averaged_rt_diff(rt_diff / 2.0);

        t_area += trace_peaks_[i].getIntensity() * averaged_rt_diff;
        t_area += averaged_int * averaged_rt_diff;
      }

    }

    return t_area;
  }

  DoubleReal MassTrace::getIntensity(bool smoothed)
  {
    // std::cout << "area old:\t" << computeFwhmAreaSmooth() << "\tarea new:\t" << computeFwhmAreaSmoothRobust() << "area old:\t" << computeFwhmArea() << "\tarea new:\t" << computeFwhmAreaRobust() << std::endl;

    if (smoothed)
    {
      return computeFwhmAreaSmoothRobust();
    }

    return computeFwhmAreaRobust();
  }

  DoubleReal MassTrace::getMaxIntensity(bool smoothed)
  {
    DoubleReal max_int(0.0);

    if (smoothed)
    {
      for (Size i = 0; i < smoothed_intensities_.size(); ++i)
      {
        if (smoothed_intensities_[i] > max_int)
        {
          max_int = smoothed_intensities_[i];
        }
      }
    }
    else
    {
      for (Size i = 0; i < trace_peaks_.size(); ++i)
      {
        if (trace_peaks_[i].getIntensity() > max_int)
        {
          max_int = trace_peaks_[i].getIntensity();
        }
      }
    }

    return max_int;
  }

  DoubleReal MassTrace::getMaxIntensity(bool smoothed) const
  {
    DoubleReal max_int(0.0);

    if (smoothed)
    {
      for (Size i = 0; i < smoothed_intensities_.size(); ++i)
      {
        if (smoothed_intensities_[i] > max_int)
        {
          max_int = smoothed_intensities_[i];
        }
      }
    }
    else
    {
      for (Size i = 0; i < trace_peaks_.size(); ++i)
      {
        if (trace_peaks_[i].getIntensity() > max_int)
        {
          max_int = trace_peaks_[i].getIntensity();
        }
      }
    }

    return max_int;
  }

  ConvexHull2D MassTrace::getConvexhull() const
  {
    ConvexHull2D::PointArrayType hull_points(trace_peaks_.size());

    Size i = 0;
    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      hull_points[i][0] = (*l_it).getRT();
      hull_points[i][1] = (*l_it).getMZ();
      ++i;
    }

    ConvexHull2D hull;
    hull.addPoints(hull_points);

    return hull;
  }

  void MassTrace::updateWeightedMeanRT()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
    }

    DoubleReal trace_area(this->computePeakArea());

    if (trace_area < std::numeric_limits<DoubleReal>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peak area equals to zero... impossible to compute weights!", String(trace_peaks_.size()));
    }

    DoubleReal wmean_rt(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      wmean_rt += ((*l_it).getIntensity() * (*l_it).getRT()) * scan_time_;
    }

    centroid_rt_ = wmean_rt / trace_area;
  }

  void MassTrace::updateSmoothedWeightedMeanRT()
  {
    if (smoothed_intensities_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
    }

    DoubleReal trace_area(0.0), wmean_rt(0.0);

    for (Size i = 0; i < smoothed_intensities_.size(); ++i)
    {
      if (smoothed_intensities_[i] > 0.0)
      {
        wmean_rt += smoothed_intensities_[i] * trace_peaks_[i].getRT();
        trace_area += smoothed_intensities_[i];
      }
    }

    if (trace_area < std::numeric_limits<DoubleReal>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peak area equals to zero... impossible to compute weights!", String(trace_peaks_.size()));
    }

    centroid_rt_ = wmean_rt / trace_area;
  }

  void MassTrace::updateSmoothedMaxRT()
  {
    if (smoothed_intensities_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
    }

    DoubleReal tmp_max(-1.0);
    Size max_idx(0);

    for (Size i = 0; i < smoothed_intensities_.size(); ++i)
    {
      if (smoothed_intensities_[i] > tmp_max)
      {
        tmp_max = smoothed_intensities_[i];
        max_idx = i;
      }
    }

    if (tmp_max <= 0.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Negative max intensity encountered!", String(tmp_max));
    }

    centroid_rt_ = trace_peaks_[max_idx].getRT();
  }

  void MassTrace::updateMedianRT()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
    }

    if (trace_peaks_.size() == 1)
    {
      centroid_rt_ = (*(trace_peaks_.begin())).getRT();

      return;
    }

    // copy mz values to temp vec
    std::vector<DoubleReal> temp_rt;

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      temp_rt.push_back((*l_it).getRT());
    }

    std::sort(temp_rt.begin(), temp_rt.end());

    Size temp_mz_size = temp_rt.size();

    if ((temp_mz_size % 2) == 0)
    {
      centroid_rt_ = (temp_rt[std::floor(temp_mz_size / 2.0) - 1] +  temp_rt[std::floor(temp_mz_size / 2.0)]) / 2;
    }
    else
    {
      centroid_rt_ = temp_rt[std::floor(temp_mz_size / 2.0)];
    }


    return;
  }

  void MassTrace::updateMedianMZ()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
    }

    if (trace_peaks_.size() == 1)
    {
      centroid_mz_ = (*(trace_peaks_.begin())).getMZ();

      return;
    }

    // copy mz values to temp vec
    std::vector<DoubleReal> temp_mz;

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      temp_mz.push_back((*l_it).getMZ());
    }

    std::sort(temp_mz.begin(), temp_mz.end());

    Size temp_mz_size = temp_mz.size();

    if ((temp_mz_size % 2) == 0)
    {
      centroid_mz_ = (temp_mz[std::floor(temp_mz_size / 2.0) - 1] +  temp_mz[std::floor(temp_mz_size / 2.0)]) / 2;
    }
    else
    {
      centroid_mz_ = temp_mz[std::floor(temp_mz_size / 2.0)];
    }

    return;
  }

  void MassTrace::updateMeanMZ()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
    }

    Size trace_size = trace_peaks_.size();

    DoubleReal sum_mz(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      sum_mz += (*l_it).getMZ();
    }

    centroid_mz_ = sum_mz / trace_size;

    return;
  }

  void MassTrace::updateWeightedMeanMZ()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
    }

    DoubleReal weighted_sum(0.0);
    DoubleReal total_weight(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      DoubleReal w_i = (*l_it).getIntensity();
      total_weight += w_i;
      weighted_sum += w_i * (*l_it).getMZ();
    }

    if (total_weight < std::numeric_limits<DoubleReal>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "All weights were equal to zero! Empty trace? Aborting...", String(total_weight));
    }

    centroid_mz_ = weighted_sum / total_weight;
  }

  void MassTrace::updateWeightedMZsd()
  {
    if (trace_peaks_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... std of MZ undefined!", String(trace_peaks_.size()));
    }

    DoubleReal weighted_sum(0.0);
    DoubleReal total_weight(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
      DoubleReal w_i = (*l_it).getIntensity();
      total_weight += w_i;
      weighted_sum += w_i * std::exp(2 * std::log(std::abs((*l_it).getMZ() - centroid_mz_)));
    }

    if (total_weight < std::numeric_limits<DoubleReal>::epsilon())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "All weights were equal to zero! Empty trace? Aborting...", String(total_weight));
    }

    centroid_sd_ = std::sqrt(weighted_sum) / std::sqrt(total_weight);

  }

} // end of MassTrace.C
