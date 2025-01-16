// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MassTrace.h>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
    // must match MassTrace::MT_QUANTMETHOD enum!
    const std::string MassTrace::names_of_quantmethod[] = {"area", "median", "max_height"};

    MassTrace::MT_QUANTMETHOD MassTrace::getQuantMethod(const String& val)
    {
      const std::string* qb = MassTrace::names_of_quantmethod;
      const std::string* qe = qb + (int)MassTrace::SIZE_OF_MT_QUANTMETHOD;
      const std::string* qm = std::find(qb, qe, val);
      return (MassTrace::MT_QUANTMETHOD)std::distance(qb, qm);
    }

    MassTrace::MassTrace(const std::list<PeakType>& trace_peaks) :
            
            trace_peaks_(),
            
            label_(),
            smoothed_intensities_()
            
    {
      trace_peaks_.reserve(trace_peaks.size());
      std::copy(trace_peaks.begin(), trace_peaks.end(), back_inserter(trace_peaks_));
    }

    MassTrace::MassTrace(const std::vector<PeakType>& trace_peaks) :
            
            trace_peaks_(trace_peaks),
            
            label_(),
            smoothed_intensities_()
            
    {
    }

    PeakType& MassTrace::operator[](const Size& mt_idx)
    {
      return trace_peaks_[mt_idx];
    }

    const PeakType& MassTrace::operator[](const Size& mt_idx) const
    {
      return trace_peaks_[mt_idx];
    }

    double MassTrace::computeSmoothedPeakArea() const
    {
      // sum all smoothed intensities in MassTrace which are non-negative
      double peak_area(0.0);

      double int_before = smoothed_intensities_[0];
      double rt_before = trace_peaks_.begin()->getRT();
      for (Size i = 1; i < smoothed_intensities_.size(); ++i)
      {
        if (smoothed_intensities_[i] > 0.0)
        {
          double rt_diff = trace_peaks_[i].getRT() - rt_before;
          peak_area += (int_before + trace_peaks_[i].getIntensity())/2 * rt_diff;
        }
        int_before = trace_peaks_[i].getIntensity();
        rt_before = trace_peaks_[i].getRT();
      }
      return peak_area;
    }

    double MassTrace::computePeakArea() const
    {
      double peak_area(0.0);

      if (trace_peaks_.empty())
      {
        return peak_area;
      }
      double int_before = trace_peaks_.begin()->getIntensity();
      double rt_before = trace_peaks_.begin()->getRT();
      for (const Peak2D& l_it : trace_peaks_)
      {
        double rt_diff = l_it.getRT() - rt_before;
        peak_area += (int_before + l_it.getIntensity())/2 * rt_diff;
        int_before = l_it.getIntensity();
        rt_before = l_it.getRT();
      }

      return peak_area;
    }


    double MassTrace::computeIntensitySum() const
    {
      return std::accumulate(trace_peaks_.begin(), trace_peaks_.end(), 0.0, 
        [](double sum, const Peak2D& peak) { return sum + peak.getIntensity(); });
    }    

    Size MassTrace::findMaxByIntPeak(bool use_smoothed_ints) const
    {
      if (use_smoothed_ints && smoothed_intensities_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
      }

      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace appears to be empty! Aborting...", String(trace_peaks_.size()));
      }

      double max_int;
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
        double act_int = use_smoothed_ints ? smoothed_intensities_[i] : trace_peaks_[i].getIntensity();

        if (act_int > max_int)
        {
          max_int = act_int;
          max_idx = i;
        }

      }

      return max_idx;
    }

    double MassTrace::estimateFWHM(bool use_smoothed_ints)
    {
      Size max_idx(this->findMaxByIntPeak(use_smoothed_ints));

      std::vector<double> tmp_ints;

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

      double half_max_int(tmp_ints[max_idx] / 2.0);

      // mass trace is empty OR no points left of apex in mass trace OR no points right of apex in mass trace
      if (tmp_ints.empty() || max_idx == 0 || max_idx == tmp_ints.size() - 1)
      {
        fwhm_start_idx_ = 0;
        fwhm_end_idx_ = 0;
        return 0;
      }

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

      // calculate RT of left side
      Size fwhm_left_top_idx = left_border + 1;
      double left_rt_bottom = trace_peaks_[left_border].getRT();
      double left_rt_top = trace_peaks_[fwhm_left_top_idx].getRT();
      double fwhm_begin_rt;
      // no points below half max -> don't extrapolate, use leftmost RT
      if (tmp_ints[left_border] > half_max_int)
      {
        fwhm_begin_rt = trace_peaks_[left_border].getRT();
      }
      else
      {
        fwhm_begin_rt = linearInterpolationAtY_(left_rt_bottom, left_rt_top,
                                                (tmp_ints[left_border]),
                                                (tmp_ints[fwhm_left_top_idx]),
                                                (half_max_int));
      }

      // calculate RT of right side
      Size fwhm_right_top_idx = right_border - 1;
      double right_rt_bottom = trace_peaks_[right_border].getRT();
      double right_rt_top = trace_peaks_[fwhm_right_top_idx].getRT();

      double fwhm_end_rt;
      // no points below half max -> don't extrapolate, use rightmost RT
      if (tmp_ints[right_border] > half_max_int)
      {
        fwhm_end_rt = trace_peaks_[right_border].getRT();
      }
      else
      {
        fwhm_end_rt = linearInterpolationAtY_(right_rt_top, right_rt_bottom,
                                              (tmp_ints[right_border]),
                                              (tmp_ints[fwhm_right_top_idx]),
                                              (half_max_int));
      }

      fwhm_ = std::fabs(fwhm_end_rt - fwhm_begin_rt);

      OPENMS_POSTCONDITION(fwhm_ <= std::fabs(trace_peaks_[right_border].getRT() - trace_peaks_[left_border].getRT()), "Interpolated fwhm is not smaller than extrapolated fwhm!")

      return fwhm_;
    }

    double MassTrace::linearInterpolationAtY_(double xA, double xB, double yA, double yB, double y_eval) const
    {
      OPENMS_PRECONDITION(yA <= y_eval && y_eval <= yB, "y_eval is not between yA and yB")
      // no solution -> return an estimate
      if (std::fabs(xA - xB) == 0 || std::fabs(yA - yB) == 0)  { return xA; }

      double xC = (xA + ((y_eval - yA) * (xB - xA) / (yB - yA)));
      OPENMS_POSTCONDITION(xA <= xC && xC <= xB, "xC is not between xA and xB");

      return xC;
    }

    /// determine if area or median is used for quantification
    /// @throw Exception::InvalidValue if SIZE_OF_MT_QUANTMETHOD is given
    void MassTrace::setQuantMethod(MassTrace::MT_QUANTMETHOD method)
    {
      if (method >= SIZE_OF_MT_QUANTMETHOD)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value of 'quant_method' cannot be 'SIZE_OF_MT_QUANTMETHOD'.", "");
      }
      quant_method_ = method;
    }

    /// check if area or median is used for quantification
    MassTrace::MT_QUANTMETHOD MassTrace::getQuantMethod() const
    {
      return quant_method_;
    }

    double MassTrace::computeFwhmAreaSmooth() const
    {
      if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
      {
        return 0;
      }

      double t_area(0.0);
      double int_before = smoothed_intensities_[fwhm_start_idx_];
      double rt_before = trace_peaks_[fwhm_start_idx_].getRT();
      // note '<=' operator, since fwhm_end_idx_ is inclusive!
      for (Size i = fwhm_start_idx_ + 1; i <= fwhm_end_idx_; ++i)
      {
        t_area += (int_before + smoothed_intensities_[i])/2 * (trace_peaks_[i].getRT() - rt_before);
        int_before = smoothed_intensities_[i];
        rt_before = trace_peaks_[i].getRT();
      }

      return t_area;
    }

    double MassTrace::computeFwhmArea() const
    {
      if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
      {
        return 0;
      }

      double t_area(0);
      double int_before = trace_peaks_[fwhm_start_idx_].getIntensity();
      double rt_before = trace_peaks_[fwhm_start_idx_].getRT();
      // note '<=' operator, since fwhm_end_idx_ is inclusive!
      for (Size i = fwhm_start_idx_ + 1; i <= fwhm_end_idx_; ++i)
      {
        t_area += (int_before + trace_peaks_[i].getIntensity())/2 * (trace_peaks_[i].getRT() - rt_before);
        int_before = trace_peaks_[i].getIntensity();
        rt_before = trace_peaks_[i].getRT();
      }

      return t_area;
    }

    double MassTrace::computeMedianIntensity_() const
    {
      // determine median of intensities
      double t_area(0);
      std::vector<double> ints;
      ints.reserve(trace_peaks_.size());
      for (Size i=0; i<trace_peaks_.size(); ++i)
      {
        ints.push_back(trace_peaks_[i].getIntensity());
      }
      sort(ints.begin(), ints.end());
      t_area = (ints.size() % 2 == 0 ? 0.5*(ints[ints.size()/2-1]+ints[ints.size()/2])
                                     : ints[ints.size()/2]);
      return t_area;
    }

    double MassTrace::getIntensity(bool smoothed) const
    {
      if (smoothed)
      { // will be removed soon
        switch (quant_method_)
        {
          case MT_QUANT_AREA:
            return computeFwhmAreaSmooth();
          case MT_QUANT_MEDIAN:
            // median quantification (e.g. for direct infusion data) should work indepentently from smoothing
            return computeMedianIntensity_();
          case MT_QUANT_HEIGHT:
            return getMaxIntensity(true);
          default:
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Member 'quant_method_' has unsupported value.", String(quant_method_));
        }

      }
      else
      {
        switch (quant_method_)
        {
          case MT_QUANT_AREA:
            return computeFwhmArea();
          case MT_QUANT_MEDIAN:
            return computeMedianIntensity_();
          case MT_QUANT_HEIGHT:
            return getMaxIntensity(false);
          default:
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Member 'quant_method_' has unsupported value.", String(quant_method_));

        }
      }
    }

    double MassTrace::getMaxIntensity(bool smoothed) const
    {
      double max_int(0.0);

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
      for (const Peak2D& l_it : trace_peaks_)
      {
        hull_points[i][0] = (l_it).getRT();
        hull_points[i][1] = (l_it).getMZ();
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
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
      }






      /* seems not to work with the way we compute the area in the code below -> as a result the RT values are outside the feature boundaries
      trace_area = this->computePeakArea();

      if (trace_area < std::numeric_limits<double>::epsilon())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peak area equals zero... impossible to compute weights!", String(trace_peaks_.size()));
      }
      the reason is because computePeakArea uses trapezoidal rule to compute the area, which is not the same as the sum of the intensities
      we could probably change the part below to also use trapezoidal rule to compute the trace area
      */

      double wmean_rt(0.0);
      double trace_area = 0;

      double rt_before = trace_peaks_[0].getRT();
      for (MassTrace::const_iterator l_it = trace_peaks_.begin() + 1; l_it != trace_peaks_.end(); ++l_it)
      {
        double rt_diff = l_it->getRT() - rt_before;                
        wmean_rt += l_it->getIntensity() * l_it->getRT() * rt_diff;
        rt_before = l_it->getRT();
        trace_area += l_it->getIntensity() * rt_diff;
      }
      centroid_rt_ = wmean_rt / trace_area;
    }

    void MassTrace::updateSmoothedWeightedMeanRT()
    {
      if (smoothed_intensities_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
      }

      double trace_area(0.0), wmean_rt(0.0);

      for (Size i = 0; i < smoothed_intensities_.size(); ++i)
      {
        if (smoothed_intensities_[i] > 0.0)
        {
          wmean_rt += smoothed_intensities_[i] * trace_peaks_[i].getRT();
          trace_area += smoothed_intensities_[i];
        }
      }

      if (trace_area < std::numeric_limits<double>::epsilon())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peak area equals to zero... impossible to compute weights!", String(trace_peaks_.size()));
      }

      centroid_rt_ = wmean_rt / trace_area;
    }

    void MassTrace::updateSmoothedMaxRT()
    {
      if (smoothed_intensities_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
      }

      double tmp_max(-1.0);
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
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Negative max intensity encountered!", String(tmp_max));
      }

      centroid_rt_ = trace_peaks_[max_idx].getRT();
    }

    void MassTrace::updateMedianRT()
    {
      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
      }

      if (trace_peaks_.size() == 1)
      {
        centroid_rt_ = (*(trace_peaks_.begin())).getRT();

        return;
      }

      // copy mz values to temp vec
      std::vector<double> temp_rt;

      for (const Peak2D& l_it : trace_peaks_)
      {
        temp_rt.push_back((l_it).getRT());
      }

      std::sort(temp_rt.begin(), temp_rt.end());

      Size temp_mz_size = temp_rt.size();
      Size mid = static_cast<Size>(temp_mz_size / 2.0);
      if ((temp_mz_size % 2) == 0)
      {
        centroid_rt_ = (temp_rt[mid - 1] +  temp_rt[mid]) / 2;
      }
      else
      {
        centroid_rt_ = temp_rt[mid];
      }


      return;
    }

    void MassTrace::updateMedianMZ()
    {
      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
      }

      if (trace_peaks_.size() == 1)
      {
        centroid_mz_ = (*(trace_peaks_.begin())).getMZ();

        return;
      }

      // copy mz values to temp vec
      std::vector<double> temp_mz;

      for (const Peak2D& l_it : trace_peaks_)
      {
        temp_mz.push_back((l_it).getMZ());
      }

      std::sort(temp_mz.begin(), temp_mz.end());

      Size temp_mz_size = temp_mz.size();
      Size mid = static_cast<Size>(temp_mz_size / 2.0);
      if ((temp_mz_size % 2) == 0)
      {
        centroid_mz_ = (temp_mz[mid - 1] +  temp_mz[mid]) / 2;
      }
      else
      {
        centroid_mz_ = temp_mz[mid];
      }

      return;
    }

    void MassTrace::updateMeanMZ()
    {
      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
      }

      Size trace_size = trace_peaks_.size();

      double sum_mz(0.0);

      for (const Peak2D& l_it : trace_peaks_)
      {
        sum_mz += (l_it).getMZ();
      }

      centroid_mz_ = sum_mz / trace_size;

      return;
    }

    void MassTrace::updateWeightedMeanMZ()
    {
      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
      }

      double weighted_sum(0.0);
      double total_weight(0.0);

      for (const Peak2D& l_it : trace_peaks_)
      {
        double w_i = (l_it).getIntensity();
        total_weight += w_i;
        weighted_sum += w_i * (l_it).getMZ();
      }

      if (total_weight < std::numeric_limits<double>::epsilon())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "All weights were equal to zero! Empty trace? Aborting...", String(total_weight));
      }

      centroid_mz_ = weighted_sum / total_weight;
    }

    void MassTrace::updateWeightedMZsd()
    {
      if (trace_peaks_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MassTrace is empty... std of MZ undefined!", String(trace_peaks_.size()));
      }

      double weighted_sum(0.0);
      double total_weight(0.0);

      for (const Peak2D& l_it : trace_peaks_)
      {
        double w_i = l_it.getIntensity();
        total_weight += w_i;
        weighted_sum += w_i * std::exp(2 * std::log(std::abs(l_it.getMZ() - centroid_mz_)));
      }

      if (total_weight < std::numeric_limits<double>::epsilon())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "All weights were equal to zero! Empty trace? Aborting...", String(total_weight));
      }

      centroid_sd_ = std::sqrt(weighted_sum) / std::sqrt(total_weight);

    }

} // end of MassTrace.cpp
