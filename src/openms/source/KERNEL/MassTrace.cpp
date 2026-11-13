// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MassTrace.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MassTraceStatistics.h>

#include <algorithm>
#include <iterator>
#include <span>
#include <utility>

namespace OpenMS
{
    // must match MassTrace::MT_QUANTMETHOD enum!
    const std::string MassTrace::names_of_quantmethod[] = {"area", "median", "max_height"};

    MassTrace::MT_QUANTMETHOD MassTrace::getQuantMethod(const std::string& val)
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

    MassTrace::MassTrace(std::vector<PeakType>&& trace_peaks) :
            trace_peaks_(std::move(trace_peaks))
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
      return MassTraceStatistics::peakArea(trace_peaks_);
    }

    double MassTrace::computeIntensitySum() const
    {
      return MassTraceStatistics::intensitySum(trace_peaks_);
    }

    Size MassTrace::findMaxByIntPeak(bool use_smoothed_ints) const
    {
      return use_smoothed_ints ? MassTraceStatistics::apexIndex(trace_peaks_, smoothed_intensities_)
                               : MassTraceStatistics::apexIndex(trace_peaks_);
    }

    double MassTrace::estimateFWHM(bool use_smoothed_ints)
    {
      const MassTraceStatistics::FWHM estimate = use_smoothed_ints ? MassTraceStatistics::estimateFWHM(trace_peaks_, smoothed_intensities_)
                                                                   : MassTraceStatistics::estimateFWHM(trace_peaks_);
      fwhm_start_idx_ = estimate.start_idx;
      fwhm_end_idx_ = estimate.end_idx;

      // no points left or right of the apex: no FWHM. As before, a FWHM stored by an earlier call is kept.
      if (estimate.start_idx == 0 && estimate.end_idx == 0)
      {
        return 0;
      }

      fwhm_ = estimate.fwhm;
      return fwhm_;
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
      OPENMS_PRECONDITION(fwhm_end_idx_ < smoothed_intensities_.size(), "MassTrace was not smoothed before")

      // note that fwhm_end_idx_ is inclusive
      const Size count = fwhm_end_idx_ - fwhm_start_idx_ + 1;
      return MassTraceStatistics::peakArea(std::span<const PeakType>(trace_peaks_).subspan(fwhm_start_idx_, count),
                                           std::span<const double>(smoothed_intensities_).subspan(fwhm_start_idx_, count));
    }

    double MassTrace::computeFwhmArea() const
    {
      if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0)
      {
        return 0;
      }

      // note that fwhm_end_idx_ is inclusive
      const Size count = fwhm_end_idx_ - fwhm_start_idx_ + 1;
      return MassTraceStatistics::peakArea(std::span<const PeakType>(trace_peaks_).subspan(fwhm_start_idx_, count));
    }

    double MassTrace::getIntensity(bool smoothed) const
    {
      switch (quant_method_)
      {
        case MT_QUANT_AREA:
          return smoothed ? computeFwhmAreaSmooth() : computeFwhmArea();
        case MT_QUANT_MEDIAN:
          // median quantification (e.g. for direct infusion data) should work independently from smoothing
          return MassTraceStatistics::medianIntensity(trace_peaks_);
        case MT_QUANT_HEIGHT:
          return getMaxIntensity(smoothed);
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Member 'quant_method_' has unsupported value.",StringUtils::toStr(quant_method_));
      }
    }

    double MassTrace::getMaxIntensity(bool smoothed) const
    {
      return smoothed ? MassTraceStatistics::maxIntensity(smoothed_intensities_)
                      : MassTraceStatistics::maxIntensity(trace_peaks_);
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
      centroid_rt_ = MassTraceStatistics::weightedMeanRT(trace_peaks_);
    }

    void MassTrace::updateSmoothedWeightedMeanRT()
    {
      centroid_rt_ = MassTraceStatistics::smoothedWeightedMeanRT(trace_peaks_, smoothed_intensities_);
    }

    void MassTrace::updateSmoothedMaxRT()
    {
      centroid_rt_ = MassTraceStatistics::smoothedMaxRT(trace_peaks_, smoothed_intensities_);
    }

    void MassTrace::updateMedianRT()
    {
      centroid_rt_ = MassTraceStatistics::medianRT(trace_peaks_);
    }

    void MassTrace::updateMedianMZ()
    {
      centroid_mz_ = MassTraceStatistics::medianMZ(trace_peaks_);
    }

    void MassTrace::updateMeanMZ()
    {
      centroid_mz_ = MassTraceStatistics::meanMZ(trace_peaks_);
    }

    void MassTrace::updateWeightedMeanMZ()
    {
      centroid_mz_ = MassTraceStatistics::weightedMeanMZ(trace_peaks_);
    }

    void MassTrace::updateWeightedMZsd()
    {
      centroid_sd_ = MassTraceStatistics::weightedMZsd(trace_peaks_, centroid_mz_);
    }

} // end of MassTrace.cpp
