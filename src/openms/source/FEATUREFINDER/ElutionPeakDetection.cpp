// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <boost/dynamic_bitset.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG_EPD

namespace OpenMS
{
  ElutionPeakDetection::ElutionPeakDetection() :
    DefaultParamHandler("ElutionPeakDetection"), ProgressLogger()
  {
    defaults_.setValue("chrom_fwhm", 5.0, "Expected full-width-at-half-maximum of chromatographic peaks (in seconds).");
    defaults_.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");

    // NOTE: the algorithm will only act upon the "fixed" value, if you would
    // like to use the "auto" setting, you will have to call filterByPeakWidth
    // yourself
    defaults_.setValue("width_filtering", "fixed", "Enable filtering of unlikely peak widths. The fixed setting filters out mass traces outside the [min_fwhm, max_fwhm] interval (set parameters accordingly!). The auto setting filters with the 5 and 95% quantiles of the peak width distribution.");
    defaults_.setValidStrings("width_filtering", {"off","fixed","auto"});
    defaults_.setValue("min_fwhm", 1.0, "Minimum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if parameter width_filtering is off or auto.", {"advanced"});
    defaults_.setValue("max_fwhm", 60.0, "Maximum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if parameter width_filtering is off or auto.", {"advanced"});

    defaults_.setValue("masstrace_snr_filtering", "false", "Apply post-filtering by signal-to-noise ratio after smoothing.", {"advanced"});
    defaults_.setValidStrings("masstrace_snr_filtering", {"true","false"});

    defaultsToParam_();
    this->setLogType(CMD);
  }

  ElutionPeakDetection::~ElutionPeakDetection() = default;

  double ElutionPeakDetection::computeMassTraceNoise(const MassTrace& tr)
  {
    // compute RMSE
    double squared_sum(0.0);
    std::vector<double> smooth_ints(tr.getSmoothedIntensities());

    for (Size i = 0; i < smooth_ints.size(); ++i)
    {
      squared_sum += (tr[i].getIntensity() - smooth_ints[i]) * (tr[i].getIntensity() - smooth_ints[i]);
    }

    double rmse(0.0);

    if (!smooth_ints.empty())
    {
      rmse = std::sqrt(squared_sum / smooth_ints.size());
    }

    return rmse;
  }

  double ElutionPeakDetection::computeMassTraceSNR(const MassTrace& tr)
  {
    double snr(0.0);

    if (tr.getSize() > 0)
    {
      double noise_area = computeMassTraceNoise(tr) * tr.getTraceLength();
      double signal_area = tr.computePeakArea();

      snr = signal_area / noise_area;
    }

    // std::cout << "snr " << snr << " ";

    return snr;
  }

  double ElutionPeakDetection::computeApexSNR(const MassTrace& tr)
  {
    double noise_level(computeMassTraceNoise(tr));

    double snr = 0;
    if (noise_level > 0.0)
    {
      double smoothed_apex_int(tr.getMaxIntensity(true));
      snr = smoothed_apex_int / noise_level;
    }

    // std::cout << "snr " << snr << " ";

    return snr;
  }

  void ElutionPeakDetection::findLocalExtrema(const MassTrace& tr, const Size& num_neighboring_peaks,
                                              std::vector<Size>& chrom_maxes, std::vector<Size>& chrom_mins) const
  {
    std::vector<double> smoothed_ints_vec(tr.getSmoothedIntensities());

    Size mt_length(smoothed_ints_vec.size());

    if (mt_length != tr.getSize())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "MassTrace was not smoothed before! Aborting...", String(smoothed_ints_vec.size()));
    }

    // first make sure that everything is cleared
    chrom_maxes.clear();
    chrom_mins.clear();

    // Remember which indices we have already used
    boost::dynamic_bitset<> used_idx(mt_length);

    // Extract RTs from the chromatogram and store them into vectors for index access
    // Store indices along with smoothed_ints to keep track of the peak order
    std::multimap<double, Size> intensity_indices;
    for (Size idx = 0; idx < mt_length; ++idx)
    {
      intensity_indices.insert(std::make_pair(smoothed_ints_vec[idx], idx));
    }

    // Step 1: Identify maxima
    for (std::multimap<double, Size>::const_iterator c_it = intensity_indices.begin(); c_it != intensity_indices.end(); ++c_it)
    {
      double ref_int = c_it->first;
      Size ref_idx = c_it->second;

      if (!(used_idx[ref_idx]) && ref_int > 0.0) 
      { // only allow unused points as seeds (potential local maximum)
        bool real_max = true;

        // Get start_idx and end_idx based on expected peak width
        Size start_idx(0);
        if (ref_idx > num_neighboring_peaks)
        {
          start_idx = ref_idx - num_neighboring_peaks;
        }
        Size end_idx = ref_idx + num_neighboring_peaks;
        if (end_idx > mt_length)
        {
          end_idx = mt_length;
        }

        // Identify putative peak between start_idx and end_idx, now check if
        // no other maxima exist within the expected boundaries (check whether
        // ref_int is higher than all smoothed intensities within the
        // boundaries).
        for (Size j = start_idx; j < end_idx; ++j)
        {
          if (j == ref_idx)
          { // skip seed
            continue;
          }

          if (used_idx[j])
          { // peak has already been collected?
            if (smoothed_ints_vec[j] > ref_int)
            { // break if higher intensity
              real_max = false;
              break;
            }
            else
            { // skip if only a low intensity peak (e.g. flanks of elution profile)
              continue;
            }
          }

          if (smoothed_ints_vec[j] > ref_int)
          {
            real_max = false;
            break;
          }
        }

        // If no other maxima exists, then add the current one to the list and
        // mark all indices as used
        if (real_max)
        {
          chrom_maxes.push_back(ref_idx);

          for (Size j = start_idx; j < end_idx; ++j)
          {
            used_idx[j] = true;
          }
        }
      }
    }

    std::sort(chrom_maxes.begin(), chrom_maxes.end());

    // Step 2: Identify minima using bisection between two maxima
    if (chrom_maxes.size() > 1)
    {
      // Keep track of two maxima
      Size left_idx(0), right_idx(1);
      while (left_idx < right_idx && right_idx < chrom_maxes.size())
      {

        // 2.1 Perform bisection between the two maxima to find potential minimum
        Size left_bound(chrom_maxes[left_idx] + 1);
        Size right_bound(chrom_maxes[right_idx] - 1);
        while ((left_bound + 1) < right_bound)
        {
          // Identify middle between two bounds
          double mid_dist((right_bound - left_bound) / 2.0);
          Size mid_element_idx(left_bound + std::floor(mid_dist));
          double mid_element_int = smoothed_ints_vec[mid_element_idx];

          // Walk to the left if the slope is positive here
          if (mid_element_int <= smoothed_ints_vec[mid_element_idx + 1])
          {
            right_bound = mid_element_idx;
          }
          // else walk to the right ... 
          else
          {
            left_bound = mid_element_idx;
          }

        }

        // 2.2 Choose minimum (either left_bound or right_bound) and get minimal RT / Intensity
        Size min_rt((smoothed_ints_vec[left_bound] < smoothed_ints_vec[right_bound]) ? left_bound : right_bound);
        double min_int(1.0);
        if (smoothed_ints_vec[min_rt] > min_int)
        {
          min_int = smoothed_ints_vec[min_rt];
        }

        // 2.3 Compute distance and intensities
        double left_max_int(smoothed_ints_vec[chrom_maxes[left_idx]]);
        double right_max_int(smoothed_ints_vec[chrom_maxes[right_idx]]);

        double left_rt(tr[chrom_maxes[left_idx]].getRT());
        double mid_rt(tr[min_rt].getRT());
        double right_rt(tr[chrom_maxes[right_idx]].getRT());

        // compute the distance from the two maxima to the new minima 
        double left_dist(std::fabs(mid_rt - left_rt));
        double right_dist(std::fabs(right_rt - mid_rt));
        double min_dist(min_fwhm_ / 2.0);

        // out debug info
#ifdef DEBUG_EPD
        std::cout << "findLocalExtrema: Identified potential minimum " << std::endl;
        std::cout << "    " << tr.getLabel() << ": left_idx,right_idx " << left_idx << "," << right_idx << 
          ":" << left_max_int << " min: " << min_int << " " << right_max_int << 
          " l " << left_rt << " r " << right_rt << " m " << mid_rt << std::endl;
        std::cout << "    Int: min " << min_int << ", left: " << left_max_int << ", right: " << right_max_int << std::endl;
        std::cout << "    Distance: min " << min_dist << ", left: " << left_dist << ", right: " << right_dist << std::endl;
#endif

        // 2.4 Decide whether to split the masstrace (introduce a minimum):
        // i)  the maxima intensity should be at least 2x above the minimum for a split
        // ii) check that splitting the trace would not create peaks smaller than min_dist 
        if (left_max_int / min_int >= 2.0
           && right_max_int / min_int >= 2.0
           && left_dist >= min_dist
           && right_dist >= min_dist)
        {
#ifdef DEBUG_EPD
        std::cout << "    -> add new minima " << ": left_idx,right_idx " << left_idx << "," << right_idx << 
          " l " << left_rt << " r " << right_rt << " m " << mid_rt << std::endl;
#endif

          chrom_mins.push_back(min_rt);
          left_idx = right_idx;
          ++right_idx;
        }
        else
        {
          // keep one of the maxima (the one with higher intensity), replace
          // the other with the next in RT
          if (left_max_int > right_max_int)
          {
            ++right_idx;
          }
          else
          {
            left_idx = right_idx;
            ++right_idx;
          }
        }
      }
    }

    return;
  }

  void ElutionPeakDetection::detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces)
  {
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    detectElutionPeaks_(mt, single_mtraces);
    return;
  }

  void ElutionPeakDetection::detectPeaks(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& single_mtraces)
  {
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    this->startProgress(0, mt_vec.size(), "elution peak detection");
    Size progress(0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < (SignedSize) mt_vec.size(); ++i)
    {
      IF_MASTERTHREAD this->setProgress(progress);

#ifdef _OPENMP
#pragma omp atomic
#endif
      ++progress;

      // push_back to 'single_mtraces' is protected, so threading is ok
      detectElutionPeaks_(mt_vec[i], single_mtraces);
    }

    this->endProgress();

    return;
  }

  void ElutionPeakDetection::filterByPeakWidth(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& filt_mtraces)
  {
    filt_mtraces.clear();

    std::multimap<double, Size> sorted_by_peakwidth;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
      sorted_by_peakwidth.insert(std::make_pair(mt_vec[i].estimateFWHM(true), i));
    }

    double mapsize(sorted_by_peakwidth.size());
    Size lower_quartile_idx(std::floor(mapsize * 0.05));
    Size upper_quartile_idx(std::floor(mapsize * 0.95));
    Size count_mt(0);

    // filter out mass traces below lower quartile and above upper quartile
    for (std::multimap<double, Size>::const_iterator m_it = sorted_by_peakwidth.begin(); m_it != sorted_by_peakwidth.end(); ++m_it)
    {
      if (count_mt >= lower_quartile_idx && count_mt <= upper_quartile_idx)
      {
        // std::cout << "pw added " << m_it->first << std::endl;
        filt_mtraces.push_back(mt_vec[m_it->second]);
      }
      ++count_mt;
    }

    std::cout << "pw low: " << filt_mtraces[0].estimateFWHM(true) << " " << " pw high: " << filt_mtraces[filt_mtraces.size() - 1].estimateFWHM(true) << std::endl;

    return;
  }

  void ElutionPeakDetection::detectElutionPeaks_(MassTrace& mt, std::vector<MassTrace>& single_mtraces)
  {

    // *********************************************************************
    // Step 1: Smooth data
    // *********************************************************************

    double scan_time(mt.getAverageMS1CycleTime());
    Size win_size = std::ceil(chrom_fwhm_ / scan_time);

    // add smoothed data (original data is still accessible)
    smoothData(mt, static_cast<Int>(win_size));

#ifdef DEBUG_EPD
    Size i = 0;
    std::cout << "*****" << std::endl;
    std::cout << "   finding elution peaks in mass traces RT "  << mt.getCentroidRT()  << " / mz " << mt.getCentroidMZ() << std::endl;
    std::cout << "   used for smoothing: win_size "  << win_size << " FWHM scan num " /* << mt.getFWHMScansNum() */ << std::endl;
    std::cout << "*****" << std::endl;
    for (MassTrace::const_iterator mt_it = mt.begin(); mt_it != mt.end(); ++mt_it)
    {
      // std::cout << mt_it->getIntensity() << " " << mt.getSmoothedIntensities()[i] << std::endl;
      ++i;
    }
    std::cout << "*****" << std::endl;
#endif

    // *********************************************************************
    // Step 2: Identify local maxima and minima
    // *********************************************************************
    std::vector<Size> maxes, mins;
    findLocalExtrema(mt, win_size / 2, maxes, mins);

#ifdef DEBUG_EPD
    std::cout << "findLocalExtrema returned: maxima " << maxes.size() << " / minima " << mins.size() << std::endl;
#endif

    // *********************************************************************
    // Step 3: Split mass trace according to detected peaks
    // *********************************************************************

    // if only one maximum exists: finished!
    if (maxes.size() == 1)
    {
      bool pw_ok = true;
      bool snr_ok = true;

      // *********************************************************************
      // Step 3.1: check mass trace length criteria (if fixed filter is enabled)
      // *********************************************************************
      if (pw_filtering_ == "fixed")
      {
        double act_fwhm(mt.estimateFWHM(true));
        if (act_fwhm < min_fwhm_ || act_fwhm > max_fwhm_)
        {
          pw_ok = false;
        }
      }

      // *********************************************************************
      // Step 3.2: check mass trace signal to noise filter criteria
      // *********************************************************************
      if (mt_snr_filtering_)
      {
        if (computeApexSNR(mt) < chrom_peak_snr_)
        {
          snr_ok = false;
        }
      }

      if (pw_ok && snr_ok)
      {
        mt.updateSmoothedMaxRT();

        if (pw_filtering_ != "fixed")
        {
          mt.estimateFWHM(true);
        }

#ifdef _OPENMP
#pragma omp critical (OPENMS_ElutionPeakDetection_mtraces)
#endif
        {
          single_mtraces.push_back(mt);
        }

      }
    }
    else if (maxes.empty())
    {
      return;
    }
    else // split mt to sub-traces
    {
      MassTrace::const_iterator cp_it = mt.begin();
      Size last_idx(0);

      // add last data point as last minimum (to grep the last chunk of the MT)
      mins.push_back(mt.getSize() - 1);

      for (Size min_idx = 0; min_idx < mins.size(); ++min_idx)
      {

        // *********************************************************************
        // Step 3.1: Create new mass trace (sub-trace between cp_it and split point)
        // *********************************************************************
        std::vector<PeakType> tmp_mt;
        std::vector<double> smoothed_tmp;
        while (last_idx <= mins[min_idx])
        {
          tmp_mt.push_back(*cp_it);
          smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
          ++cp_it;
          ++last_idx;
        }

        // Create new mass trace, copy smoothed intensities
        MassTrace new_mt(tmp_mt);
        new_mt.setSmoothedIntensities(smoothed_tmp);

        // check filter criteria
        bool pw_ok = true;
        bool snr_ok = true;

        // *********************************************************************
        // Step 3.2: check mass trace length criteria (if fixed filter is enabled)
        // *********************************************************************
        if (pw_filtering_ == "fixed")
        {
          double act_fwhm(new_mt.estimateFWHM(true));
          if (act_fwhm < min_fwhm_ || act_fwhm > max_fwhm_)
          {
            pw_ok = false;
          }
        }

        // *********************************************************************
        // Step 3.3: check mass trace signal to noise filter criteria
        // *********************************************************************
        if (mt_snr_filtering_)
        {
          if (computeApexSNR(mt) < chrom_peak_snr_)
          {
            snr_ok = false;
          }
        }

        if (pw_ok && snr_ok)
        {
          // set label of sub-trace
          new_mt.setLabel(mt.getLabel() + "." + String(min_idx + 1));
          new_mt.updateSmoothedMaxRT();
          new_mt.updateWeightedMeanMZ();
          new_mt.updateWeightedMZsd();
          new_mt.setQuantMethod(mt.getQuantMethod());
          if (pw_filtering_ != "fixed")
          {
            new_mt.estimateFWHM(true);
          }

#ifdef _OPENMP
#pragma omp critical (OPENMS_ElutionPeakDetection_mtraces)
#endif
          {
            single_mtraces.push_back(new_mt);
          }
        }
      }

    }
    return;
  }

  void ElutionPeakDetection::smoothData(MassTrace& mt, int win_size) const
  {
    // alternative smoothing using SavitzkyGolay
    // looking at the unit test, this method gives better fits than lowess smoothing
    // reference paper uses lowess smoothing

    MSSpectrum spectrum;
    for (Size i = 0; i != mt.getSize(); ++i)
    {
      spectrum.push_back(Peak1D(mt[i].getRT(), mt[i].getIntensity()));
    }

    SavitzkyGolayFilter sg;
    Param param;
    param.setValue("polynomial_order", 2);
    param.setValue("frame_length", std::max(3, win_size)); // frame length must be at least polynomial_order+1, otherwise SG will fail
    sg.setParameters(param);
    sg.filter(spectrum);
    MSSpectrum::iterator iter = spectrum.begin();
    std::vector<double> smoothed_intensities;
    for (; iter != spectrum.end(); ++iter)
    {
      smoothed_intensities.push_back(iter->getIntensity());
    }
    mt.setSmoothedIntensities(smoothed_intensities);
    //alternative end

    // std::cout << "win_size elution: " << scan_time << " " << win_size << std::endl;

    // if there is no previous FWHM estimation... do it now
    //    if (win_size == 0)
    //    {
    //        mt.estimateFWHM(false); // estimate FWHM
    //        win_size = mt.getFWHMScansNum();
    //    }

    // use one global window size for all mass traces to smooth
    //  std::vector<double> rts, ints;
    //
    //  for (MassTrace::const_iterator c_it = mt.begin(); c_it != mt.end(); ++c_it)
    //  {
    //      rts.push_back(c_it->getRT());
    //      ints.push_back(c_it->getIntensity());
    //  }
    //  LowessSmoothing lowess_smooth;
    //  Param lowess_params;
    //  lowess_params.setValue("window_size", win_size);
    //  lowess_smooth.setParameters(lowess_params);
    //  std::vector<double> smoothed_data;
    //  lowess_smooth.smoothData(rts, ints, smoothed_data);
    //  mt.setSmoothedIntensities(smoothed_data);
  }

  void ElutionPeakDetection::updateMembers_()
  {
    chrom_fwhm_ = (double)param_.getValue("chrom_fwhm");
    chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
    min_fwhm_ = (double)param_.getValue("min_fwhm");
    max_fwhm_ = (double)param_.getValue("max_fwhm");

    pw_filtering_ = param_.getValue("width_filtering").toString();
    mt_snr_filtering_ = param_.getValue("masstrace_snr_filtering").toBool();
  }

} //namespace OpenMS
