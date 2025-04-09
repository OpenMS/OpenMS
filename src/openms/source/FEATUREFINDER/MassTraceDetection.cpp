// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>

#include <OpenMS/MATH/StatisticFunctions.h>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
    MassTraceDetection::MassTraceDetection() :
            DefaultParamHandler("MassTraceDetection"), ProgressLogger()
    {
      defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass deviation (in ppm).");
      defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are removed as noise.");
      defaults_.setValue("chrom_peak_snr", 3.0, "Minimum intensity above noise_threshold_int (signal-to-noise) a peak should have to be considered an apex.");

      defaults_.setValue("reestimate_mt_sd", "true", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
      defaults_.setValidStrings("reestimate_mt_sd", {"true","false"});

      defaults_.setValue("quant_method", String(MassTrace::names_of_quantmethod[0]), "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.");
      defaults_.setValidStrings("quant_method", std::vector<std::string>(MassTrace::names_of_quantmethod, MassTrace::names_of_quantmethod +(int)MassTrace::SIZE_OF_MT_QUANTMETHOD));

      // advanced parameters
      defaults_.setValue("trace_termination_criterion", "outlier", "Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold.", {"advanced"});
      defaults_.setValidStrings("trace_termination_criterion", {"outlier","sample_rate"});
      defaults_.setValue("trace_termination_outliers", 5, "Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.", {"advanced"});

      defaults_.setValue("min_sample_rate", 0.5, "Minimum fraction of scans along the mass trace that must contain a peak.", {"advanced"});
      defaults_.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).", {"advanced"});
      defaults_.setValue("max_trace_length", -1.0, "Maximum expected length of a mass trace (in seconds). Set to a negative value to disable maximal length check during mass trace detection.", {"advanced"});

      defaultsToParam_();

      this->setLogType(CMD);
    }

    MassTraceDetection::~MassTraceDetection() = default;

    MassTraceDetection::Apex::Apex(double intensity, Size scan_idx, Size peak_idx):
      intensity(intensity),
      scan_idx(scan_idx),
      peak_idx(peak_idx)
    {}

    void MassTraceDetection::updateIterativeWeightedMeanMZ(const double& added_mz,
                                                           const double& added_int, double& centroid_mz, double& prev_counter,
                                                           double& prev_denom)
    {
      double new_weight(added_int);
      double new_mz(added_mz);

      double counter_tmp(1 + (new_weight * new_mz) / prev_counter);
      double denom_tmp(1 + (new_weight) / prev_denom);
      centroid_mz *= (counter_tmp / denom_tmp);
      prev_counter *= counter_tmp;
      prev_denom *= denom_tmp;

      return;
    }

    void MassTraceDetection::run(PeakMap::ConstAreaIterator& begin,
                                 PeakMap::ConstAreaIterator& end,
                                 std::vector<MassTrace>& found_masstraces)
    {
      PeakMap map;
      MSSpectrum current_spectrum;

      if (begin == end)
      {
        return;
      }

      for (; begin != end; ++begin)
      {
        // AreaIterator points on novel spectrum?
        if (begin.getRT() != current_spectrum.getRT())
        {
          // save new spectrum in map
          if (current_spectrum.getRT() != -1)
          {
            map.addSpectrum(current_spectrum);
          }
          current_spectrum.clear(false);
          current_spectrum.setRT(begin.getRT());
        }
        current_spectrum.push_back(*begin);
      }
      map.addSpectrum(current_spectrum);

      run(map, found_masstraces);
    }

// update function for FTL method

    void updateMeanEstimate(const double& x_t, double& mean_t, Size t)
    {
      mean_t +=  (1.0 / ((double)t + 1.0)) * (x_t - mean_t);
    }

    void updateSDEstimate(const double& x_t, const double& mean_t, double& sd_t, Size t)
    {
      double i(t);
      sd_t = (i / (i + 1)) * sd_t + (i / (i + 1) * (i + 1)) * (x_t - mean_t) * (x_t - mean_t);
      // std::cerr << "func:  " << tmp << " " << i << std::endl;
    }

    void updateWeightedSDEstimate(const PeakType& p, const double& mean_t1, double& sd_t, double& last_weights_sum)
    {
      double denom = last_weights_sum * sd_t * sd_t + p.getIntensity() * (p.getMZ() - mean_t1) * (p.getMZ() - mean_t1);
      double weights_sum = last_weights_sum + p.getIntensity();

      double tmp_sd = std::sqrt(denom / weights_sum);

      if (tmp_sd > std::numeric_limits<double>::epsilon())
      {
        sd_t = tmp_sd;
      }

      last_weights_sum = weights_sum;
    }

    void updateWeightedSDEstimateRobust(const PeakType& p, const double& mean_t1, double& sd_t, double& last_weights_sum)
    {
      double denom1 = std::log(last_weights_sum) + 2 * std::log(sd_t);
      double denom2 = std::log(p.getIntensity()) + 2 * std::log(std::abs(p.getMZ() - mean_t1));
      double denom = std::sqrt(std::exp(denom1) + std::exp(denom2));
      double weights_sum = last_weights_sum + p.getIntensity();
      double tmp_sd = denom / std::sqrt(weights_sum);

      if (tmp_sd > std::numeric_limits<double>::epsilon())
      {
        sd_t = tmp_sd;
      }

      last_weights_sum = weights_sum;
    }

    void computeWeightedSDEstimate(std::list<PeakType> tmp, const double& mean_t, double& sd_t, const double& /* lower_sd_bound */)
    {
      double denom(0.0), weights_sum(0.0);

      for (std::list<PeakType>::const_iterator l_it = tmp.begin(); l_it != tmp.end(); ++l_it)
      {
        denom += l_it->getIntensity() * (l_it->getMZ() - mean_t) * (l_it->getMZ() - mean_t);
        weights_sum += l_it->getIntensity();
      }

      double tmp_sd = std::sqrt(denom / (weights_sum));

      // std::cout << "tmp_sd" << tmp_sd << std::endl;

      if (tmp_sd > std::numeric_limits<double>::epsilon())
      {
        sd_t = tmp_sd;
      }

      return;
    }


    void MassTraceDetection::run(const PeakMap& input_exp, std::vector<MassTrace>& found_masstraces, const Size max_traces)
    {
      // make sure the output vector is empty
      found_masstraces.clear();

      // gather all peaks that are potential chromatographic peak apices
      //   - use work_exp for actual work (remove peaks below noise threshold)
      //   - store potential apices in chrom_apices
      PeakMap work_exp;
      std::vector<Apex> chrom_apices;

      Size total_peak_count(0);
      std::vector<Size> spec_offsets;
      spec_offsets.push_back(0);

      Size spectra_count(0);

      // *********************************************************** //
      //  Step 1: Detecting potential chromatographic apices
      // *********************************************************** //
      for (const MSSpectrum& it : input_exp)
      {
        // check if this is a MS1 survey scan
        if (it.getMSLevel() != 1)
        {
          continue;
        }
        std::vector<Size> indices_passing;
        for (Size peak_idx = 0; peak_idx < it.size(); ++peak_idx)
        {
          double tmp_peak_int((it)[peak_idx].getIntensity());
          if (tmp_peak_int > noise_threshold_int_)
          {
            // Assume that noise_threshold_int_ contains the noise level of the
            // data and we want to be chrom_peak_snr times above the noise level
            // --> add this peak as possible chromatographic apex
            if (tmp_peak_int > chrom_peak_snr_ * noise_threshold_int_)
            {
              chrom_apices.emplace_back(tmp_peak_int, spectra_count, indices_passing.size());
            }
            indices_passing.push_back(peak_idx);
            ++total_peak_count;
          }
        }
        PeakMap::SpectrumType tmp_spec(it);
        tmp_spec.select(indices_passing);
        work_exp.addSpectrum(tmp_spec);
        spec_offsets.push_back(spec_offsets.back() + tmp_spec.size());
        ++spectra_count;
      }

      if (spectra_count < 3)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Input map consists of too few MS1 spectra (less than 3!). Aborting...", String(spectra_count));
      }

      // discard last spectrum's offset
      spec_offsets.pop_back();

      std::sort(chrom_apices.begin(), chrom_apices.end(),
                [](const Apex & a,
                    const Apex & b) -> bool
      {
        return a.intensity < b.intensity;
      });

      // *********************************************************************
      // Step 2: start extending mass traces beginning with the apex peak (go
      // through all peaks in order of decreasing intensity)
      // *********************************************************************
      run_(chrom_apices, total_peak_count, work_exp, spec_offsets, found_masstraces, max_traces);

      return;
    } // end of MassTraceDetection::run

    void MassTraceDetection::run_(const std::vector<Apex>& chrom_apices,
                                  const Size total_peak_count,
                                  const PeakMap& work_exp,
                                  const std::vector<Size>& spec_offsets,
                                  std::vector<MassTrace>& found_masstraces,
                                  const Size max_traces)
    {
      boost::dynamic_bitset<> peak_visited(total_peak_count);
      Size trace_number(1);

      // check presence of FWHM meta data
      int fwhm_meta_idx(-1);
      Size fwhm_meta_count(0);
      for (Size i = 0; i < work_exp.size(); ++i)
      {
        if (!work_exp[i].getFloatDataArrays().empty() &&
            work_exp[i].getFloatDataArrays()[0].getName() == "FWHM_ppm")
        {
          if (work_exp[i].getFloatDataArrays()[0].size() != work_exp[i].size())
          { // float data should always have the same size as the corresponding array
            throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, work_exp[i].size());
          }
          fwhm_meta_idx = 0;
          ++fwhm_meta_count;
        }
      }
      if (fwhm_meta_count > 0 && fwhm_meta_count != work_exp.size())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("FWHM meta arrays are expected to be missing or present for all MS spectra [") + fwhm_meta_count + "/" + work_exp.size() + "].");
      }


      this->startProgress(0, total_peak_count, "mass trace detection");
      Size peaks_detected(0);

      for (auto m_it = chrom_apices.crbegin(); m_it != chrom_apices.crend(); ++m_it)
      {
        Size apex_scan_idx(m_it->scan_idx);
        Size apex_peak_idx(m_it->peak_idx);

        if (peak_visited[spec_offsets[apex_scan_idx] + apex_peak_idx])
        {
          continue;
        }

        Peak2D apex_peak;
        apex_peak.setRT(work_exp[apex_scan_idx].getRT());
        apex_peak.setMZ(work_exp[apex_scan_idx][apex_peak_idx].getMZ());
        apex_peak.setIntensity(work_exp[apex_scan_idx][apex_peak_idx].getIntensity());

        Size trace_up_idx(apex_scan_idx);
        Size trace_down_idx(apex_scan_idx);

        std::list<PeakType> current_trace;
        current_trace.push_back(apex_peak);
        std::vector<double> fwhms_mz; // peak-FWHM meta values of collected peaks

        // Initialization for the iterative version of weighted m/z mean calculation
        double centroid_mz(apex_peak.getMZ());
        double prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
        double prev_denom(apex_peak.getIntensity());

        updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

        std::vector<std::pair<Size, Size> > gathered_idx;
        gathered_idx.emplace_back(apex_scan_idx, apex_peak_idx);
        if (fwhm_meta_idx != -1)
        {
          fwhms_mz.push_back(work_exp[apex_scan_idx].getFloatDataArrays()[fwhm_meta_idx][apex_peak_idx]);
        }

        Size up_hitting_peak(0), down_hitting_peak(0);
        Size up_scan_counter(0), down_scan_counter(0);

        bool toggle_up = true, toggle_down = true;

        Size conseq_missed_peak_up(0), conseq_missed_peak_down(0);
        Size max_consecutive_missing(trace_termination_outliers_);

        double current_sample_rate(1.0);
        // Size min_scans_to_consider(std::floor((min_sample_rate_ /2)*10));
        Size min_scans_to_consider(5);

        // double outlier_ratio(0.3);

        // double ftl_mean(centroid_mz);
        double ftl_sd((centroid_mz / 1e6) * mass_error_ppm_);
        double intensity_so_far(apex_peak.getIntensity());

        while (((trace_down_idx > 0) && toggle_down) ||
               ((trace_up_idx < work_exp.size() - 1) && toggle_up)
                )
        {
          // *********************************************************** //
          // Step 2.1 MOVE DOWN in RT dim
          // *********************************************************** //
          if ((trace_down_idx > 0) && toggle_down)
          {
            const MSSpectrum& spec_trace_down = work_exp[trace_down_idx - 1];
            if (!spec_trace_down.empty())
            {
              Size next_down_peak_idx = spec_trace_down.findNearest(centroid_mz);
              double next_down_peak_mz = spec_trace_down[next_down_peak_idx].getMZ();
              double next_down_peak_int = spec_trace_down[next_down_peak_idx].getIntensity();

              double right_bound = centroid_mz + 3 * ftl_sd;
              double left_bound = centroid_mz - 3 * ftl_sd;

              if ((next_down_peak_mz <= right_bound) &&
                  (next_down_peak_mz >= left_bound) &&
                  !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx]
                      )
              {
                Peak2D next_peak;
                next_peak.setRT(spec_trace_down.getRT());
                next_peak.setMZ(next_down_peak_mz);
                next_peak.setIntensity(next_down_peak_int);

                current_trace.push_front(next_peak);
                // FWHM average
                if (fwhm_meta_idx != -1)
                {
                  fwhms_mz.push_back(spec_trace_down.getFloatDataArrays()[fwhm_meta_idx][next_down_peak_idx]);
                }
                // Update the m/z mean of the current trace as we added a new peak
                updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
                gathered_idx.emplace_back(trace_down_idx - 1, next_down_peak_idx);

                // Update the m/z variance dynamically
                if (reestimate_mt_sd_)           //  && (down_hitting_peak+1 > min_flank_scans))
                {
                  // if (ftl_t > min_fwhm_scans)
                  {
                    updateWeightedSDEstimateRobust(next_peak, centroid_mz, ftl_sd, intensity_so_far);
                  }
                }

                ++down_hitting_peak;
                conseq_missed_peak_down = 0;
              }
              else
              {
                ++conseq_missed_peak_down;
              }

            }
            --trace_down_idx;
            ++down_scan_counter;

            // trace termination criterion: max allowed number of
            // consecutive outliers reached OR cancel extension if
            // sampling_rate falls below min_sample_rate_
            if (trace_termination_criterion_ == "outlier")
            {
              if (conseq_missed_peak_down > max_consecutive_missing)
              {
                toggle_down = false;
              }
            }
            else if (trace_termination_criterion_ == "sample_rate")
            {
              current_sample_rate = (double)(down_hitting_peak + up_hitting_peak + 1) /
                                    (double)(down_scan_counter + up_scan_counter + 1);
              if (down_scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
              {
                // std::cout << "stopping down..." << std::endl;
                toggle_down = false;
              }
            }
          }

          // *********************************************************** //
          // Step 2.2 MOVE UP in RT dim
          // *********************************************************** //
          if ((trace_up_idx < work_exp.size() - 1) && toggle_up)
          {
            const MSSpectrum& spec_trace_up = work_exp[trace_up_idx + 1];
            if (!spec_trace_up.empty())
            {
              Size next_up_peak_idx = spec_trace_up.findNearest(centroid_mz);
              double next_up_peak_mz = spec_trace_up[next_up_peak_idx].getMZ();
              double next_up_peak_int = spec_trace_up[next_up_peak_idx].getIntensity();

              double right_bound = centroid_mz + 3 * ftl_sd;
              double left_bound = centroid_mz - 3 * ftl_sd;

              if ((next_up_peak_mz <= right_bound) &&
                  (next_up_peak_mz >= left_bound) &&
                  !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
              {
                Peak2D next_peak;
                next_peak.setRT(spec_trace_up.getRT());
                next_peak.setMZ(next_up_peak_mz);
                next_peak.setIntensity(next_up_peak_int);

                current_trace.push_back(next_peak);
                if (fwhm_meta_idx != -1)
                {
                  fwhms_mz.push_back(spec_trace_up.getFloatDataArrays()[fwhm_meta_idx][next_up_peak_idx]);
                }
                // Update the m/z mean of the current trace as we added a new peak
                updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
                gathered_idx.emplace_back(trace_up_idx + 1, next_up_peak_idx);

                // Update the m/z variance dynamically
                if (reestimate_mt_sd_)           //  && (up_hitting_peak+1 > min_flank_scans))
                {
                  // if (ftl_t > min_fwhm_scans)
                  {
                    updateWeightedSDEstimateRobust(next_peak, centroid_mz, ftl_sd, intensity_so_far);
                  }
                }

                ++up_hitting_peak;
                conseq_missed_peak_up = 0;

              }
              else
              {
                ++conseq_missed_peak_up;
              }

            }

            ++trace_up_idx;
            ++up_scan_counter;

            if (trace_termination_criterion_ == "outlier")
            {
              if (conseq_missed_peak_up > max_consecutive_missing)
              {
                toggle_up = false;
              }
            }
            else if (trace_termination_criterion_ == "sample_rate")
            {
              current_sample_rate = (double)(down_hitting_peak + up_hitting_peak + 1) / (double)(down_scan_counter + up_scan_counter + 1);

              if (up_scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
              {
                // std::cout << "stopping up" << std::endl;
                toggle_up = false;
              }
            }


          }

        }

        // std::cout << "current sr: " << current_sample_rate << std::endl;
        double num_scans(down_scan_counter + up_scan_counter + 1 - conseq_missed_peak_down - conseq_missed_peak_up);

        double mt_quality((double)current_trace.size() / (double)num_scans);
        // std::cout << "mt quality: " << mt_quality << std::endl;
        double rt_range(std::fabs(current_trace.rbegin()->getRT() - current_trace.begin()->getRT()));

        // *********************************************************** //
        // Step 2.3 check if minimum length and quality of mass trace criteria are met
        // *********************************************************** //
        bool max_trace_criteria = (max_trace_length_ < 0.0 || rt_range < max_trace_length_);
        if (rt_range >= min_trace_length_ && max_trace_criteria && mt_quality >= min_sample_rate_)
        {
          // std::cout << "T" << trace_number << "\t" << mt_quality << std::endl;

          // mark all peaks as visited
          for (Size i = 0; i < gathered_idx.size(); ++i)
          {
            peak_visited[spec_offsets[gathered_idx[i].first] +  gathered_idx[i].second] = true;
          }

          // create new MassTrace object and store collected peaks from list current_trace
          MassTrace new_trace(current_trace);
          new_trace.updateWeightedMeanRT();
          new_trace.updateWeightedMeanMZ();
          if (!fwhms_mz.empty())
          {
            new_trace.fwhm_mz_avg = Math::median(fwhms_mz.begin(), fwhms_mz.end());
          }
          new_trace.setQuantMethod(quant_method_);
          //new_trace.setCentroidSD(ftl_sd);
          new_trace.updateWeightedMZsd();
          new_trace.setLabel("T" + String(trace_number));
          ++trace_number;

          found_masstraces.push_back(new_trace);

          peaks_detected += new_trace.getSize();
          this->setProgress(peaks_detected);

          // check if we already reached the (optional) maximum number of traces
          if (max_traces > 0 && found_masstraces.size() == max_traces)
          {
            break;
          }
        }
      }

      this->endProgress();

    }

    void MassTraceDetection::updateMembers_()
    {
      mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
      noise_threshold_int_ = (double)param_.getValue("noise_threshold_int");
      chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
      quant_method_ = MassTrace::getQuantMethod((String)param_.getValue("quant_method").toString());

      trace_termination_criterion_ = (String)param_.getValue("trace_termination_criterion").toString();
      trace_termination_outliers_ = (Size)param_.getValue("trace_termination_outliers");
      min_sample_rate_ = (double)param_.getValue("min_sample_rate");
      min_trace_length_ = (double)param_.getValue("min_trace_length");
      max_trace_length_ = (double)param_.getValue("max_trace_length");
      reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
    }

}
