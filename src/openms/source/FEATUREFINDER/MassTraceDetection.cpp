// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>

#include <OpenMS/MATH/StatisticFunctions.h>

#include <boost/dynamic_bitset.hpp>

#include <OpenMS/KERNEL/SpectrumHelper.h>

#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
    MassTraceDetection::MassTraceDetection() :
            DefaultParamHandler("MassTraceDetection"), ProgressLogger()
    {
      defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass deviation (in ppm).");
      defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are removed as noise.");
      defaults_.setValue("chrom_peak_snr", 3.0, "Minimum intensity above noise_threshold_int (signal-to-noise) a peak should have to be considered an apex.");
      defaults_.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/k0).");

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

    void MassTraceDetection::updateIterativeWeightedMean_(const double& added_value,
                                                         const double& added_intensity,
                                                         double& centroid_value,
                                                         double& prev_counter,
                                                         double& prev_denom)
    {
      double new_weight = added_intensity;
      double new_val = added_value;

      double counter_tmp = 1.0 + (new_weight * new_val) / prev_counter;
      double denom_tmp = 1.0 + (new_weight) / prev_denom;

      centroid_value *= (counter_tmp / denom_tmp);
      prev_counter *= counter_tmp;
      prev_denom *= denom_tmp;
    }

    // detect presence of ion mobility data and locate all relevant meta array indices.
    void MassTraceDetection::getIMIndices_(
      const PeakMap& spectra,
      int& fwhm_meta_idx, bool& has_fwhm_mz,
      int& im_idx, bool& has_centroid_im,
      int& im_fwhm_idx, bool& has_fwhm_im
    ) const
    {
      for (const auto& spec : spectra)
      {
        const auto& fda = spec.getFloatDataArrays();
        if (!fda.empty())
        {
          auto it_fwhm = getDataArrayByName(fda, Constants::UserParam::FWHM_MZ_ppm);
          auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
          auto it_imf = getDataArrayByName(fda, Constants::UserParam::FWHM_IM);

          if (it_fwhm != fda.end())
          {
            fwhm_meta_idx = std::distance(fda.begin(), it_fwhm);
            has_fwhm_mz = true;
          }
          if (it_im != fda.end())
          {
            im_idx = std::distance(fda.begin(), it_im);
            has_centroid_im = true;
          }
          if (it_imf != fda.end())
          {
            im_fwhm_idx = std::distance(fda.begin(), it_imf);
            has_fwhm_im = true;
          }
          break;
        }
      }

      // validate that all meta arrays are consistently present or absent across all spectra
      auto validate_meta_array = [&](const String& name, int idx) {
        if (idx == -1) return;

        Size valid_count = 0;
        for (const auto& spec : spectra)
        {
          const auto& fda = spec.getFloatDataArrays();
          if (!fda.empty() && idx < (int)fda.size() && fda[idx].getName() == name)
          {
            if (fda[idx].size() != spec.size())
            {
              throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, spec.size());
            }
            ++valid_count;
          }
        }

        if (valid_count > 0 && valid_count != spectra.size())
        {
          throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                name + " meta arrays must be consistently present or absent across all MS spectra ["
                                                  + String(valid_count) + "/" + String(spectra.size()) + "].");
        }
      };

      validate_meta_array(Constants::UserParam::FWHM_MZ_ppm, fwhm_meta_idx);
      validate_meta_array(Constants::UserParam::ION_MOBILITY, im_idx);
      validate_meta_array(Constants::UserParam::FWHM_IM, im_fwhm_idx);
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

      std::stable_sort(chrom_apices.begin(), chrom_apices.end(),
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

    MassTraceDetection::PeakCandidate MassTraceDetection::findBestPeak_(
        const MSSpectrum& spectrum,
        double centroid_mz,
        double ftl_sd,
        double centroid_im) const
    {
      PeakCandidate candidate;
      
      if (spectrum.empty())
      {
        return candidate;
      }

      double right_bound = centroid_mz + 3 * ftl_sd;
      double left_bound = centroid_mz - 3 * ftl_sd;

      if (!has_centroid_im_)
      {
        // Standard LC-MS data: find peak closest to target m/z
        candidate.idx = spectrum.findNearest(centroid_mz);
        candidate.mz = spectrum[candidate.idx].getMZ();
        candidate.intensity = spectrum[candidate.idx].getIntensity();
        candidate.found = true;
      }
      else
      {
        // LC-IMS-MS data: find best peak considering both m/z and ion mobility
        double right_bound_im = centroid_im + ion_mobility_tolerance_;
        double left_bound_im = centroid_im - ion_mobility_tolerance_;

        auto left_bound_it = spectrum.MZBegin(left_bound);
        Size idx_left = left_bound_it - spectrum.begin();
        auto right_bound_it = spectrum.MZEnd(right_bound);
        Size idx_right = right_bound_it - spectrum.begin();

        // Search within m/z window for peak with closest ion mobility match
        for (Size i = idx_left; i < idx_right; ++i)
        {
          double im_value = spectrum.getFloatDataArrays()[ion_mobility_idx_][i];
          
          if (im_value >= left_bound_im && im_value <= right_bound_im)
          {
            if (!candidate.found || std::abs(spectrum[i].getMZ() - centroid_mz) < std::abs(candidate.mz - centroid_mz))
            {
              candidate.idx = i;
              candidate.mz = spectrum[i].getMZ();
              candidate.intensity = spectrum[i].getIntensity();
              candidate.im = im_value;
              candidate.found = true;
            }
          }
        }
      }

      return candidate;
    }

    bool MassTraceDetection::isPeakAcceptable_(
        const PeakCandidate& candidate,
        double centroid_mz,
        double ftl_sd,
        double centroid_im,
        Size spectrum_idx,
        const std::vector<Size>& spec_offsets,
        const boost::dynamic_bitset<>& peak_visited) const
    {
      if (!candidate.found)
      {
        return false;
      }

      double right_bound = centroid_mz + 3 * ftl_sd;
      double left_bound = centroid_mz - 3 * ftl_sd;

      // Peak must fall within m/z tolerance window (±3 standard deviations)
      if (!((candidate.mz <= right_bound) && (candidate.mz >= left_bound)))
      {
        return false;
      }

      // For ion mobility data, peak must also fall within IM tolerance window
      if (has_centroid_im_)
      {
        double right_bound_im = centroid_im + ion_mobility_tolerance_;
        double left_bound_im = centroid_im - ion_mobility_tolerance_;
        
        if (candidate.im < left_bound_im || candidate.im > right_bound_im)
        {
          return false;
        }
      }

      // Peak must not have been used in another trace already
      return !peak_visited[spec_offsets[spectrum_idx] + candidate.idx];
    }

    void MassTraceDetection::processPeak_(
        const PeakCandidate& candidate,
        const MSSpectrum& spectrum,
        std::list<PeakType>& current_trace,
        std::vector<std::pair<Size, Size>>& gathered_idx,
        std::vector<double>& fwhms_mz,
        std::vector<double>& fwhms_im,
        double& centroid_mz,
        double& centroid_im,
        double& prev_counter,
        double& prev_denom,
        double& prev_counter_im,
        double& prev_denom_im,
        double& ftl_sd,
        double& intensity_so_far,
        Size spectrum_idx,
        bool is_upward_extension)
    {
      Peak2D next_peak;
      next_peak.setRT(spectrum.getRT());
      next_peak.setMZ(candidate.mz);
      next_peak.setIntensity(candidate.intensity);

      // Add peak to the growing mass trace
      if (is_upward_extension)
      {
        current_trace.push_back(next_peak);
      }
      else
      {
        current_trace.push_front(next_peak);
      }

      // Update trace centroid m/z with intensity-weighted average
      updateIterativeWeightedMean_(candidate.mz, candidate.intensity, centroid_mz, prev_counter, prev_denom);
      gathered_idx.emplace_back(spectrum_idx, candidate.idx);

      // Update ion mobility centroid if available
      if (has_centroid_im_)
      {
        updateIterativeWeightedMean_(candidate.im, candidate.intensity, centroid_im, prev_counter_im, prev_denom_im);
      }

      // Collect FWHM metadata for trace quality assessment
      if (has_fwhm_mz_)
      {
        fwhms_mz.push_back(spectrum.getFloatDataArrays()[fwhm_meta_idx_][candidate.idx]);
      }
      if (has_fwhm_im_)
      {
        fwhms_im.push_back(spectrum.getFloatDataArrays()[im_fwhm_idx_][candidate.idx]);
      }

      // Dynamically adjust m/z tolerance based on observed variance
      if (reestimate_mt_sd_)
      {
        updateWeightedSDEstimateRobust(next_peak, centroid_mz, ftl_sd, intensity_so_far);
      }
    }


    bool MassTraceDetection::isTraceValid_(
        const std::list<PeakType>& trace,
        Size total_scans_visited,
        Size consecutive_missed_down,
        Size consecutive_missed_up) const
    {
      double rt_range = std::fabs(trace.rbegin()->getRT() - trace.begin()->getRT());
      
      // Check length criteria
      if (rt_range < min_trace_length_)
      {
        return false;
      }
      
      if (max_trace_length_ >= 0.0 && rt_range > max_trace_length_)
      {
        return false;
      }

      // Check quality (sample rate)
      Size adjusted_scans = total_scans_visited - consecutive_missed_down - consecutive_missed_up;
      double mt_quality = static_cast<double>(trace.size()) / static_cast<double>(adjusted_scans);
      
      return mt_quality >= min_sample_rate_;
    }

    void MassTraceDetection::run_(const std::vector<Apex>& chrom_apices,
                                  const Size total_peak_count,
                                  const PeakMap& work_exp,
                                  const std::vector<Size>& spec_offsets,
                                  std::vector<MassTrace>& found_masstraces,
                                  const Size max_traces)
    {
      boost::dynamic_bitset<> peak_visited(total_peak_count);
      Size trace_number(1);

      // Detect ion mobility and FWHM metadata arrays in the dataset
      getIMIndices_(work_exp,
                    fwhm_meta_idx_, has_fwhm_mz_,
                    ion_mobility_idx_, has_centroid_im_,
                    im_fwhm_idx_, has_fwhm_im_);

      this->startProgress(0, total_peak_count, "mass trace detection");
      Size peaks_detected(0);

      for (auto m_it = chrom_apices.crbegin(); m_it != chrom_apices.crend(); ++m_it)
      {
        Size apex_scan_idx(m_it->scan_idx);
        Size apex_peak_idx(m_it->peak_idx);

        if (peak_visited[spec_offsets[apex_scan_idx] + apex_peak_idx]) { continue; }

        Peak2D apex_peak;
        apex_peak.setRT(work_exp[apex_scan_idx].getRT());
        apex_peak.setMZ(work_exp[apex_scan_idx][apex_peak_idx].getMZ());
        apex_peak.setIntensity(work_exp[apex_scan_idx][apex_peak_idx].getIntensity());

        std::list<PeakType> current_trace;
        current_trace.push_back(apex_peak);
        std::vector<double> fwhms_mz; // peak-FWHM meta values of collected peaks
        std::vector<double> fwhms_im; // peak-FWHM ion mobility peak FWHM of collected peaks

        // Initialization for the iterative version of weighted m/z mean calculation
        double centroid_mz(apex_peak.getMZ());
        double prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
        double prev_denom(apex_peak.getIntensity());

        updateIterativeWeightedMean_(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

        // Initialization for the iterative version of weighted ion mobility mean calculation
        double centroid_im(-1);
        double prev_counter_im(-1);
        double prev_denom_im(-1);
        if (has_centroid_im_)
        {
          centroid_im = work_exp[apex_scan_idx].getFloatDataArrays()[ion_mobility_idx_][apex_peak_idx];
          prev_counter_im = apex_peak.getIntensity() * centroid_im;
          prev_denom_im = apex_peak.getIntensity();
          updateIterativeWeightedMean_(work_exp[apex_scan_idx].getFloatDataArrays()[ion_mobility_idx_][apex_peak_idx],
                                        apex_peak.getIntensity(), centroid_im, prev_counter_im, prev_denom_im);
        }

        std::vector<std::pair<Size, Size>> gathered_idx;
        gathered_idx.emplace_back(apex_scan_idx, apex_peak_idx);
        
        if (has_fwhm_mz_) { fwhms_mz.push_back(work_exp[apex_scan_idx].getFloatDataArrays()[fwhm_meta_idx_][apex_peak_idx]); }
        if (has_fwhm_im_) { fwhms_im.push_back(work_exp[apex_scan_idx].getFloatDataArrays()[im_fwhm_idx_][apex_peak_idx]); }

        TraceExtensionState down_state, up_state;
        Size trace_down_idx(apex_scan_idx);
        Size trace_up_idx(apex_scan_idx);

        double ftl_sd((centroid_mz / 1e6) * mass_error_ppm_);
        double intensity_so_far(apex_peak.getIntensity());

        while ((trace_down_idx > 0 && down_state.active) || (trace_up_idx < work_exp.size() - 1 && up_state.active))
        {
          // *********************************************************** //
          // Step 2.1 MOVE DOWN in RT dim
          // *********************************************************** //
          if (trace_down_idx > 0 && down_state.active)
          {
            const MSSpectrum& spec_trace_down = work_exp[trace_down_idx - 1];
            
            // Only process spectra that contain peaks
            if (!spec_trace_down.empty())
            {
              PeakCandidate candidate = findBestPeak_(spec_trace_down, centroid_mz, ftl_sd, centroid_im);
              
              if (isPeakAcceptable_(candidate, centroid_mz, ftl_sd, centroid_im, trace_down_idx - 1, spec_offsets, peak_visited))
              {
                processPeak_(candidate, spec_trace_down, current_trace, gathered_idx,
                            fwhms_mz, fwhms_im, centroid_mz, centroid_im,
                            prev_counter, prev_denom, prev_counter_im, prev_denom_im,
                            ftl_sd, intensity_so_far, trace_down_idx - 1, false);
                ++down_state.hitting_peak_count;
                down_state.consecutive_missed = 0;
              }
              else
              {
                ++down_state.consecutive_missed;
              }
            }
            // Empty spectra don't affect termination counters
            
            // Move to next spectrum regardless of whether peak was found
            --trace_down_idx;
            ++down_state.scan_counter;
            
            // Apply termination criteria based on user configuration
            if (trace_termination_criterion_ == "outlier")
            {
              if (down_state.consecutive_missed > trace_termination_outliers_)
              {
                down_state.active = false;
              }
            }
            else if (trace_termination_criterion_ == "sample_rate")
            {
              Size min_scans_to_consider = 5;
              Size total_hits = down_state.hitting_peak_count + up_state.hitting_peak_count + 1; // +1 for apex
              Size total_scans = down_state.scan_counter + up_state.scan_counter + 1;
              double current_sample_rate = static_cast<double>(total_hits) / static_cast<double>(total_scans);
              
              if (down_state.scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
              {
                down_state.active = false;
              }
            }
          }
          // *********************************************************** //
          // Step 2.2 MOVE UP in RT dim
          // *********************************************************** //
          if (trace_up_idx < work_exp.size() - 1 && up_state.active)
          {
            const MSSpectrum& spec_trace_up = work_exp[trace_up_idx + 1];
            
            // Only process spectra that contain peaks
            if (!spec_trace_up.empty())
            {
              PeakCandidate candidate = findBestPeak_(spec_trace_up, centroid_mz, ftl_sd, centroid_im);
              
              if (isPeakAcceptable_(candidate, centroid_mz, ftl_sd, centroid_im, trace_up_idx + 1, spec_offsets, peak_visited))
              {
                processPeak_(candidate, spec_trace_up, current_trace, gathered_idx,
                            fwhms_mz, fwhms_im, centroid_mz, centroid_im,
                            prev_counter, prev_denom, prev_counter_im, prev_denom_im,
                            ftl_sd, intensity_so_far, trace_up_idx + 1, true);
                ++up_state.hitting_peak_count;
                up_state.consecutive_missed = 0;
              }
              else
              {
                ++up_state.consecutive_missed;
              }
            }
            // Empty spectra don't affect termination counters
            
            // Move to next spectrum regardless of whether peak was found
            ++trace_up_idx;
            ++up_state.scan_counter;
            
            // Apply termination criteria based on user configuration
            if (trace_termination_criterion_ == "outlier")
            {
              if (up_state.consecutive_missed > trace_termination_outliers_)
              {
                up_state.active = false;
              }
            }
            else if (trace_termination_criterion_ == "sample_rate")
            {
              Size min_scans_to_consider = 5;
              Size total_hits = down_state.hitting_peak_count + up_state.hitting_peak_count + 1; // +1 for apex
              Size total_scans = down_state.scan_counter + up_state.scan_counter + 1;
              double current_sample_rate = static_cast<double>(total_hits) / static_cast<double>(total_scans);
              
              if (up_state.scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
              {
                up_state.active = false;
              }
            }
          }
        }

        // *********************************************************** //
        // Step 2.3 check if minimum length and quality of mass trace criteria are met
        // *********************************************************** //
        Size total_scans = down_state.scan_counter + up_state.scan_counter + 1;
        if (isTraceValid_(current_trace, total_scans, down_state.consecutive_missed, up_state.consecutive_missed))
        {
          // mark all peaks as visited
          for (Size i = 0; i < gathered_idx.size(); ++i)
          {
            peak_visited[spec_offsets[gathered_idx[i].first] + gathered_idx[i].second] = true;
          }

          // create new MassTrace object and store collected peaks from list current_trace
          MassTrace new_trace(current_trace);
          new_trace.updateWeightedMeanRT();
          new_trace.updateWeightedMeanMZ();
          if (!fwhms_mz.empty()) { new_trace.fwhm_mz_avg = Math::median(fwhms_mz.begin(), fwhms_mz.end()); }
          if (!fwhms_im.empty()) { new_trace.fwhm_im_avg = Math::median(fwhms_im.begin(), fwhms_im.end()); }
          if (has_centroid_im_) { new_trace.setCentroidIM(centroid_im); }

          new_trace.setQuantMethod(quant_method_);
          new_trace.updateWeightedMZsd();
          new_trace.setLabel("T" + String(trace_number));
          ++trace_number;

          found_masstraces.push_back(new_trace);

          peaks_detected += new_trace.getSize();
          this->setProgress(peaks_detected);

          // check if we already reached the (optional) maximum number of traces
          if (max_traces > 0 && found_masstraces.size() == max_traces) { break; }
        }
      }
      this->endProgress();
    }


    void MassTraceDetection::updateMembers_()
    {
      mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
      noise_threshold_int_ = (double)param_.getValue("noise_threshold_int");
      chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
      ion_mobility_tolerance_ = (double)param_.getValue("ion_mobility_tolerance");
      quant_method_ = MassTrace::getQuantMethod((String)param_.getValue("quant_method").toString());

      trace_termination_criterion_ = (String)param_.getValue("trace_termination_criterion").toString();
      trace_termination_outliers_ = (Size)param_.getValue("trace_termination_outliers");
      min_sample_rate_ = (double)param_.getValue("min_sample_rate");
      min_trace_length_ = (double)param_.getValue("min_trace_length");
      max_trace_length_ = (double)param_.getValue("max_trace_length");
      reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
    }

}
