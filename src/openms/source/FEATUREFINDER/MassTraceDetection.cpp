// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <algorithm>
#include <boost/pending/disjoint_sets.hpp>
#include <exception>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <OpenMS/KERNEL/SpectrumHelper.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

namespace OpenMS
{

    MassTraceDetection::MassTraceDetection() :
            DefaultParamHandler("MassTraceDetection"), ProgressLogger()
    {
      defaults_.setValue(
        "mass_error_ppm", 20.0,
        "Initial m/z standard deviation in ppm. Adaptive searches use three standard deviations, additionally bounded by max_apex_deviation_ppm.");
      defaults_.setValue(
        "max_apex_deviation_ppm", 20.0,
        "Hard m/z deviation in ppm from the fixed seed/apex, independent of centroid and variance updates. Set to 0 for unbounded legacy behavior.",
        {"advanced"});
      defaults_.setMinFloat("max_apex_deviation_ppm", 0.0);
      defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are removed as noise.");
      defaults_.setValue("chrom_peak_snr", 3.0, "Minimum intensity above noise_threshold_int (signal-to-noise) a peak should have to be considered an apex.");
      defaults_.setValue("ion_mobility_tolerance", 0.01, "Allowed ion mobility deviation (in 1/K0 units, e.g., 0.01-0.05). For CCS data, use larger values (e.g., 2-10 square angstroms).");

      defaults_.setValue("reestimate_mt_sd", "true", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
      defaults_.setValidStrings("reestimate_mt_sd", {"true","false"});

      defaults_.setValue("quant_method",std::string(MassTrace::names_of_quantmethod[0]), "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace.");
      defaults_.setValidStrings("quant_method", std::vector<std::string>(MassTrace::names_of_quantmethod, MassTrace::names_of_quantmethod +(int)MassTrace::SIZE_OF_MT_QUANTMETHOD));

      // advanced parameters
      defaults_.setValue("trace_termination_criterion", "outlier", "Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold.", {"advanced"});
      defaults_.setValidStrings("trace_termination_criterion", {"outlier","sample_rate"});
      defaults_.setValue("trace_termination_outliers", 5,
                         "Maximum consecutive nonempty working spectra without an acceptable peak; extension stops on the next miss. Empty working "
                         "spectra do not count as misses.",
                         {"advanced"});

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
      auto validate_meta_array = [&](const std::string& name, int idx) {
        if (idx == -1) return;

        Size valid_count = 0;
        for (const auto& spec : spectra)
        {
          const auto& fda = spec.getFloatDataArrays();
          if (!fda.empty() && idx < (int)fda.size() && fda[idx].getName() == name)
          {
            if (fda[idx].size() != spec.size())
            {
              throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, spec.size(), "FloatDataArray size does not match spectrum size");
            }
            ++valid_count;
          }
        }

        if (valid_count > 0 && valid_count != spectra.size())
        {
          throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                name + " meta arrays must be consistently present or absent across all MS spectra ["
                                                  + StringUtils::toStr(valid_count) + "/" + StringUtils::toStr(spectra.size()) + "].");
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

      OPENMS_PRECONDITION(std::is_sorted(input_exp.begin(), input_exp.end(), PeakMap::SpectrumType::RTLess()),
       "MassTraceDetection::run() requires input spectra to be sorted by retention time. Call MSExperiment::sortSpectra() first.");

      // Check IM unit and warn if CCS data with small tolerance (single-threaded check at algorithm start)
      for (const auto& spec : input_exp)
      {
        if (spec.containsIMData())
        {
          const auto [im_data_index, im_unit] = spec.getIMData();
          if (im_unit == DriftTimeUnit::CCS && ion_mobility_tolerance_ < 1.0)
          {
            OPENMS_LOG_WARN << "Warning: Ion mobility data is in CCS units (square angstroms), but "
                            << "ion_mobility_tolerance = " << ion_mobility_tolerance_
                            << " appears to be set for 1/K0 data. "
                            << "For CCS data, consider using larger values (e.g., 2-10)." << '\n';
          }
          break; // Found IM data, no need to check further spectra
        }
      }

      // gather all peaks that are potential chromatographic peak apices
      //   - use work_exp for actual work (remove peaks below noise threshold)
      //   - store potential apices in chrom_apices
      PeakMap work_exp;
      std::vector<Apex> chrom_apices;

      Size total_peak_count(0);
      std::vector<Size> spec_offsets;
      spec_offsets.reserve(input_exp.size() + 1); // Pre-allocate for performance
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
                                      "Input map consists of too few MS1 spectra (less than 3!). Aborting...",StringUtils::toStr(spectra_count));
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

    MassTraceDetection::PeakCandidate
    MassTraceDetection::findBestPeak_(const MSSpectrum& spectrum, double centroid_mz, double ftl_sd, double apex_mz, double centroid_im) const
    {
      PeakCandidate candidate;
      
      if (spectrum.empty())
      {
        return candidate;
      }

      double right_bound = centroid_mz + 3 * ftl_sd;
      double left_bound = centroid_mz - 3 * ftl_sd;
      if (max_apex_deviation_ppm_ > 0.0)
      {
        const double apex_tolerance = apex_mz * max_apex_deviation_ppm_ * 1e-6;
        left_bound = std::max(left_bound, apex_mz - apex_tolerance);
        right_bound = std::min(right_bound, apex_mz + apex_tolerance);
      }
      if (left_bound > right_bound) { return candidate; }

      if (!has_centroid_im_)
      {
        // Standard LC-MS data: find peak closest to target m/z
        const double target = std::clamp(centroid_mz, left_bound, right_bound);
        const Int index = spectrum.findNearest(target, target - left_bound, right_bound - target);
        if (index < 0) { return candidate; }
        candidate.idx = static_cast<Size>(index);
        candidate.mz = spectrum[candidate.idx].getMZ();
        candidate.intensity = spectrum[candidate.idx].getIntensity();
        candidate.found = candidate.mz >= left_bound && candidate.mz <= right_bound;
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

    bool MassTraceDetection::isPeakAcceptable_(const PeakCandidate& candidate,
                                               double centroid_mz,
                                               double ftl_sd,
                                               double centroid_im,
                                               Size spectrum_idx,
                                               const std::vector<Size>& spec_offsets,
                                               const std::vector<unsigned char>& peak_visited) const
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

    std::vector<std::vector<MassTraceDetection::Apex>> MassTraceDetection::partitionApices_(const std::vector<Apex>& chrom_apices,
                                                                                            const PeakMap& work_exp,
                                                                                            const std::vector<Size>& spec_offsets,
                                                                                            Size peak_count) const
    {
      const double gate = 2.0 * max_apex_deviation_ppm_ * 1e-6 + 16.0 * std::numeric_limits<double>::epsilon();
      if (gate >= 1.0) { return {chrom_apices}; }

      for (const auto& spectrum : work_exp)
      {
        for (const auto& peak : spectrum)
        {
          if (! std::isfinite(peak.getMZ()) || peak.getMZ() <= 0.0)
          {
            OPENMS_LOG_WARN << "Mass-trace component parallelism requires finite, positive m/z values; using serial extraction.\n";
            return {chrom_apices};
          }
        }
      }

      boost::disjoint_sets_with_storage<> connected(peak_count);
      std::vector<Size> nonempty_scans;
      nonempty_scans.reserve(work_exp.size());

      for (Size scan = 0; scan < work_exp.size(); ++scan)
      {
        const auto& spectrum = work_exp[scan];
        if (spectrum.empty()) { continue; }

        for (Size peak = 1; peak < spectrum.size(); ++peak)
        {
          const double lower = std::nextafter(spectrum[peak].getMZ() * (1.0 - gate), -std::numeric_limits<double>::infinity());
          if (spectrum[peak - 1].getMZ() >= lower) { connected.union_set(spec_offsets[scan] + peak - 1, spec_offsets[scan] + peak); }
        }

        // Empty scans do not count as misses. A hit can follow outliers + 1 nonempty scans.
        const Size first = nonempty_scans.size() > trace_termination_outliers_ ? nonempty_scans.size() - trace_termination_outliers_ - 1 : 0;
        for (Size previous = first; previous < nonempty_scans.size(); ++previous)
        {
          const Size other_scan = nonempty_scans[previous];
          const auto& other = work_exp[other_scan];
          Size left = 0, right = 0;
          for (Size peak = 0; peak < spectrum.size(); ++peak)
          {
            const double mz = spectrum[peak].getMZ();
            // Symmetric 2H gate relative to max(m/z); outward rounding includes cap endpoints.
            const double lower = std::nextafter(mz * (1.0 - gate), -std::numeric_limits<double>::infinity());
            const double upper = std::nextafter(mz / (1.0 - gate), std::numeric_limits<double>::infinity());
            while (left < other.size() && other[left].getMZ() < lower)
            {
              ++left;
            }
            right = std::max(right, left);
            while (right < other.size() && other[right].getMZ() <= upper)
            {
              ++right;
            }
            for (Size candidate = left; candidate < right; ++candidate)
            {
              connected.union_set(spec_offsets[scan] + peak, spec_offsets[other_scan] + candidate);
            }
          }
        }
        nonempty_scans.push_back(scan);
      }

      std::vector<std::vector<Apex>> components;
      std::vector<Size> component_index(peak_count, peak_count);
      for (const auto& apex : chrom_apices)
      {
        const Size root = connected.find_set(spec_offsets[apex.scan_idx] + apex.peak_idx);
        if (component_index[root] == peak_count)
        {
          component_index[root] = components.size();
          components.emplace_back();
        }
        components[component_index[root]].push_back(apex);
      }
      std::stable_sort(components.begin(), components.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
      return components;
    }

    void MassTraceDetection::run_(const std::vector<Apex>& chrom_apices,
                                  const Size total_peak_count,
                                  const PeakMap& work_exp,
                                  const std::vector<Size>& spec_offsets,
                                  std::vector<MassTrace>& found_masstraces,
                                  const Size max_traces)
    {
      // Unlike vector<bool>, distinct bytes can be written by independent components.
      std::vector<unsigned char> peak_visited(total_peak_count, 0);

      // Detect ion mobility and FWHM metadata arrays in the dataset
      fwhm_meta_idx_ = ion_mobility_idx_ = im_fwhm_idx_ = -1;
      has_fwhm_mz_ = has_centroid_im_ = has_fwhm_im_ = false;
      getIMIndices_(work_exp,
                    fwhm_meta_idx_, has_fwhm_mz_,
                    ion_mobility_idx_, has_centroid_im_,
                    im_fwhm_idx_, has_fwhm_im_);

      this->startProgress(0, total_peak_count, "mass trace detection");
      bool parallel = false;
#ifdef _OPENMP
      parallel = omp_get_max_threads() > 1 && omp_in_parallel() == 0 && max_apex_deviation_ppm_ > 0.0 && trace_termination_criterion_ == OUTLIER
                 && max_traces == 0;
#endif
      if (parallel && ! chrom_apices.empty())
      {
        const auto components = partitionApices_(chrom_apices, work_exp, spec_offsets, total_peak_count);
        std::vector<std::vector<MassTrace>> results(components.size());
        std::vector<std::vector<Apex>> accepted_apices(components.size());
        std::vector<std::exception_ptr> errors(components.size());
        Size peaks_detected = 0;
#pragma omp parallel for schedule(dynamic, 1)
        for (SignedSize component = 0; component < static_cast<SignedSize>(components.size()); ++component)
        {
          try
          {
            runComponent_(components[component], work_exp, spec_offsets, peak_visited, results[component], &accepted_apices[component], 0, false);
          }
          catch (...)
          {
            errors[component] = std::current_exception();
          }
          Size collected = 0;
          for (const auto& trace : results[component])
          {
            collected += trace.getSize();
          }
#pragma omp critical(MassTraceDetection_progress)
          {
            peaks_detected += collected;
            this->setProgress(peaks_detected);
          }
        }
        for (const auto& error : errors)
        {
          if (error)
          {
            this->endProgress();
            std::rethrow_exception(error);
          }
        }
        Size trace_count = 0;
        for (const auto& result : results)
        {
          trace_count += result.size();
        }
        // MassTrace has no move operations; sort handles rather than copying all trace peaks.
        std::vector<std::pair<Apex, const MassTrace*>> traces;
        traces.reserve(trace_count);
        for (Size component = 0; component < results.size(); ++component)
        {
          for (Size trace = 0; trace < results[component].size(); ++trace)
          {
            traces.emplace_back(accepted_apices[component][trace], &results[component][trace]);
          }
        }
        std::sort(traces.begin(), traces.end(), [](const auto& a, const auto& b) {
          if (a.first.intensity != b.first.intensity) { return a.first.intensity > b.first.intensity; }
          if (a.first.scan_idx != b.first.scan_idx) { return a.first.scan_idx > b.first.scan_idx; }
          return a.first.peak_idx > b.first.peak_idx;
        });
        found_masstraces.reserve(traces.size());
        for (const auto& trace : traces)
        {
          found_masstraces.push_back(*trace.second);
        }
      }
      else
      {
        runComponent_(chrom_apices, work_exp, spec_offsets, peak_visited, found_masstraces, nullptr, max_traces, true);
      }
      for (Size trace = 0; trace < found_masstraces.size(); ++trace)
      {
        found_masstraces[trace].setLabel("T" + StringUtils::toStr(trace + 1));
      }
      this->endProgress();
    }

    void MassTraceDetection::runComponent_(const std::vector<Apex>& chrom_apices,
                                           const PeakMap& work_exp,
                                           const std::vector<Size>& spec_offsets,
                                           std::vector<unsigned char>& peak_visited,
                                           std::vector<MassTrace>& found_masstraces,
                                           std::vector<Apex>* accepted_apices,
                                           Size max_traces,
                                           bool report_progress)
    {
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
              PeakCandidate candidate = findBestPeak_(spec_trace_down, centroid_mz, ftl_sd, apex_peak.getMZ(), centroid_im);

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
            if (trace_termination_criterion_ == OUTLIER)
            {
              if (down_state.consecutive_missed > trace_termination_outliers_)
              {
                down_state.active = false;
              }
            }
            else if (trace_termination_criterion_ == SAMPLE_RATE)
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
              PeakCandidate candidate = findBestPeak_(spec_trace_up, centroid_mz, ftl_sd, apex_peak.getMZ(), centroid_im);

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
            if (trace_termination_criterion_ == OUTLIER)
            {
              if (up_state.consecutive_missed > trace_termination_outliers_)
              {
                up_state.active = false;
              }
            }
            else if (trace_termination_criterion_ == SAMPLE_RATE)
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
          peaks_detected += new_trace.getSize();
          found_masstraces.push_back(new_trace);
          if (accepted_apices != nullptr) { accepted_apices->push_back(*m_it); }
          if (report_progress) { this->setProgress(peaks_detected); }

          // check if we already reached the (optional) maximum number of traces
          if (max_traces > 0 && found_masstraces.size() == max_traces) { break; }
        }
      }
    }


    void MassTraceDetection::updateMembers_()
    {
      mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
      max_apex_deviation_ppm_ = (double)param_.getValue("max_apex_deviation_ppm");
      if (! std::isfinite(max_apex_deviation_ppm_))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "max_apex_deviation_ppm must be finite");
      }
      noise_threshold_int_ = (double)param_.getValue("noise_threshold_int");
      chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
      ion_mobility_tolerance_ = (double)param_.getValue("ion_mobility_tolerance");
      quant_method_ = MassTrace::getQuantMethod(param_.getValue("quant_method").toString());

      std::string criterion_str = param_.getValue("trace_termination_criterion").toString();
      if (criterion_str == "outlier")
      {
        trace_termination_criterion_ = OUTLIER;
      }
      else if (criterion_str == "sample_rate")
      {
        trace_termination_criterion_ = SAMPLE_RATE;
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid trace_termination_criterion: " + criterion_str);
      }
      
      trace_termination_outliers_ = (Size)param_.getValue("trace_termination_outliers");
      min_sample_rate_ = (double)param_.getValue("min_sample_rate");
      min_trace_length_ = (double)param_.getValue("min_trace_length");
      max_trace_length_ = (double)param_.getValue("max_trace_length");
      reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
    }

}
