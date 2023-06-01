// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Erhan Kenar, Holger Franken, Tristan Aretz, Manuel Zschaebitz $
// --------------------------------------------------------------------------

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

bool isInRange(int value, int lowerBound, int upperBound) {
    return value >= lowerBound && value <= upperBound;
}

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

    MassTraceDetection::Apex::Apex(PeakMap& map, const Size scan_idx, const Size peak_idx):
      map_(map),
      scan_idx_(scan_idx),
      peak_idx_(peak_idx)
    {}

    double MassTraceDetection::Apex::getMZ() const
    {
      return map_.get()[scan_idx_][peak_idx_].getMZ();
    }

    double MassTraceDetection::Apex::getRT() const
    {
      return map_.get()[scan_idx_].getRT();
    }

    double MassTraceDetection::Apex::getIntensity() const
    {
      return map_.get()[scan_idx_][peak_idx_].getIntensity();
    }
    
    MassTraceDetection::NextIndex::NextIndex(const std::vector<Apex>& data, const Size total_peak_count, const std::vector<Size>& spec_offsets, const double mass_error_ppm):
      data_(data),
      spec_offsets_(spec_offsets),
      peak_visited_(total_peak_count),
      current_Apex_(0),
      mass_error_ppm_(mass_error_ppm)
    {
    }

    void MassTraceDetection::NextIndex::setNumberOfThreads(const Size thread_num)
    {
      lock_list_.resize(thread_num);
    }

    /// checks if for the current Apex a lock exist in the lock_list, returns a bool
    bool MassTraceDetection::NextIndex::isConflictingApex(const MassTraceDetection::Apex a) const
    {
      for (const double & i : lock_list_)
      {
        if (isInRange(i, a.getMZ() - (3 * findOffset_(a.getMZ(), mass_error_ppm_)), a.getMZ() + (3 * findOffset_(a.getMZ(), mass_error_ppm_))))
        {
          return true;
        }
      }
      return false;
    }

    /// checks if a peak with a certain scan and peak index was already visited, returns a bool
    bool MassTraceDetection::NextIndex::isVisited(const Size scan_idx, const Size peak_idx) const
    {
      return peak_visited_[spec_offsets_[scan_idx] +  peak_idx];
    }
    
    /// returns the index of the Next Free Index, ignoring peaks that might be marked visited from a thread that is currently working
    Size MassTraceDetection::NextIndex::getNextFreeIndex()
    {
      Size result{};
      //isConflictingApex muss anders behandelt werden als isVisited
      #ifdef _OPENMP
      #pragma omp critical (look_for_lock)
      #endif
      {
        while((current_Apex_ < data_.size()) &&
              (*this).isVisited(data_[current_Apex_].scan_idx_, data_[current_Apex_].peak_idx_))
        {
          ++current_Apex_;
        }
        if (current_Apex_ < data_.size())
        {
          while((*this).isConflictingApex(data_[current_Apex_]))
          {
            #ifdef _OPENMP
            int threadnum = omp_get_thread_num();
            #else
            int threadnum = 0;
            #endif
          }
          #ifdef _OPENMP
          int threadnum = omp_get_thread_num();
          #else
          int threadnum = 0;
          #endif
          lock_list_[threadnum] = data_[current_Apex_].getMZ();
        }
        result = current_Apex_;
        ++current_Apex_;
      }
      return result; //return apex 
    }

    /// removes the lock for the Apex of the current thread
    void MassTraceDetection::NextIndex::setApexAsProcessed()
    {
      #ifdef _OPENMP
      lock_list_[omp_get_thread_num()] = 0;
      #else
      lock_list_[0] = 0;
      #endif
    }
    
    /// removes the lock for the Apex of the current thread and marks the gathered peaks as visited
    void MassTraceDetection::NextIndex::setApexAsProcessed(const std::vector<std::pair<Size, Size> >& gathered_idx)
    {
      #ifdef _OPENMP
      #pragma omp critical (remove_lock_from_vec_2)
      #endif
      {
        for (Size i = 0; i < gathered_idx.size(); ++i)
        {
          peak_visited_[spec_offsets_[gathered_idx[i].first] +  gathered_idx[i].second] = true;
        } // sequence important because a lock may only be removed once found peaks have been marked as visited
        (*this).setApexAsProcessed();
      }
    }

    void MassTraceDetection::updateIterativeWeightedMeanMZ(const double added_mz,
                                                           const double added_int,
                                                           double& centroid_mz,
                                                           double& prev_counter,
                                                           double& prev_denom)
    {
      const double nominater(1 + (added_int * added_mz) / prev_counter);
      const double denominator(1 + (added_int) / prev_denom);

      centroid_mz *= (nominater / denominator);
      prev_counter *= nominater;
      prev_denom *= denominator;

      return;
    }


    void MassTraceDetection::run(PeakMap::ConstAreaIterator& begin,
                                 PeakMap::ConstAreaIterator& end,
                                 std::vector<MassTrace>& found_masstraces)
    {
      if (begin == end)
      {
        return;
      }

      PeakMap map;
      MSSpectrum current_spectrum;

      while (begin != end)
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
        ++begin;
      }

      map.addSpectrum(std::move(current_spectrum));

      run(map, found_masstraces);
    }


    void updateWeightedSDEstimateRobust(const PeakType& p, const double& mean_t1, double& sd_t, double& last_weights_sum)
    {
      double denom = std::sqrt
      (
        std::exp(std::log(last_weights_sum) + 2 * std::log(sd_t)) +
        std::exp(std::log(p.getIntensity()) + 2 * std::log(std::abs(p.getMZ() - mean_t1)))
      );
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
      std::vector<Size> spec_offsets;

      Size total_peak_count(0);
      Size spectra_count(0);

      spec_offsets.push_back(0);

      // *********************************************************** //
      //  Step 1: Detecting potential chromatographic apices
      // *********************************************************** //

      for (const MSSpectrum& it : input_exp)
      {
        // check if this is a MS1 survey scan
        if (it.getMSLevel() != 1 || it.empty())
        {
          continue;
        }
        std::vector<Size> indices_passing;
        
        for(Size peak_idx = 0; peak_idx < it.size(); ++peak_idx)
        {
          if (it[peak_idx].getIntensity() > noise_threshold_int_)
          {
            // Assume that noise_threshold_int_ contains the noise level of the
            // data and we want to be chrom_peak_snr times above the noise level
            // --> add this peak as possible chromatographic apex
            if (it[peak_idx].getIntensity() > chrom_peak_snr_ * noise_threshold_int_)
            {
              chrom_apices.emplace_back(work_exp, spectra_count, indices_passing.size());
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
                [&chrom_apices](const Apex & a,
                    const Apex & b) -> bool
      { 
        return a.getIntensity() > b.getIntensity();
      });

      // *********************************************************************
      // Step 2: start extending mass traces beginning with the apex peak (go
      // through all peaks in order of decreasing intensity)
      // *********************************************************************
      run_(chrom_apices, total_peak_count, work_exp, spec_offsets, found_masstraces, max_traces);

      return;
    } // end of MassTraceDetection::run

    /// calculates the offset depending off the peaks mz and the mass_error_ppm, returns double
    double MassTraceDetection::findOffset_(const double centroid_mz, const double mass_error_ppm_)
    {
      return (3 * Math::ppmToMass(mass_error_ppm_,centroid_mz));
    }


    bool MassTraceDetection::checkFWHMMetaData_(const PeakMap& work_exp)
    {
      Size fwhm_meta_count(0);

      for (Size i = 0; (i < work_exp.size()); ++i)
      {
        if ((!work_exp[i].getFloatDataArrays().empty()) &&
            (work_exp[i].getFloatDataArrays()[0].getName() == "FWHM_ppm"))
        {
          if (work_exp[i].getFloatDataArrays()[0].size() != work_exp[i].size())
          { // float data should always have the same size as the corresponding array
            throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, work_exp[i].size());
          }
          ++fwhm_meta_count;
        }
      }
      if (fwhm_meta_count == 0) 
      {
        return false;
      }
      else if (fwhm_meta_count == work_exp.size())
      {
        return true;
      }
      else 
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("FWHM meta arrays are expected to be missing or present for all MS spectra [") + fwhm_meta_count + "/" + work_exp.size() + "].");
      }
    }

    void MassTraceDetection::run_(std::vector<Apex>& chrom_apices,
                                  const Size total_peak_count,
                                  const PeakMap& work_exp,
                                  const std::vector<Size>& spec_offsets,
                                  std::vector<MassTrace>& found_masstraces,
                                  const Size max_traces)
    {         
      ///check for FWHM meta data & corrupted data
      bool fwhm_meta_idx = checkFWHMMetaData_(work_exp);

      Size trace_number(1);
      Size peaks_detected(0);
      
      this->startProgress(0, total_peak_count, "mass trace detection");
      
      NextIndex nextIndex(chrom_apices, total_peak_count, spec_offsets, mass_error_ppm_);
      #ifdef _OPENMP
      #pragma omp parallel 
      #endif
      {
        #ifdef _OPENMP
        {
          #pragma omp single
          nextIndex.setNumberOfThreads(omp_get_num_threads()); // resize lock_list_ to the size of used threads
        }
        #else 
          nextIndex.setNumberOfThreads(1); 
        #endif
      } 
    
      // start parallel region 
      #ifdef _OPENMP
      #pragma omp parallel
      #endif
      while((max_traces < 1) || (found_masstraces.size() < max_traces)) // only here for the max_traces threshold
      {
        Size index{};

     
        index = nextIndex.getNextFreeIndex();

        if(index >= chrom_apices.size()) break; // break while loop if index 

        if(nextIndex.isVisited(chrom_apices[index].scan_idx_, chrom_apices[index].peak_idx_)) 
        {
          nextIndex.setApexAsProcessed();
          continue;
        } // Important because a previous free index from getNextFreeIndex() can be marked as visited later setApexAsProcessed, especially if it was locked due to the mz and rt range
        
        Peak2D apex_peak;
        apex_peak.setRT(chrom_apices[index].getRT());
        apex_peak.setMZ(chrom_apices[index].getMZ());
        apex_peak.setIntensity(chrom_apices[index].getIntensity());

        Size trace_up_idx(chrom_apices[index].scan_idx_);
        Size trace_down_idx(chrom_apices[index].scan_idx_);

        std::deque<PeakType> current_trace;
        current_trace.push_back(apex_peak);
        std::vector<double> fwhms_mz; // peak-FWHM meta values of collected peaks

        // Initialization for the iterative version of weighted m/z mean calculation
        double centroid_mz(apex_peak.getMZ());
        double prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
        double prev_denom(apex_peak.getIntensity());

        updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

        std::vector<std::pair<Size, Size> > gathered_idx;
        gathered_idx.emplace_back(chrom_apices[index].scan_idx_,chrom_apices[index].peak_idx_);
        if (fwhm_meta_idx)
        {
          fwhms_mz.push_back(work_exp[chrom_apices[index].scan_idx_].getFloatDataArrays()[fwhm_meta_idx][chrom_apices[index].peak_idx_]);
        }

        Size up_hitting_peak(0), down_hitting_peak(0);
        Size up_scan_counter(0), down_scan_counter(0);

        bool toggle_up = true, toggle_down = true;

        Size conseq_missed_peak_up(0), conseq_missed_peak_down(0);
        Size max_consecutive_missing(trace_termination_outliers_);

        double current_sample_rate(1.0);
        // Size min_scans_to_consider(std::floor((min_sample_rate_ /2)*10));
        Size min_scans_to_consider(5);

        double ftl_sd(Math::ppmToMass(mass_error_ppm_, centroid_mz));
        double intensity_so_far(apex_peak.getIntensity());

        while (((trace_down_idx > 0) && toggle_down) ||
              ((trace_up_idx < work_exp.size() - 1) && toggle_up) 
              )
        {
  
          // *********************************************************** //
          // Step 2.1 MOVE DOWN in RT dim
          // *********************************************************** //
          if ((trace_down_idx > 0) && 
              toggle_down)
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
                  (!nextIndex.isVisited(trace_down_idx - 1, next_down_peak_idx))
                      )
              {

                Peak2D next_peak;
                next_peak.setRT(spec_trace_down.getRT());
                next_peak.setMZ(next_down_peak_mz);
                next_peak.setIntensity(next_down_peak_int);

                current_trace.push_front(next_peak);
                // FWHM average
                if (fwhm_meta_idx)
                {
                  fwhms_mz.push_back(spec_trace_down.getFloatDataArrays()[0][next_down_peak_idx]);
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
                toggle_down = false;
              }
            }
          }

          // *********************************************************** //
          // Step 2.2 MOVE UP in RT dim
          // *********************************************************** //
          if ((trace_up_idx < work_exp.size() - 1) && 
              toggle_up)
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
                  (!nextIndex.isVisited(trace_up_idx + 1, next_up_peak_idx)))
                  // !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
              {
                Peak2D next_peak;
                next_peak.setRT(spec_trace_up.getRT());
                next_peak.setMZ(next_up_peak_mz);
                next_peak.setIntensity(next_up_peak_int);

                current_trace.push_back(next_peak);
                if (fwhm_meta_idx)
                {
                  fwhms_mz.push_back(spec_trace_up.getFloatDataArrays()[0][next_up_peak_idx]);
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
                toggle_up = false;
              }
            }
          }
        }



        double num_scans(down_scan_counter + up_scan_counter + 1 - conseq_missed_peak_down - conseq_missed_peak_up);

        double mt_quality((double)current_trace.size() / (double)num_scans);

        double rt_range(std::fabs(current_trace.rbegin()->getRT() - current_trace.begin()->getRT()));

        // *********************************************************** //
        // Step 2.3 check if minimum length and quality of mass trace criteria are met
        // *********************************************************** //
        bool max_trace_criteria = (max_trace_length_ < 0.0 || rt_range < max_trace_length_);
        if (rt_range >= min_trace_length_ && max_trace_criteria && mt_quality >= min_sample_rate_ )
        {
              nextIndex.setApexAsProcessed(gathered_idx); // free lock and gathered peaks as visited
              MassTrace new_trace(current_trace.begin(), current_trace.end());
              new_trace.updateWeightedMeanRT();
              new_trace.updateWeightedMeanMZ();

              if (!fwhms_mz.empty())
              {
                new_trace.fwhm_mz_avg = Math::median(fwhms_mz.begin(), fwhms_mz.end());
              }
              new_trace.setQuantMethod(quant_method_);
              //new_trace.setCentroidSD(ftl_sd);
              new_trace.updateWeightedMZsd();

              // critical for safe push_back, also safe addtion of global variables
              #ifdef _OPENMP
              #pragma omp critical (add_trace)
              #endif
              {
                new_trace.setLabel("T" + String(trace_number));
                ++trace_number;
                found_masstraces.push_back(new_trace);
                peaks_detected += new_trace.getSize();
                
                this->setProgress(peaks_detected);
              }
          } 
          else
          {
            nextIndex.setApexAsProcessed();
          } // else case to free lock if trace criteriums not matched, Apex processsed wihtout the collected peaks as visited 
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