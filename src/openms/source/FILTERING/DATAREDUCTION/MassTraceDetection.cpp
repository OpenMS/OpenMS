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

#include <omp.h>

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
      return map_[scan_idx_][peak_idx_].getMZ();
    }

    double MassTraceDetection::Apex::getRT() const
    {
      return map_[scan_idx_].getRT();
    }

    double MassTraceDetection::Apex::getIntensity() const
    {
      return map_[scan_idx_][peak_idx_].getIntensity();
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
      lock_list_2_.resize(thread_num);
    }

    bool MassTraceDetection::NextIndex::isConflictingApex(const MassTraceDetection::Apex a) const
    {
      // std::cout << "Here1?\n";
      for (const double & i : lock_list_2_)
      {
        // std::cout << "Here2?\n";
        if (isInRange(i, a.getMZ() - (2 * findOffset_(a.getMZ(), mass_error_ppm_)), a.getMZ() + (2 * findOffset_(a.getMZ(), mass_error_ppm_))))
        {
          // std::cout << "Here3?\n";
          return true;
        }
        // if(i != 0) return true;
      }
      return false;
    }

    bool MassTraceDetection::NextIndex::isVisited(const Size scan_idx, const Size peak_idx) const
    {
      return peak_visited_[spec_offsets_[scan_idx] +  peak_idx];
    }
    
    Size MassTraceDetection::NextIndex::getNextFreeIndex()
    {
      Size result{};
      //isConflictingApex muss anders behandelt werden als isVisited
      #pragma omp critical (look_for_lock)
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
            // std::cout << "Found conflict\n";
            // usleep(1000);
            // std::cout << lock_list_2_.size() << '\n';
          }
          
          // std::cout << "After while\n";
          int threadnum = omp_get_thread_num();
          lock_list_2_[threadnum] = data_[current_Apex_].getMZ();
          // result = current_Apex_;
          // ++current_Apex_;
        }
        result = current_Apex_;
        ++current_Apex_;
      }
      return result; //next apex 
    }

    void MassTraceDetection::NextIndex::setApexAsProcessed()
    {
      // #pragma omp critical (remove_lock_from_vec)
      // {
        // auto it = std::find(lock_list_2_.begin(), lock_list_2_.end(), data_[index].getMZ());
        // if (it != lock_list_2_.end()) 
        // {
        //   // Remove the element
        //   lock_list_2_.erase(it);
        //   std::cout << "Removed conflict\n";
        // }
        lock_list_2_[omp_get_thread_num()] = 0;
        // std::cout << lock_list_2_[0] << ' ' << lock_list_2_[1] << '\n';
      // }
    }
    
    void MassTraceDetection::NextIndex::setApexAsProcessed(const std::vector<std::pair<Size, Size> >& gathered_idx)
    {
      #pragma omp critical (remove_lock_from_vec_2)
      {
        for (Size i = 0; i < gathered_idx.size(); ++i)
        {
          // Size peaks_index = spec_offsets_[gathered_idx[i].first] +  gathered_idx[i].second;
          // peak_visited_.at(spec_offsets_[ gathered_idx[i].first ] +  gathered_idx[i].second) = true;
          peak_visited_[spec_offsets_[gathered_idx[i].first] +  gathered_idx[i].second] = true;
        } // reihenfolge wichtig
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

      ///used variables
      // boost::dynamic_bitset<> peak_visited(total_peak_count);

      Size trace_number(1);
      Size peaks_detected(0);
      
      this->startProgress(0, total_peak_count, "mass trace detection");
      
      NextIndex relevant(chrom_apices, total_peak_count, spec_offsets, mass_error_ppm_); // I think relevant is a weird name
      #pragma omp parallel 
      {
        #pragma omp single
        relevant.setNumberOfThreads(omp_get_num_threads()); // todo: function 'setNumberOfThreads()' implementieren
      } 

      // for(Size i{}; i < relevant.lock_list_2_.size(); ++i)
      // {
      //   relevant.lock_list_2_[i] = 0;
      // }

      // Size index{};
    
      #pragma omp parallel
      while((max_traces < 1) || (found_masstraces.size() < max_traces))
      {
        Size index{};

     
        index = relevant.getNextFreeIndex();

        if(index >= chrom_apices.size()) break;

        if(relevant.isVisited(chrom_apices[index].scan_idx_, chrom_apices[index].peak_idx_)) 
        {
          relevant.setApexAsProcessed();
          continue;
        }
        
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
          fwhms_mz.push_back(work_exp[chrom_apices[index].scan_idx_].getFloatDataArrays()[0][chrom_apices[index].peak_idx_]);
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
              toggle_down &&
              (!relevant.isVisited(relevant.data_[index].scan_idx_, relevant.data_[index].peak_idx_))
              )
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
                  (!relevant.isVisited(trace_down_idx - 1, next_down_peak_idx))
                  // !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx]
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
              toggle_up &&
              (!relevant.isVisited(chrom_apices[index].scan_idx_, chrom_apices[index].peak_idx_))
              )
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
                  (!relevant.isVisited(trace_up_idx + 1, next_up_peak_idx)))
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
            //create new MassTrace object and store collected peaks from list current_trace
              // for (Size i = 0; i < gathered_idx.size(); ++i)
              // {
              //   peak_visited[spec_offsets[gathered_idx[i].first] +  gathered_idx[i].second] = true;
              // }
              relevant.setApexAsProcessed(gathered_idx);
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

              #pragma omp critical (add_trace)
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
            relevant.setApexAsProcessed();
          }
          // if(peaks_detected > chrom_apices.size()) break;
          // check if we already reached the (optional) maximum number of traces
          // if (max_traces > 0 && found_masstraces.size() == max_traces)
          // {
          //   continue;
          // }

      }
        // if (index == 0)
        // {
        // std::cout 
        //   << index << " :\n " 
        //   << trace_number-1 << " : " 
        //   << new_trace.getSize() << " : "
        //   << found_masstraces.size() << "\n"
        //   << found_masstraces[found_masstraces.size() -1].getCentroidMZ() << " : "
        //   << found_masstraces[found_masstraces.size() -1].getCentroidRT() << " : "
        //   << found_masstraces[found_masstraces.size() -1].getCentroidSD() << " \n";
        // } 

    // Size ca = chrom_apices.size()/2;
    // Size fm = found_masstraces.size()/2;
    // std::cout << "Daten zum testen \n"
    //   << "Anzahl an chrom apices: " << chrom_apices.size()
    //   << "\n Anzahl an Traces: " << found_masstraces.size()
    //   << "\n Einzelne Peaks zum ueberpruefen RT|MZ|IN|scan_id|peak_id \n"
    //   << "\nPeak index: 0 " << chrom_apices[0].getRT() << "|" << chrom_apices[0].getMZ() << "|" << chrom_apices[0].getIntensity() << "|" << chrom_apices[0].scan_idx_ << "|" << chrom_apices[0].peak_idx_
    //   << "\nPeak index: " << ca << " " << chrom_apices[ca].getRT() << "|" << chrom_apices[ca].getMZ() << "|" << chrom_apices[ca].getIntensity() << "|" << chrom_apices[ca].scan_idx_ << "|" << chrom_apices[ca].peak_idx_
    //   << "\nPeak index: " << chrom_apices.size()-1 << " " << chrom_apices[chrom_apices.size()-1].getRT() << "|" << chrom_apices[chrom_apices.size()-1].getMZ() << "|" << chrom_apices[chrom_apices.size()-1].getIntensity() << "|" << chrom_apices[chrom_apices.size()-1].scan_idx_ << "|" << chrom_apices[chrom_apices.size()-1].peak_idx_
    //   << "\n Einzelne Traces zum ueberpruefen RT|MZ|SD|Size|Label \n"
    //   << "\nPeak index: 0 " << found_masstraces[0].getIntensity(false) << "|" << found_masstraces[0].getCentroidMZ() << "|" << found_masstraces[0].getCentroidSD() << "|" << found_masstraces[0].getSize() << "|" << found_masstraces[0].getLabel()
    //   << "\nPeak index: " << fm << " " << found_masstraces[fm].getCentroidRT() << "|" << found_masstraces[fm].getCentroidMZ() << "|" << found_masstraces[fm].getCentroidSD() << "|" << found_masstraces[fm].getSize() << "|" << found_masstraces[fm].getLabel()
    //   << "\nPeak index: " << found_masstraces.size()-1 << " " << found_masstraces[found_masstraces.size()-1].getCentroidRT() << "|" << found_masstraces[found_masstraces.size()-1].getCentroidMZ() << "|" << found_masstraces[found_masstraces.size()-1].getCentroidSD() << "|" << found_masstraces[found_masstraces.size()-1].getSize() << "|" << found_masstraces[found_masstraces.size()-1].getLabel()
    // << std::endl;
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

// SA1_subset_verysmall 



//       // *********************************************************************
//       // Step 2: start extending mass traces beginning with the apex peak (go
//       // through all peaks in order of decreasing intensity)
//       // *********************************************************************
//       run_(chrom_apices, total_peak_count, work_exp, spec_offsets, found_masstraces, max_traces);

//       return;
//     } // end of MassTraceDetection::run

//     void MassTraceDetection::run_(const std::vector<Apex>& chrom_apices,
//                                   const Size total_peak_count,
//                                   const PeakMap& work_exp,
//                                   const std::vector<Size>& spec_offsets,
//                                   std::vector<MassTrace>& found_masstraces,
//                                   const Size max_traces)
//     {
//       boost::dynamic_bitset<> peak_visited(total_peak_count);
//       Size trace_number(1);

//       // check presence of FWHM meta data
//       int fwhm_meta_idx(-1);
//       Size fwhm_meta_count(0);
//       for (Size i = 0; i < work_exp.size(); ++i)
//       {
//         if (!work_exp[i].getFloatDataArrays().empty() &&
//             work_exp[i].getFloatDataArrays()[0].getName() == "FWHM_ppm")
//         {
//           if (work_exp[i].getFloatDataArrays()[0].size() != work_exp[i].size())
//           { // float data should always have the same size as the corresponding array
//             throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, work_exp[i].size());
//           }
//           fwhm_meta_idx = 0;
//           ++fwhm_meta_count;
//         }
//       }
//       if (fwhm_meta_count > 0 && fwhm_meta_count != work_exp.size())
//       {
//         throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
//                                       String("FWHM meta arrays are expected to be missing or present for all MS spectra [") + fwhm_meta_count + "/" + work_exp.size() + "].");
//       }


//       this->startProgress(0, total_peak_count, "mass trace detection");
//       Size peaks_detected(0);

//       for (auto m_it = chrom_apices.crbegin(); m_it != chrom_apices.crend(); ++m_it)
//       {
//         Size apex_scan_idx(m_it->scan_idx);
//         Size apex_peak_idx(m_it->peak_idx);

//         if (peak_visited[spec_offsets[apex_scan_idx] + apex_peak_idx])
//         {
//           continue;
//         }

//         Peak2D apex_peak;
//         apex_peak.setRT(work_exp[apex_scan_idx].getRT());
//         apex_peak.setMZ(work_exp[apex_scan_idx][apex_peak_idx].getMZ());
//         apex_peak.setIntensity(work_exp[apex_scan_idx][apex_peak_idx].getIntensity());

//         Size trace_up_idx(apex_scan_idx);
//         Size trace_down_idx(apex_scan_idx);

//         std::list<PeakType> current_trace;
//         current_trace.push_back(apex_peak);
//         std::vector<double> fwhms_mz; // peak-FWHM meta values of collected peaks

//         // Initialization for the iterative version of weighted m/z mean calculation
//         double centroid_mz(apex_peak.getMZ());
//         double prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
//         double prev_denom(apex_peak.getIntensity());

//         updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

//         std::vector<std::pair<Size, Size> > gathered_idx;
//         gathered_idx.emplace_back(apex_scan_idx, apex_peak_idx);
//         if (fwhm_meta_idx != -1)
//         {
//           fwhms_mz.push_back(work_exp[apex_scan_idx].getFloatDataArrays()[fwhm_meta_idx][apex_peak_idx]);
//         }

//         Size up_hitting_peak(0), down_hitting_peak(0);
//         Size up_scan_counter(0), down_scan_counter(0);

//         bool toggle_up = true, toggle_down = true;

//         Size conseq_missed_peak_up(0), conseq_missed_peak_down(0);
//         Size max_consecutive_missing(trace_termination_outliers_);

//         double current_sample_rate(1.0);
//         // Size min_scans_to_consider(std::floor((min_sample_rate_ /2)*10));
//         Size min_scans_to_consider(5);

//         // double outlier_ratio(0.3);

//         // double ftl_mean(centroid_mz);
//         double ftl_sd((centroid_mz / 1e6) * mass_error_ppm_);
//         double intensity_so_far(apex_peak.getIntensity());

//         while (((trace_down_idx > 0) && toggle_down) ||
//                ((trace_up_idx < work_exp.size() - 1) && toggle_up)
//                 )
//         {
//           // *********************************************************** //
//           // Step 2.1 MOVE DOWN in RT dim
//           // *********************************************************** //
//           if ((trace_down_idx > 0) && toggle_down)
//           {
//             const MSSpectrum& spec_trace_down = work_exp[trace_down_idx - 1];
//             if (!spec_trace_down.empty())
//             {
//               Size next_down_peak_idx = spec_trace_down.findNearest(centroid_mz);
//               double next_down_peak_mz = spec_trace_down[next_down_peak_idx].getMZ();
//               double next_down_peak_int = spec_trace_down[next_down_peak_idx].getIntensity();

//               double right_bound = centroid_mz + 3 * ftl_sd;
//               double left_bound = centroid_mz - 3 * ftl_sd;

//               if ((next_down_peak_mz <= right_bound) &&
//                   (next_down_peak_mz >= left_bound) &&
//                   !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx]
//                       )
//               {
//                 Peak2D next_peak;
//                 next_peak.setRT(spec_trace_down.getRT());
//                 next_peak.setMZ(next_down_peak_mz);
//                 next_peak.setIntensity(next_down_peak_int);

//                 current_trace.push_front(next_peak);
//                 // FWHM average
//                 if (fwhm_meta_idx != -1)
//                 {
//                   fwhms_mz.push_back(spec_trace_down.getFloatDataArrays()[fwhm_meta_idx][next_down_peak_idx]);
//                 }
//                 // Update the m/z mean of the current trace as we added a new peak
//                 updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
//                 gathered_idx.emplace_back(trace_down_idx - 1, next_down_peak_idx);

//                 // Update the m/z variance dynamically
//                 if (reestimate_mt_sd_)           //  && (down_hitting_peak+1 > min_flank_scans))
//                 {
//                   // if (ftl_t > min_fwhm_scans)
//                   {
//                     updateWeightedSDEstimateRobust(next_peak, centroid_mz, ftl_sd, intensity_so_far);
//                   }
//                 }

//                 ++down_hitting_peak;
//                 conseq_missed_peak_down = 0;
//               }
//               else
//               {
//                 ++conseq_missed_peak_down;
//               }

//             }
//             --trace_down_idx;
//             ++down_scan_counter;

//             // trace termination criterion: max allowed number of
//             // consecutive outliers reached OR cancel extension if
//             // sampling_rate falls below min_sample_rate_
//             if (trace_termination_criterion_ == "outlier")
//             {
//               if (conseq_missed_peak_down > max_consecutive_missing)
//               {
//                 toggle_down = false;
//               }
//             }
//             else if (trace_termination_criterion_ == "sample_rate")
//             {
//               current_sample_rate = (double)(down_hitting_peak + up_hitting_peak + 1) /
//                                     (double)(down_scan_counter + up_scan_counter + 1);
//               if (down_scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
//               {
//                 // std::cout << "stopping down..." << std::endl;
//                 toggle_down = false;
//               }
//             }
//           }

//           // *********************************************************** //
//           // Step 2.2 MOVE UP in RT dim
//           // *********************************************************** //
//           if ((trace_up_idx < work_exp.size() - 1) && toggle_up)
//           {
//             const MSSpectrum& spec_trace_up = work_exp[trace_up_idx + 1];
//             if (!spec_trace_up.empty())
//             {
//               Size next_up_peak_idx = spec_trace_up.findNearest(centroid_mz);
//               double next_up_peak_mz = spec_trace_up[next_up_peak_idx].getMZ();
//               double next_up_peak_int = spec_trace_up[next_up_peak_idx].getIntensity();

//               double right_bound = centroid_mz + 3 * ftl_sd;
//               double left_bound = centroid_mz - 3 * ftl_sd;

//               if ((next_up_peak_mz <= right_bound) &&
//                   (next_up_peak_mz >= left_bound) &&
//                   !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
//               {
//                 Peak2D next_peak;
//                 next_peak.setRT(spec_trace_up.getRT());
//                 next_peak.setMZ(next_up_peak_mz);
//                 next_peak.setIntensity(next_up_peak_int);

//                 current_trace.push_back(next_peak);
//                 if (fwhm_meta_idx != -1)
//                 {
//                   fwhms_mz.push_back(spec_trace_up.getFloatDataArrays()[fwhm_meta_idx][next_up_peak_idx]);
//                 }
//                 // Update the m/z mean of the current trace as we added a new peak
//                 updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
//                 gathered_idx.emplace_back(trace_up_idx + 1, next_up_peak_idx);

//                 // Update the m/z variance dynamically
//                 if (reestimate_mt_sd_)           //  && (up_hitting_peak+1 > min_flank_scans))
//                 {
//                   // if (ftl_t > min_fwhm_scans)
//                   {
//                     updateWeightedSDEstimateRobust(next_peak, centroid_mz, ftl_sd, intensity_so_far);
//                   }
//                 }

//                 ++up_hitting_peak;
//                 conseq_missed_peak_up = 0;

//               }
//               else
//               {
//                 ++conseq_missed_peak_up;
//               }

//             }

//             ++trace_up_idx;
//             ++up_scan_counter;

//             if (trace_termination_criterion_ == "outlier")
//             {
//               if (conseq_missed_peak_up > max_consecutive_missing)
//               {
//                 toggle_up = false;
//               }
//             }
//             else if (trace_termination_criterion_ == "sample_rate")
//             {
//               current_sample_rate = (double)(down_hitting_peak + up_hitting_peak + 1) / (double)(down_scan_counter + up_scan_counter + 1);

//               if (up_scan_counter > min_scans_to_consider && current_sample_rate < min_sample_rate_)
//               {
//                 // std::cout << "stopping up" << std::endl;
//                 toggle_up = false;
//               }
//             }


//           }

//         }

//         // std::cout << "current sr: " << current_sample_rate << std::endl;
//         double num_scans(down_scan_counter + up_scan_counter + 1 - conseq_missed_peak_down - conseq_missed_peak_up);

//         double mt_quality((double)current_trace.size() / (double)num_scans);
//         // std::cout << "mt quality: " << mt_quality << std::endl;
//         double rt_range(std::fabs(current_trace.rbegin()->getRT() - current_trace.begin()->getRT()));

//         // *********************************************************** //
//         // Step 2.3 check if minimum length and quality of mass trace criteria are met
//         // *********************************************************** //
//         bool max_trace_criteria = (max_trace_length_ < 0.0 || rt_range < max_trace_length_);
//         if (rt_range >= min_trace_length_ && max_trace_criteria && mt_quality >= min_sample_rate_)
//         {
//           // std::cout << "T" << trace_number << "\t" << mt_quality << std::endl;

//           // mark all peaks as visited
//           for (Size i = 0; i < gathered_idx.size(); ++i)
//           {
//             peak_visited[spec_offsets[gathered_idx[i].first] +  gathered_idx[i].second] = true;
//           }

//           // create new MassTrace object and store collected peaks from list current_trace
//           MassTrace new_trace(current_trace);
//           new_trace.updateWeightedMeanRT();
//           new_trace.updateWeightedMeanMZ();
//           if (!fwhms_mz.empty())
//           {
//             new_trace.fwhm_mz_avg = Math::median(fwhms_mz.begin(), fwhms_mz.end());
//           }
//           new_trace.setQuantMethod(quant_method_);
//           //new_trace.setCentroidSD(ftl_sd);
//           new_trace.updateWeightedMZsd();
//           new_trace.setLabel("T" + String(trace_number));
//           ++trace_number;

//           found_masstraces.push_back(new_trace);

//           peaks_detected += new_trace.getSize();
//           this->setProgress(peaks_detected);

//           // check if we already reached the (optional) maximum number of traces
//           if (max_traces > 0 && found_masstraces.size() == max_traces)
//           {
//             break;
//           }
//         }
//       }

//       this->endProgress();

//     }

//     void MassTraceDetection::updateMembers_()
//     {
//       mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
//       noise_threshold_int_ = (double)param_.getValue("noise_threshold_int");
//       chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
//       quant_method_ = MassTrace::getQuantMethod((String)param_.getValue("quant_method").toString());

//       trace_termination_criterion_ = (String)param_.getValue("trace_termination_criterion").toString();
//       trace_termination_outliers_ = (Size)param_.getValue("trace_termination_outliers");
//       min_sample_rate_ = (double)param_.getValue("min_sample_rate");
//       min_trace_length_ = (double)param_.getValue("min_trace_length");
//       max_trace_length_ = (double)param_.getValue("max_trace_length");
//       reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
//     }
// }