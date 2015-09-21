// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>

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
    defaults_.setValidStrings("reestimate_mt_sd", ListUtils::create<String>("true,false"));

    // advanced parameters
    defaults_.setValue("trace_termination_criterion", "outlier", "Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("trace_termination_criterion", ListUtils::create<String>("outlier,sample_rate"));
    defaults_.setValue("trace_termination_outliers", 5, "Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.", ListUtils::create<String>("advanced"));

    defaults_.setValue("min_sample_rate", 0.5, "Minimum fraction of scans along the mass trace that must contain a peak.", ListUtils::create<String>("advanced"));
    defaults_.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_trace_length", 300.0, "Maximum expected length of a mass trace (in seconds).", ListUtils::create<String>("advanced"));

    defaultsToParam_();

    this->setLogType(CMD);
  }

  MassTraceDetection::~MassTraceDetection()
  {
  }

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

  void MassTraceDetection::run(MSExperiment<Peak1D>::ConstAreaIterator& begin,
                               MSExperiment<Peak1D>::ConstAreaIterator& end, std::vector<MassTrace>&
                               found_masstraces)
  {
    MSExperiment<Peak1D> map;
    MSSpectrum<Peak1D> current_spectrum;

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

  void updateWeightedSDEstimate(PeakType p, const double& mean_t1, double& sd_t, double& last_weights_sum)
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

  void updateWeightedSDEstimateRobust(PeakType p, const double& mean_t1, double& sd_t, double& last_weights_sum)
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

  double computeLoss(const double& x_t, const double& mean_t, const double& sd_t)
  {
    return ((x_t - mean_t) * (x_t - mean_t)) / (2 * sd_t * sd_t) + 0.5 * std::log(sd_t * sd_t);
  }

  void MassTraceDetection::run(const MSExperiment<Peak1D>& input_exp, std::vector<MassTrace>& found_masstraces)
  {
    // make sure the output vector is empty
    found_masstraces.clear();

    // gather all peaks that are potential chromatographic peak apeces
    //   - use work_exp for actual work (remove peaks below noise threshold)
    //   - store potential apices in chrom_apeces
    MSExperiment<Peak1D> work_exp;
    MapIdxSortedByInt chrom_apeces;

    Size peak_count(0);
    std::vector<Size> spec_offsets;
    spec_offsets.push_back(0);

    Size spectra_count(0);

    // *********************************************************** //
    //  Step 1: Detecting potential chromatographic apices
    // *********************************************************** //
    for (Size scan_idx = 0; scan_idx < input_exp.size(); ++scan_idx)
    {
      // check if this is a MS1 survey scan
      if (input_exp[scan_idx].getMSLevel() == 1)
      {
        double scan_rt = input_exp[scan_idx].getRT();
        MSSpectrum<Peak1D> tmp_spec;
        Size spec_peak_idx = 0;
        tmp_spec.setRT(scan_rt);

        for (Size peak_idx = 0; peak_idx < input_exp[scan_idx].size(); ++peak_idx)
        {
          double tmp_peak_int(input_exp[scan_idx][peak_idx].getIntensity());
          if (tmp_peak_int > noise_threshold_int_)
          {
            // Assume that noise_threshold_int_ contains the noise level of the
            // data and we want to be chrom_peak_snr times above the noise
            // level
            // -> add this peak as possible chromatographic apex
            tmp_spec.push_back(input_exp[scan_idx][peak_idx]);
            if (tmp_peak_int > chrom_peak_snr_ * noise_threshold_int_)
            {
              chrom_apeces.insert(std::make_pair(tmp_peak_int, std::make_pair(scan_idx, spec_peak_idx)));
            }
            ++peak_count;
            ++spec_peak_idx;
          }
        }
        work_exp.addSpectrum(tmp_spec);
        spec_offsets.push_back(spec_offsets[spec_offsets.size() - 1] + tmp_spec.size());
        ++spectra_count;
      }
    }

    if (spectra_count < 3)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                    "Input map consists of too few spectra (less than 3!). Aborting...", String(spectra_count));
    }

    // discard last spectrum's offset
    spec_offsets.pop_back();

    // *********************************************************************
    // Step 2: start extending mass traces beginning with the apex peak (go
    // through all peaks in order of decreasing intensity)
    // *********************************************************************
    run_(chrom_apeces, peak_count, work_exp, spec_offsets, found_masstraces);

    return;
  } // end of MassTraceDetection::run

  void MassTraceDetection::run_(const MapIdxSortedByInt& chrom_apeces, Size peak_count, 
                                const MSExperiment<Peak1D> & work_exp, 
                                const std::vector<Size>& spec_offsets,
                                std::vector<MassTrace> & found_masstraces)
  {
    // Size min_flank_scans(3);
    boost::dynamic_bitset<> peak_visited(peak_count);
    Size trace_number(1);

    this->startProgress(0, peak_count, "mass trace detection");
    Size peaks_detected(0);

    for (MapIdxSortedByInt::const_reverse_iterator m_it = chrom_apeces.rbegin(); m_it != chrom_apeces.rend(); ++m_it)
    {
      Size apex_scan_idx(m_it->second.first);
      Size apex_peak_idx(m_it->second.second);

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

      // Initialization for the iterative version of weighted m/z mean calculation
      double centroid_mz(apex_peak.getMZ());
      double prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
      double prev_denom(apex_peak.getIntensity());

      updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

      std::vector<std::pair<Size, Size> > gathered_idx;
      gathered_idx.push_back(std::make_pair(apex_scan_idx, apex_peak_idx));

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

        // double centroid_mz = current_trace.getCentroidMZ();

        // *********************************************************** //
        // Step 2.1 MOVE DOWN in RT dim
        // *********************************************************** //
        if ((trace_down_idx > 0) && toggle_down)
        {
          if (!work_exp[trace_down_idx - 1].empty())
          {
            Size next_down_peak_idx = work_exp[trace_down_idx - 1].findNearest(centroid_mz);
            double next_down_peak_mz = work_exp[trace_down_idx - 1][next_down_peak_idx].getMZ();
            double next_down_peak_int = work_exp[trace_down_idx - 1][next_down_peak_idx].getIntensity();

            double right_bound, left_bound;

            //                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
            //                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;

            right_bound = centroid_mz + 3 * ftl_sd;
            left_bound = centroid_mz - 3 * ftl_sd;

            //                  std::cout << "down: " << centroid_mz << " "<<  ftl_sd << std::endl;

            // Size left_next_idx = work_exp[trace_down_idx - 1].findNearest(left_bound);
            // Size right_next_idx = work_exp[trace_down_idx - 1].findNearest(right_bound);

            // double left_mz(work_exp[trace_down_idx - 1][left_next_idx].getMZ());
            // double right_mz(work_exp[trace_down_idx - 1][right_next_idx].getMZ());

            if ((next_down_peak_mz <= right_bound) &&
                (next_down_peak_mz >= left_bound) &&
                !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx]
                )
            {
              Peak2D next_peak;
              next_peak.setRT(work_exp[trace_down_idx - 1].getRT());
              next_peak.setMZ(next_down_peak_mz);
              next_peak.setIntensity(next_down_peak_int);

              current_trace.push_front(next_peak);

              // Update the m/z mean of the current trace as we added a new peak
              updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
              gathered_idx.push_back(std::make_pair(trace_down_idx - 1, next_down_peak_idx));

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
          if (!work_exp[trace_up_idx + 1].empty())
          {
            Size next_up_peak_idx = work_exp[trace_up_idx + 1].findNearest(centroid_mz);
            double next_up_peak_mz = work_exp[trace_up_idx + 1][next_up_peak_idx].getMZ();
            double next_up_peak_int = work_exp[trace_up_idx + 1][next_up_peak_idx].getIntensity();

            double right_bound, left_bound;

            //                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
            //                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;


            right_bound = centroid_mz + 3 * ftl_sd;
            left_bound = centroid_mz - 3 * ftl_sd;


            if ((next_up_peak_mz <= right_bound) &&
                (next_up_peak_mz >= left_bound) &&
                !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
            {
              Peak2D next_peak;
              next_peak.setRT(work_exp[trace_up_idx + 1].getRT());
              next_peak.setMZ(next_up_peak_mz);
              next_peak.setIntensity(next_up_peak_int);

              current_trace.push_back(next_peak);

              // Update the m/z mean of the current trace as we added a new peak
              updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
              gathered_idx.push_back(std::make_pair(trace_up_idx + 1, next_up_peak_idx));

              // Update the m/z variance dynamically
              if (reestimate_mt_sd_)           //  && (up_hitting_peak+1 > min_flank_scans))
              {
                // if (ftl_t > min_fwhm_scans)
                {
                  //computeWeightedSDEstimate(current_trace, centroid_mz, ftl_sd, lower_sd_bound);
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
      if (rt_range >= min_trace_length_ && rt_range < max_trace_length_ && mt_quality >= min_sample_rate_)
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

        //new_trace.setCentroidSD(ftl_sd);
        new_trace.updateWeightedMZsd();

        new_trace.setLabel("T" + String(trace_number));

        peaks_detected += new_trace.getSize();
        this->setProgress(peaks_detected);
        found_masstraces.push_back(new_trace);
        ++trace_number;
      }
    }

    this->endProgress();

  }

  void MassTraceDetection::updateMembers_()
  {
    mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
    noise_threshold_int_ = (double)param_.getValue("noise_threshold_int");
    chrom_peak_snr_ = (double)param_.getValue("chrom_peak_snr");
    // chrom_fwhm_ = (double)param_.getValue("chrom_fwhm");

    trace_termination_criterion_ = (String)param_.getValue("trace_termination_criterion");
    trace_termination_outliers_ = (Size)param_.getValue("trace_termination_outliers");
    min_sample_rate_ = (double)param_.getValue("min_sample_rate");
    min_trace_length_ = (double)param_.getValue("min_trace_length");
    max_trace_length_ = (double)param_.getValue("max_trace_length");
    reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
  }

}
