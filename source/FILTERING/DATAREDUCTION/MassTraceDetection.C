// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
    // defaults_.setValue( "name" , 1 , "descript" );
    defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass deviation (in ppm).");
    defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are removed as noise.");
    defaults_.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");

    defaults_.setValue("reestimate_mt_sd", "true", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
    defaults_.setValidStrings("reestimate_mt_sd", StringList::create(("true,false")));

    // advanced parameters
    defaults_.setValue("min_sample_rate", 0.5, "Minimum fraction of scans along the mass trace that must contain a peak.", StringList::create("advanced"));
    defaults_.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).", StringList::create("advanced"));
    defaults_.setValue("max_trace_length", 300.0, "Minimum expected length of a mass trace (in seconds).", StringList::create("advanced"));




    defaultsToParam_();

    this->setLogType(CMD);
}

MassTraceDetection::~MassTraceDetection()
{

}

void MassTraceDetection::updateIterativeWeightedMeanMZ(const DoubleReal & added_mz, const DoubleReal & added_int, DoubleReal & centroid_mz, DoubleReal & prev_counter, DoubleReal & prev_denom)
{
    DoubleReal new_weight(added_int);
    DoubleReal new_mz(added_mz);

    DoubleReal counter_tmp(1 + (new_weight * new_mz) / prev_counter);
    DoubleReal denom_tmp(1 + (new_weight) / prev_denom);
    centroid_mz *= (counter_tmp / denom_tmp);
    prev_counter *= counter_tmp;
    prev_denom *= denom_tmp;

    return;
}

void MassTraceDetection::run(MSExperiment<Peak1D>::ConstAreaIterator & begin, MSExperiment<Peak1D>::ConstAreaIterator & end, std::vector<MassTrace> & found_masstraces)
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

void updateMeanEstimate(const DoubleReal & x_t, DoubleReal & mean_t, Size t)
{
    DoubleReal tmp(mean_t);

    tmp = mean_t + (1 / (t + 1)) * (x_t - mean_t);

    mean_t = tmp;
}

void updateSDEstimate(const DoubleReal & x_t, const DoubleReal & mean_t, DoubleReal & sd_t, Size t)
{
    DoubleReal tmp(sd_t);
    DoubleReal i(t);


    tmp = (i / (i + 1)) * sd_t + (i / (i + 1) * (i + 1)) * (x_t - mean_t) * (x_t - mean_t);

    sd_t = tmp;
    // std::cerr << "func:  " << tmp << " " << i << std::endl;
}

void updateWeightedSDEstimate(PeakType p, const DoubleReal & mean_t1, DoubleReal & sd_t, DoubleReal & last_weights_sum)
{
    DoubleReal denom(0.0), weights_sum(0.0);

    denom = last_weights_sum * sd_t * sd_t + p.getIntensity() * (p.getMZ() - mean_t1) * (p.getMZ() - mean_t1);
    weights_sum = last_weights_sum + p.getIntensity();

    DoubleReal tmp_sd = std::sqrt(denom / weights_sum);

    if (tmp_sd > std::numeric_limits<DoubleReal>::epsilon())
    {
        sd_t = tmp_sd;
    }

    last_weights_sum = weights_sum;
}

void updateWeightedSDEstimateRobust(PeakType p, const DoubleReal & mean_t1, DoubleReal & sd_t, DoubleReal & last_weights_sum)
{
    DoubleReal denom(0.0), denom1(0.0), denom2(0.0), weights_sum(0.0);

    denom1 = std::log(last_weights_sum) + 2 * std::log(sd_t);
    denom2 = std::log(p.getIntensity()) + 2 * std::log(std::abs(p.getMZ() - mean_t1));

    denom = std::sqrt(std::exp(denom1) + std::exp(denom2));
    weights_sum = last_weights_sum + p.getIntensity();

    DoubleReal tmp_sd = denom / std::sqrt(weights_sum);

    if (tmp_sd > std::numeric_limits<DoubleReal>::epsilon())
    {
        sd_t = tmp_sd;
    }

    last_weights_sum = weights_sum;
}

void computeWeightedSDEstimate(std::list<PeakType> tmp, const DoubleReal & mean_t, DoubleReal & sd_t, const DoubleReal & lower_sd_bound)
{
    DoubleReal denom(0.0), weights_sum(0.0);

    for (std::list<PeakType>::const_iterator l_it = tmp.begin(); l_it != tmp.end(); ++l_it)
    {
        denom += l_it->getIntensity() * (l_it->getMZ() - mean_t) * (l_it->getMZ() - mean_t);
        weights_sum += l_it->getIntensity();
    }

    DoubleReal tmp_sd = std::sqrt(denom / (weights_sum));

    // std::cout << "tmp_sd" << tmp_sd << std::endl;

    if (tmp_sd > std::numeric_limits<DoubleReal>::epsilon())
    {
        sd_t = tmp_sd;
    }

    return;
}

DoubleReal computeLoss(const DoubleReal & x_t, const DoubleReal & mean_t, const DoubleReal & sd_t)
{
    return ((x_t - mean_t) * (x_t - mean_t)) / (2 * sd_t * sd_t) + 0.5 * std::log(sd_t * sd_t);
}

void MassTraceDetection::run(const MSExperiment<Peak1D> & input_exp, std::vector<MassTrace> & found_masstraces)
{
    // make sure the output vector is empty
    found_masstraces.clear();

    // gather all peaks that are potential chromatographic peak apeces
    typedef std::multimap<DoubleReal, std::pair<Size, Size> > MapIdxSortedByInt;
    MSExperiment<Peak1D> work_exp;
    MapIdxSortedByInt chrom_apeces;

    Size peak_count(0);
    std::vector<Size> spec_offsets;
    spec_offsets.push_back(0);

    Size spectra_count(0);

    // this->startProgress(0, input_exp.size(), "Detect potential chromatographic apeces...");
    for (Size scan_idx = 0; scan_idx < input_exp.size(); ++scan_idx)
    {
        // this->setProgress(scan_idx);

        // check if this is a MS1 survey scan
        if (input_exp[scan_idx].getMSLevel() == 1)
        {
            DoubleReal scan_rt = input_exp[scan_idx].getRT();
            MSSpectrum<Peak1D> tmp_spec;
            Size spec_peak_idx = 0;

            tmp_spec.setRT(scan_rt);

            for (Size peak_idx = 0; peak_idx < input_exp[scan_idx].size(); ++peak_idx)
            {
                DoubleReal tmp_peak_int(input_exp[scan_idx][peak_idx].getIntensity());

                if (tmp_peak_int > noise_threshold_int_)
                {
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
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input map consists of too few spectra (less than 3!). Aborting...", String(spectra_count));
    }


    Size min_fwhm_scans(3);
    DoubleReal scan_time(std::fabs(input_exp[input_exp.size() - 1].getRT() - input_exp[0].getRT()) / input_exp.size());

    // std::cout << "scan_time: " << scan_time << " " << min_trace_length_ << std::endl;

    Size scan_nums(std::floor(min_trace_length_ / scan_time));

    if (scan_nums > min_fwhm_scans)
    {
        min_fwhm_scans = scan_nums;
    }


    Size min_flank_scans(3);

    // discard last spectrum's offset
    spec_offsets.pop_back();

    boost::dynamic_bitset<> peak_visited(peak_count);

    // start extending mass traces beginning with the apex peak

    Size trace_number(1);

    this->startProgress(0, peak_count, "mass trace detection");
    Size peaks_detected(0);

    for (MapIdxSortedByInt::reverse_iterator m_it = chrom_apeces.rbegin(); m_it != chrom_apeces.rend(); ++m_it)
    {
        Size apex_scan_idx(m_it->second.first);
        Size apex_peak_idx(m_it->second.second);

        if (peak_visited[spec_offsets[apex_scan_idx] + apex_peak_idx])
            continue;

        Peak2D apex_peak;
        apex_peak.setRT(work_exp[apex_scan_idx].getRT());
        apex_peak.setMZ(work_exp[apex_scan_idx][apex_peak_idx].getMZ());
        apex_peak.setIntensity(work_exp[apex_scan_idx][apex_peak_idx].getIntensity());

        Size trace_up_idx(apex_scan_idx);
        Size trace_down_idx(apex_scan_idx);

        std::list<PeakType> current_trace;
        current_trace.push_back(apex_peak);

        // Initialization for the iterative version of weighted m/z mean calculation
        DoubleReal centroid_mz(apex_peak.getMZ());
        DoubleReal prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
        DoubleReal prev_denom(apex_peak.getIntensity());

        updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

        std::vector<std::pair<Size, Size> > gathered_idx;
        gathered_idx.push_back(std::make_pair(apex_scan_idx, apex_peak_idx));

        Size up_hitting_peak(0), down_hitting_peak(0);
        Size up_scan_counter(0), down_scan_counter(0);

        bool toggle_up = true, toggle_down = true;

        Size conseq_missed_peak_up(0), conseq_missed_peak_down(0);
        Size MAX_CONSEQ_MISSING(5);

        // DoubleReal outlier_ratio(0.3);

        DoubleReal ftl_mean(centroid_mz), ftl_sd((centroid_mz / 1000000) * mass_error_ppm_);
        DoubleReal intensity_so_far(apex_peak.getIntensity());

        while (((trace_down_idx > 0) && toggle_down) || ((trace_up_idx < work_exp.size() - 1) && toggle_up))
        {

            // DoubleReal centroid_mz = current_trace.getCentroidMZ();

            // try to go downwards in RT
            if (((trace_down_idx > 0) && toggle_down))
            {
                try
                {
                    Size next_down_peak_idx = work_exp[trace_down_idx - 1].findNearest(centroid_mz);
                    DoubleReal next_down_peak_mz = work_exp[trace_down_idx - 1][next_down_peak_idx].getMZ();
                    DoubleReal next_down_peak_int = work_exp[trace_down_idx - 1][next_down_peak_idx].getIntensity();

                    DoubleReal right_bound, left_bound;

                    //                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                    //                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;

                    right_bound = centroid_mz + 3 * ftl_sd;
                    left_bound = centroid_mz - 3 * ftl_sd;

                    // std::cout << "down: " << centroid_mz << " "<<  ftl_sd << std::endl;

                    Size left_next_idx = work_exp[trace_down_idx - 1].findNearest(left_bound);
                    Size right_next_idx = work_exp[trace_down_idx - 1].findNearest(right_bound);

                    DoubleReal left_mz(work_exp[trace_down_idx - 1][left_next_idx].getMZ());
                    DoubleReal right_mz(work_exp[trace_down_idx - 1][right_next_idx].getMZ());


                    if ((next_down_peak_mz <= right_bound) && (next_down_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx])
                    {
                        Peak2D next_peak;
                        next_peak.setRT(work_exp[trace_down_idx - 1].getRT());
                        next_peak.setMZ(next_down_peak_mz);
                        next_peak.setIntensity(next_down_peak_int);

                        current_trace.push_front(next_peak);

                        updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
                        gathered_idx.push_back(std::make_pair(trace_down_idx - 1, next_down_peak_idx));



                        if (reestimate_mt_sd_) //  && (down_hitting_peak+1 > min_flank_scans))
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
                catch (...)
                {
                    // std::cerr << "findNearest() ran into troubles..." << std::endl;
                }
                --trace_down_idx;
                ++down_scan_counter;


                if (conseq_missed_peak_down > MAX_CONSEQ_MISSING)
                {
                    toggle_down = false;
                }

                //                                if (down_scan_counter > min_flank_scans)
                //                                {
                //                                    //std::cout << "down ratio: " << (DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter << std::endl;
                //                                    DoubleReal trace_sampling_ratio_down((DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter);

                //                                    if (trace_sampling_ratio_down < min_sample_rate_)
                //                                    {
                //                                        toggle_down = false;
                //                                    }
                //                                }

            }

            //}
            // *********************************************************** //
            // MOVE UP in RT dim
            // *********************************************************** //

            if (((trace_up_idx < work_exp.size() - 1) && toggle_up))
            {
                try
                {
                    Size next_up_peak_idx = work_exp[trace_up_idx + 1].findNearest(centroid_mz);
                    DoubleReal next_up_peak_mz = work_exp[trace_up_idx + 1][next_up_peak_idx].getMZ();
                    DoubleReal next_up_peak_int = work_exp[trace_up_idx + 1][next_up_peak_idx].getIntensity();

                    DoubleReal right_bound, left_bound;

                    //                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                    //                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;


                    right_bound = centroid_mz + 3 * ftl_sd;
                    left_bound = centroid_mz - 3 * ftl_sd;


                    if ((next_up_peak_mz <= right_bound) && (next_up_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
                    {
                        Peak2D next_peak;
                        next_peak.setRT(work_exp[trace_up_idx + 1].getRT());
                        next_peak.setMZ(next_up_peak_mz);
                        next_peak.setIntensity(next_up_peak_int);

                        current_trace.push_back(next_peak);

                        updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
                        gathered_idx.push_back(std::make_pair(trace_up_idx + 1, next_up_peak_idx));


                        if (reestimate_mt_sd_) //  && (up_hitting_peak+1 > min_flank_scans))
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
                catch (...)
                {
                    //  std::cerr << "findNearest() ran into troubles..." << std::endl;
                }

                ++trace_up_idx;
                ++up_scan_counter;

                if (conseq_missed_peak_up > MAX_CONSEQ_MISSING)
                {
                    toggle_up = false;
                }

                //                                if (up_scan_counter > min_flank_scans)
                //                                {
                //                                    //std::cout << "down ratio: " << (DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter << std::endl;
                //                                    DoubleReal trace_sampling_ratio_up((DoubleReal)up_hitting_peak/(DoubleReal)up_scan_counter);

                //                                    if (trace_sampling_ratio_up < min_sample_rate_)
                //                                    {
                //                                        toggle_up = false;
                //                                    }
                //                                }


            }

        }

        DoubleReal num_scans(down_scan_counter + up_scan_counter + 1 - conseq_missed_peak_down - conseq_missed_peak_up);

        DoubleReal mt_quality((DoubleReal)current_trace.size() / (DoubleReal)num_scans);
        DoubleReal rt_range(std::fabs(current_trace.rbegin()->getRT() - current_trace.begin()->getRT()));

        // check if minimum length and quality of mass trace criteria are met
        if (rt_range >= min_trace_length_ && rt_range < max_trace_length_ && mt_quality >= min_sample_rate_)
        {
            // std::cout << "T" << trace_number << "\t" << mt_quality << std::endl;

            // mark all peaks as visited
            for (Size i = 0; i < gathered_idx.size(); ++i)
            {
                peak_visited[spec_offsets[gathered_idx[i].first] +  gathered_idx[i].second] = true;
            }

            String tr_num;
            std::stringstream read_in;
            read_in << trace_number;
            tr_num = read_in.str();

            // create new MassTrace object and store collected peaks from list current_trace
            MassTrace new_trace(current_trace, scan_time);
            new_trace.updateWeightedMeanRT();
            new_trace.updateWeightedMeanMZ();

            new_trace.setCentroidSD(ftl_sd);

            new_trace.setLabel("T" + tr_num);

            peaks_detected += new_trace.getSize();
            this->setProgress(peaks_detected);
            found_masstraces.push_back(new_trace);
            ++trace_number;
        }
    }

    this->endProgress();

    return;
} // end of MassTraceDetection::run

void MassTraceDetection::updateMembers_()
{
    mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
    noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
    chrom_peak_snr_ = (DoubleReal)param_.getValue("chrom_peak_snr");
    // chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");

    min_sample_rate_ = (DoubleReal)param_.getValue("min_sample_rate");
    min_trace_length_ = (DoubleReal)param_.getValue("min_trace_length");
    max_trace_length_ = (DoubleReal)param_.getValue("max_trace_length");
    reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
}

}
