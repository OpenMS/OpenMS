// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
MassTraceDetection::MassTraceDetection()
    : DefaultParamHandler("MassTraceDetection"), ProgressLogger()
{
    // defaults_.setValue( "name" , 1 , "descript" );
    defaults_.setValue("mass_error_ppm" , 20.0 , "Allowed mass deviation (in ppm).");
    defaults_.setValue("noise_threshold_int" , 10.0 , "Intensity threshold below which peaks are removed as noise.");
    defaults_.setValue("chrom_peak_snr" , 3.0 , "Minimum signal-to-noise a mass trace should have.");

    defaults_.setValue("reestimate_mt_sd", "true", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
    defaults_.setValidStrings("reestimate_mt_sd", StringList::create(("true,false")));

    // advanced parameters
    defaults_.setValue("min_sample_rate" , 0.7 , "Minimum fraction of scans along the mass trace that must contain a peak.", StringList::create("advanced"));
    defaults_.setValue("min_peak_width" , 3.0 , "Minimum expected length of a mass trace (in seconds).", StringList::create("advanced"));



    defaultsToParam_();

    this->setLogType(CMD);
}

MassTraceDetection::~MassTraceDetection()
{

}


void MassTraceDetection::updateIterativeWeightedMeanMZ(const DoubleReal& added_mz, const DoubleReal& added_int, DoubleReal& centroid_mz, DoubleReal& prev_counter, DoubleReal& prev_denom)
{
    DoubleReal new_weight(added_int);
    DoubleReal new_mz(added_mz);

    DoubleReal counter_tmp(1 + (new_weight*new_mz)/prev_counter);
    DoubleReal denom_tmp(1 + (new_weight)/prev_denom);
    centroid_mz *= (counter_tmp/denom_tmp);
    prev_counter *= counter_tmp;
    prev_denom *= denom_tmp;

    return ;
}


void MassTraceDetection::filterByPeakWidth(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& filt_mtraces)
{
    std::multimap<Size, Size> histo_map;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        mt_vec[i].estimateFWHM(false);
        Size fwhm(mt_vec[i].getFWHMScansNum());

        histo_map.insert(std::make_pair(fwhm, i));
    }

    // compute median peak width
    std::vector<DoubleReal> pw_vec;
    std::vector<Size> pw_idx_vec;

    for (std::multimap<Size, Size>::const_iterator c_it = histo_map.begin(); c_it != histo_map.end(); ++c_it)
    {
        pw_vec.push_back(c_it->first);
        pw_idx_vec.push_back(c_it->second);
    }

    // Size pw_vec_size = pw_vec.size();
    DoubleReal pw_median(Math::median(pw_vec.begin(), pw_vec.end(), true));


    // compute median of absolute deviances (MAD)
    std::vector<DoubleReal> abs_devs;

    for (Size pw_i = 0; pw_i < pw_vec.size(); ++pw_i)
    {
        abs_devs.push_back(std::fabs(pw_vec[pw_i] - pw_median));
    }

    // Size abs_devs_size = abs_devs.size();
    DoubleReal pw_mad(Math::median(abs_devs.begin(), abs_devs.end(), false));

    DoubleReal lower_pw_bound(0.0);

    if (pw_median - 2*pw_mad > 0.0)
    {
        lower_pw_bound = pw_median - 2*pw_mad;
    }

    // DoubleReal upper_pw_bound(std::floor(pw_median + 2*pw_mad));

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        // set to lowest peak width according to distribution
        if (pw_vec[i] < lower_pw_bound)
        {
            if (mt_vec[pw_idx_vec[i]].getSize() >= lower_pw_bound)
            {
                if (mt_vec[pw_idx_vec[i]].getSize() >= pw_median)
                {
                    mt_vec[pw_idx_vec[i]].setFWHMScansNum((Size)pw_median);
                }
                else
                {
                    mt_vec[pw_idx_vec[i]].setFWHMScansNum((Size)lower_pw_bound); // override "false" pw estimation
                }

                filt_mtraces.push_back(mt_vec[pw_idx_vec[i]]);
            }
        }
        else
        {
            filt_mtraces.push_back(mt_vec[pw_idx_vec[i]]);
        }
    }

    return ;
}

void MassTraceDetection::run(MSExperiment<Peak1D>::ConstAreaIterator& begin, MSExperiment<Peak1D>::ConstAreaIterator& end, std::vector<MassTrace>& found_masstraces)
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
                map.push_back(current_spectrum);
            }
            current_spectrum.clear(false);
            current_spectrum.setRT(begin.getRT());
        }
        current_spectrum.push_back(*begin);
    }
    map.push_back(current_spectrum);

    run(map, found_masstraces);
}


// update function for FTL method

void updateMeanEstimate(const DoubleReal& x_t, DoubleReal& mean_t, Size t)
{
    DoubleReal tmp(mean_t);

    tmp = mean_t + (1/(t + 1))*(x_t - mean_t);

    mean_t = tmp;
}

void updateSDEstimate(const DoubleReal& x_t, const DoubleReal& mean_t, DoubleReal& sd_t, Size t)
{
    DoubleReal tmp(sd_t);
    DoubleReal i(t);


    tmp = (i/(i+1))*sd_t + (i/(i+1)*(i+1))*(x_t - mean_t)*(x_t - mean_t);

    sd_t = tmp;
    // std::cerr << "func:  " << tmp << " " << i << std::endl;
}

void  updateWeightedSDEstimate(PeakType p, const DoubleReal& mean_t1, DoubleReal& sd_t, DoubleReal& last_weights_sum)
{
    DoubleReal denom(0.0), weights_sum(0.0);

    denom = last_weights_sum*sd_t*sd_t + p.getIntensity()*(p.getMZ() - mean_t1)*(p.getMZ() - mean_t1);
    weights_sum = last_weights_sum + p.getIntensity();

    DoubleReal tmp_sd = std::sqrt(denom/weights_sum);

    if (tmp_sd > std::numeric_limits<DoubleReal>::epsilon())
    {
        sd_t = tmp_sd;
    }

    last_weights_sum = weights_sum;
}


void computeWeightedSDEstimate(std::list<PeakType> tmp, const DoubleReal& mean_t, DoubleReal& sd_t, const DoubleReal& lower_sd_bound)
{
    DoubleReal denom(0.0), weights_sum(0.0);

    for (std::list<PeakType>::const_iterator l_it = tmp.begin(); l_it != tmp.end(); ++l_it)
    {
        denom += l_it->getIntensity()*(l_it->getMZ() - mean_t)*(l_it->getMZ() - mean_t);
        weights_sum += l_it->getIntensity();
    }

    DoubleReal tmp_sd = std::sqrt(denom/(weights_sum));

    std::cout << "tmp_sd" << tmp_sd << std::endl;

    if (tmp_sd > std::numeric_limits<DoubleReal>::epsilon())
    {
        sd_t = tmp_sd;
    }

    return ;
}


DoubleReal computeLoss(const DoubleReal& x_t, const DoubleReal& mean_t, const DoubleReal& sd_t)
{
    return ((x_t - mean_t)*(x_t - mean_t))/(2*sd_t*sd_t) + 0.5*std::log(sd_t*sd_t);
}


void MassTraceDetection::run(const MSExperiment<Peak1D>& input_exp, std::vector<MassTrace>& found_masstraces)
{
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

                    if (tmp_peak_int > chrom_peak_snr_*noise_threshold_int_)
                    {
                        chrom_apeces.insert(std::make_pair(tmp_peak_int, std::make_pair(scan_idx, spec_peak_idx)));
                    }
                    ++peak_count;
                    ++spec_peak_idx;
                }
            }

            work_exp.push_back(tmp_spec);
            spec_offsets.push_back(spec_offsets[spec_offsets.size() - 1] + tmp_spec.size());

            ++spectra_count;
        }
    }

    if (spectra_count < 3)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input map consists of too few spectra (less than 3!). Aborting...", String(spectra_count));
    }


    DoubleReal scan_time(std::fabs(input_exp[input_exp.size() - 1].getRT() - input_exp[0].getRT())/input_exp.size());
    Size min_fwhm_scans(3);

    Size scan_nums(std::floor(min_peak_width_/scan_time));

    if (scan_nums > min_fwhm_scans)
    {
        min_fwhm_scans = scan_nums;
    }


    // std::cout << "min_fwhm_scans: " << min_fwhm_scans << " fwhm: " << chrom_fwhm_ << " time: " << scan_time << std::endl;

    // this->endProgress();

    // discard last spectrum's offset
    spec_offsets.pop_back();

    boost::dynamic_bitset<> peak_visited(peak_count);

    // start extending mass traces beginning with the apex peak

    // Size min_datapoints = std::floor(chrom_fwhm_/(scan_rt_diff*2));
    Size min_data_points(3);

    // Size min_trace_quality = std::floor(2*min_datapoints*min_sample_rate_);

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

        DoubleReal half_max_int(work_exp[apex_scan_idx][apex_peak_idx].getIntensity()/2.0);

        Size trace_up_idx(apex_scan_idx);
        Size trace_down_idx(apex_scan_idx);

        // MassTrace current_trace;
        std::list<PeakType> current_trace;
        current_trace.push_back(apex_peak);

        // Initialization for the iterative version of weighted m/z mean calculation
        DoubleReal centroid_mz(apex_peak.getMZ());
        DoubleReal prev_counter(apex_peak.getIntensity() * apex_peak.getMZ());
        DoubleReal prev_denom(apex_peak.getIntensity());

        updateIterativeWeightedMeanMZ(apex_peak.getMZ(), apex_peak.getIntensity(), centroid_mz, prev_counter, prev_denom);

        std::vector<std::pair<Size, Size> > gathered_idx;
        gathered_idx.push_back(std::make_pair(apex_scan_idx, apex_peak_idx));

        Size peak_count_downward(0);
        Size peak_count_upward(0);

        Size up_hitting_peak(1), down_hitting_peak(1);
        Size up_scan_counter(0), down_scan_counter(0);

        Size fwhm_counter_down(0), fwhm_counter_up(0);

        bool toggle_up = true, toggle_down = true;

        DoubleReal int_midpoint_down(apex_peak.getIntensity()), int_midpoint_up(apex_peak.getIntensity());

        Size outliers_up(0), outliers_down(0);
        Size OUTLIER_BOUND_DOWN(10), OUTLIER_BOUND_UP(10);

        Size conseq_up(0), conseq_down(0);
        Size max_conseq_missing(5);

        DoubleReal outlier_ratio(0.3);

        DoubleReal ftl_mean(centroid_mz), ftl_sd((centroid_mz/1000000)*mass_error_ppm_);
        DoubleReal lower_sd_bound((centroid_mz/1000000)*mass_error_ppm_);

        // fixed_mean(centroid_mz), fixed_sd((centroid_mz/1000000)*mass_error_ppm_);
        DoubleReal strategy_loss(0.0);
        Size ftl_t(0);
        DoubleReal intensity_so_far(apex_peak.getIntensity());

        while (((trace_down_idx > 0) && toggle_down) || ((trace_up_idx < work_exp.size()-1) && toggle_up)) {

            // DoubleReal centroid_mz = current_trace.getCentroidMZ();

            // try to go downwards in RT
            if (((trace_down_idx > 0) && toggle_down)) {
                try
                {
                    Size next_down_peak_idx = work_exp[trace_down_idx - 1].findNearest(centroid_mz);
                    DoubleReal next_down_peak_mz = work_exp[trace_down_idx - 1][next_down_peak_idx].getMZ();
                    DoubleReal next_down_peak_int = work_exp[trace_down_idx - 1][next_down_peak_idx].getIntensity();

                    DoubleReal right_bound, left_bound;

                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;

                    //                    if (ftl_t < 5)
                    //                    {
                    //                        //                        right_bound = ftl_mean + (ftl_mean/1000000)*mass_error_ppm_;
                    //                        //                        left_bound = ftl_mean - (ftl_mean/1000000)*mass_error_ppm_;

                    //                        ftl_sd = ((ftl_mean/1000000)*mass_error_ppm_)/3.0;
                    //                        std::cout << std::setprecision(10) << "starting sd: " << ftl_sd << std::endl;
                    //                    }

                    // right_bound = centroid_mz + 3*ftl_sd;
                    // left_bound = centroid_mz - 3*ftl_sd;

                    // std::cout << "down: " << centroid_mz << " "<<  ftl_sd << std::endl;

                    Size left_next_idx = work_exp[trace_down_idx - 1].findNearest(left_bound);
                    Size right_next_idx = work_exp[trace_down_idx - 1].findNearest(right_bound);

                    DoubleReal left_mz(work_exp[trace_down_idx - 1][left_next_idx].getMZ());
                    DoubleReal right_mz(work_exp[trace_down_idx - 1][right_next_idx].getMZ());


                    //                    if ((left_next_idx < right_next_idx))
                    //                    {
                    //                        std::cout << "left: " << left_next_idx << " right: " << right_next_idx << std::endl;
                    //                    }




                    if ((next_down_peak_mz <= right_bound) && (next_down_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx])
                    {
                        Peak2D next_peak;
                        next_peak.setRT(work_exp[trace_down_idx - 1].getRT());
                        next_peak.setMZ(next_down_peak_mz);
                        next_peak.setIntensity(next_down_peak_int);

                        current_trace.push_front(next_peak);

                        updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
                        gathered_idx.push_back(std::make_pair(trace_down_idx - 1, next_down_peak_idx));

                        ++ftl_t;
                        // updateMeanEstimate(next_down_peak_mz, ftl_mean, ftl_t);
                        // updateSDEstimate(next_down_peak_mz, ftl_mean, ftl_sd, ftl_t);

                        // std::cout << "res: " << reestimate_mt_sd_ << std::endl;

                        if (reestimate_mt_sd_) {
                            // if (ftl_t > min_fwhm_scans)
                            {
                                // computeWeightedSDEstimate(current_trace, centroid_mz, ftl_sd, lower_sd_bound);
                                updateWeightedSDEstimate(next_peak, centroid_mz, ftl_sd, intensity_so_far);
                            }
                        }


                        // std::cout << "S: " << computeLoss(next_down_peak_mz, centroid_mz, ftl_sd) << " against fixed: " << computeLoss(next_down_peak_mz, fixed_mean, fixed_sd);
                        // strategy_loss += computeLoss(next_down_peak_mz, centroid_mz, ftl_sd);

                        //                        DoubleReal fixed_mean(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_mean += l_it->getMZ();
                        //                        }

                        //                        fixed_mean /= current_trace.size();

                        //                        DoubleReal fixed_sd(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_sd += (l_it->getMZ() - fixed_mean)*(l_it->getMZ() - fixed_mean);
                        //                        }

                        //                        fixed_sd = std::sqrt(fixed_sd/current_trace.size());


                        //                        DoubleReal fixed_loss(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_loss += computeLoss(l_it->getMZ(), fixed_mean, fixed_sd);
                        //                        }

                        //                        std::cout << "strategy loss: " << strategy_loss << " fixed loss: " << fixed_loss << std::endl;


                        // std::cout << (int_midpoint_down + next_down_peak_int)/2.0 << std::endl;

                        DoubleReal new_midpoint((int_midpoint_down + next_down_peak_int)/2.0);

                        if (new_midpoint > half_max_int)
                        {
                            // fwhm_down = false;
                            int_midpoint_down = new_midpoint;
                            ++fwhm_counter_down;
                        }


                        DoubleReal sd_right = centroid_mz + 3*ftl_sd;
                        DoubleReal sd_left = centroid_mz - 3*ftl_sd;


                        if (next_down_peak_mz < sd_left || next_down_peak_mz > sd_right)
                        {
                            ++outliers_down;
                        }

                        ++peak_count_downward;
                        ++down_hitting_peak;


                        conseq_down = 0;
                    }
                    else
                    {
                        ++outliers_down;

                        ++conseq_down;
                        // std::cout << "outlier down! " << std::endl;
                    }


                }
                catch(...)
                {
                    // std::cerr << "findNearest() ran into troubles..." << std::endl;
                }
                --trace_down_idx;
                ++down_scan_counter;


                //                if (down_scan_counter > min_data_points) {
                //                    DoubleReal sample_rate_down = (DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter;

                //                if (outlier_ratio*DoubleReal(peak_count_downward) > OUTLIER_BOUND_DOWN)
                //                {
                //                    OUTLIER_BOUND_DOWN = outlier_ratio*DoubleReal(peak_count_downward);
                //                }

                //                if (outliers_down > OUTLIER_BOUND_DOWN) {
                //                    toggle_down = false;
                //                }

                if (conseq_down > max_conseq_missing) {
                    toggle_down = false;
                }
                //                }
            }

            //}
            // *********************************************************** //
            // MOVE UP in RT dim
            // *********************************************************** //

            if (((trace_up_idx < work_exp.size()-1) && toggle_up))
            {
                try
                {
                    Size next_up_peak_idx = work_exp[trace_up_idx + 1].findNearest(centroid_mz);
                    DoubleReal next_up_peak_mz = work_exp[trace_up_idx + 1][next_up_peak_idx].getMZ();
                    DoubleReal next_up_peak_int = work_exp[trace_up_idx + 1][next_up_peak_idx].getIntensity();

                    DoubleReal right_bound, left_bound;

                    right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                    left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;



                    //                    if (ftl_t < 5)
                    //                    {
                    //                        //                        right_bound = ftl_mean + (ftl_mean/1000000)*mass_error_ppm_;
                    //                        //                        left_bound = ftl_mean - (ftl_mean/1000000)*mass_error_ppm_;

                    //                        ftl_sd = ((ftl_mean/1000000)*mass_error_ppm_)/3.0;
                    //                        std::cout << std::setprecision(10) << "starting sd: " << ftl_sd << std::endl;
                    //                    }

                    // right_bound = centroid_mz + 3*ftl_sd;
                    // left_bound = centroid_mz - 3*ftl_sd;

                    // std::cout << "up: " << centroid_mz << " " << ftl_sd << std::endl;

                    if ((next_up_peak_mz <= right_bound) && (next_up_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx])
                    {
                        Peak2D next_peak;
                        next_peak.setRT(work_exp[trace_up_idx + 1].getRT());
                        next_peak.setMZ(next_up_peak_mz);
                        next_peak.setIntensity(next_up_peak_int);

                        current_trace.push_back(next_peak);

                        updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
                        gathered_idx.push_back(std::make_pair(trace_up_idx + 1, next_up_peak_idx));

                        ++ftl_t;
                        // updateMeanEstimate(next_up_peak_mz, ftl_mean, ftl_t);
                        // updateSDEstimate(next_up_peak_mz, ftl_mean, ftl_sd, ftl_t);

                        if (reestimate_mt_sd_) {
                            // if (ftl_t > min_fwhm_scans)
                            {
                                //computeWeightedSDEstimate(current_trace, centroid_mz, ftl_sd, lower_sd_bound);
                                updateWeightedSDEstimate(next_peak, centroid_mz, ftl_sd, intensity_so_far);
                            }
                        }


                        // strategy_loss += computeLoss(next_up_peak_mz, centroid_mz, ftl_sd);

                        //                        DoubleReal fixed_mean(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_mean += l_it->getMZ();
                        //                        }

                        //                        fixed_mean /= current_trace.size();

                        //                        DoubleReal fixed_sd(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_sd += (l_it->getMZ() - fixed_mean)*(l_it->getMZ() - fixed_mean);
                        //                        }

                        //                        fixed_sd = std::sqrt(fixed_sd/current_trace.size());


                        //                        DoubleReal fixed_loss(0.0);

                        //                        for (std::list<PeakType>::const_iterator l_it = current_trace.begin(); l_it != current_trace.end(); ++l_it)
                        //                        {
                        //                            fixed_loss += computeLoss(l_it->getMZ(), fixed_mean, fixed_sd);
                        //                        }

                        //                        std::cout << "strategy loss: " << strategy_loss << " fixed loss: " << fixed_loss << std::endl;


                        DoubleReal new_midpoint((int_midpoint_up + next_up_peak_int)/2.0);

                        if (new_midpoint > half_max_int)
                        {
                            // fwhm_up = false;
                            int_midpoint_up = new_midpoint;
                            ++fwhm_counter_up;
                        }

                        ++peak_count_upward;
                        ++up_hitting_peak;

                        DoubleReal sd_right = centroid_mz + 3*ftl_sd;
                        DoubleReal sd_left = centroid_mz - 3*ftl_sd;


                        if (next_up_peak_mz < sd_left || next_up_peak_mz > sd_right)
                        {
                            ++outliers_up;
                        }

                        conseq_up = 0;

                    }
                    else
                    {
                        ++outliers_up;
                        ++conseq_up;
                        // std::cout << "outlier up" << std::endl;
                    }

                }
                catch(...)
                {
                    //  std::cerr << "findNearest() ran into troubles..." << std::endl;
                }

                ++trace_up_idx;
                ++up_scan_counter;


                //                if (up_scan_counter > min_data_points) {
                //                    DoubleReal sample_rate_up = (DoubleReal)up_hitting_peak/(DoubleReal)up_scan_counter;

                //                if (outlier_ratio*DoubleReal(peak_count_upward) > OUTLIER_BOUND_UP)
                //                {
                //                    OUTLIER_BOUND_UP = outlier_ratio*DoubleReal(peak_count_upward);
                //                }


                //                if (outliers_up > OUTLIER_BOUND_UP)
                //                {
                //                    toggle_up = false;
                //                }
                //                }

                if (conseq_up > max_conseq_missing)
                {
                    toggle_up = false;
                }


            }

        }

        DoubleReal num_scans(down_scan_counter+up_scan_counter+1 - conseq_up - conseq_down);


        //if (current_trace.size() >= min_fwhm_scans)
        if (current_trace.size() >= std::floor(num_scans*min_sample_rate_) && num_scans >= min_fwhm_scans)
        {

            // std::cout << "minscans: " << num_scans << " peaks: " << current_trace.size() << std::endl;

            // std::cout << "CURR: " << current_trace.size() << " " << down_scan_counter + up_scan_counter +1 << std::endl;
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
            new_trace.setFWHMScansNum(min_fwhm_scans);

            peaks_detected += new_trace.getSize();
            this->setProgress(peaks_detected);
            found_masstraces.push_back(new_trace);
            ++trace_number;
        }
    }

    // std::vector<MassTrace> tmp_mt_vec;

    // filterByPeakWidth(found_masstraces, tmp_mt_vec);

    // std::cout << "result: " << found_masstraces.size() << " filt: " << tmp_mt_vec.size() << std::endl;
    // found_masstraces = tmp_mt_vec;


    this->endProgress();
    // std::cout << found_masstraces.size() << " traces found" << std::endl;

    return ;
} // end of MassTraceDetection::run

void MassTraceDetection::updateMembers_()
{
    mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
    noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");

    // chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");

    min_sample_rate_ = (DoubleReal)param_.getValue("min_sample_rate");
    min_peak_width_ = (DoubleReal)param_.getValue("min_peak_width");
    reestimate_mt_sd_ = param_.getValue("reestimate_mt_sd").toBool();
}

}
