// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
    defaults_.setValue("chrom_apex_snt" , 3.0 , "Minimum signal-to-noise a mass trace should have.");
    defaults_.setValue("chrom_fwhm" , 0.0 , "Allows filtering of mass traces with peak width (in seconds) less than this threshold. Disabled by default (set to 0.0).", StringList::create("advanced"));
    defaults_.setValue("min_sample_rate" , 0.5 , "Minimum fraction of scans along the mass trace that must contain a peak.", StringList::create("advanced"));

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

        //        Size pw_vec_size = pw_vec.size();
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
        return ;
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

                    if (tmp_peak_int > chrom_apex_snt_*noise_threshold_int_)
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

    if (spectra_count < 5)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input map consists of too few spectra (less than 5!). Aborting...", String(spectra_count));
    }

        // this->endProgress();

        // discard last spectrum's offset
        spec_offsets.pop_back();

        boost::dynamic_bitset<> peak_visited(peak_count);

        // start extending mass traces beginning with the apex peak

        // Size min_datapoints = std::floor(chrom_fwhm_/(scan_rt_diff*2));
    Size min_data_points(2);

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
            Size up_scan_counter(1), down_scan_counter(1);

            Size fwhm_counter_down(0), fwhm_counter_up(0);

            bool toggle_up = true, toggle_down = true;

            DoubleReal int_midpoint_down(apex_peak.getIntensity()), int_midpoint_up(apex_peak.getIntensity());

            while (((trace_down_idx > 0) && toggle_down) || ((trace_up_idx < work_exp.size()-1) && toggle_up)) {

            // DoubleReal centroid_mz = current_trace.getCentroidMZ();

                // try to go downwards in RT
                if (((trace_down_idx > 0) && toggle_down)) {
                    try
                    {
                        Size next_down_peak_idx = work_exp[trace_down_idx - 1].findNearest(centroid_mz);
                        DoubleReal next_down_peak_mz = work_exp[trace_down_idx - 1][next_down_peak_idx].getMZ();
                        DoubleReal next_down_peak_int = work_exp[trace_down_idx - 1][next_down_peak_idx].getIntensity();

                        DoubleReal right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                        DoubleReal left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;


                        if ((next_down_peak_mz <= right_bound) && (next_down_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_down_idx - 1] + next_down_peak_idx]) {
                            Peak2D next_peak;
                            next_peak.setRT(work_exp[trace_down_idx - 1].getRT());
                            next_peak.setMZ(next_down_peak_mz);
                            next_peak.setIntensity(next_down_peak_int);

                        current_trace.push_front(next_peak);

                        updateIterativeWeightedMeanMZ(next_down_peak_mz, next_down_peak_int, centroid_mz, prev_counter, prev_denom);
                            gathered_idx.push_back(std::make_pair(trace_down_idx - 1, next_down_peak_idx));

                            // std::cout << (int_midpoint_down + next_down_peak_int)/2.0 << std::endl;

                            DoubleReal new_midpoint((int_midpoint_down + next_down_peak_int)/2.0);

                            if (new_midpoint > half_max_int)
                            {
                                // fwhm_down = false;
                                int_midpoint_down = new_midpoint;
                                ++fwhm_counter_down;
                            }

                            ++peak_count_downward;
                            ++down_hitting_peak;
                        }


                    }
                    catch(...)
                    {
                        // std::cerr << "findNearest() ran into troubles..." << std::endl;
                    }
                    --trace_down_idx;
                    ++down_scan_counter;


                    if (down_scan_counter > min_data_points) {
                        DoubleReal sample_rate_down = (DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter;

                        if (sample_rate_down < min_sample_rate_) {
                            toggle_down = false;
                        }
                        }
                    }

                //}
                // *********************************************************** //
                // MOVE UP in RT dim
                // *********************************************************** //

                if (((trace_up_idx < work_exp.size()-1) && toggle_up)) {
                    try
                    {
                        Size next_up_peak_idx = work_exp[trace_up_idx + 1].findNearest(centroid_mz);
                        DoubleReal next_up_peak_mz = work_exp[trace_up_idx + 1][next_up_peak_idx].getMZ();
                        DoubleReal next_up_peak_int = work_exp[trace_up_idx + 1][next_up_peak_idx].getIntensity();

                        DoubleReal right_bound = centroid_mz + (centroid_mz/1000000)*mass_error_ppm_;
                        DoubleReal left_bound = centroid_mz - (centroid_mz/1000000)*mass_error_ppm_;

                        if ((next_up_peak_mz <= right_bound) && (next_up_peak_mz >= left_bound) && !peak_visited[spec_offsets[trace_up_idx + 1] + next_up_peak_idx]) {
                            Peak2D next_peak;
                            next_peak.setRT(work_exp[trace_up_idx + 1].getRT());
                            next_peak.setMZ(next_up_peak_mz);
                            next_peak.setIntensity(next_up_peak_int);

                        current_trace.push_back(next_peak);

                        updateIterativeWeightedMeanMZ(next_up_peak_mz, next_up_peak_int, centroid_mz, prev_counter, prev_denom);
                            gathered_idx.push_back(std::make_pair(trace_up_idx + 1, next_up_peak_idx));

                            DoubleReal new_midpoint((int_midpoint_up + next_up_peak_int)/2.0);

                            if (new_midpoint > half_max_int)
                            {
                                // fwhm_up = false;
                                int_midpoint_up = new_midpoint;
                                ++fwhm_counter_up;
                            }

                            ++peak_count_upward;
                            ++up_hitting_peak;

                        }

                    }
                    catch(...)
                    {
                        //  std::cerr << "findNearest() ran into troubles..." << std::endl;
                    }

                    ++trace_up_idx;
                    ++up_scan_counter;


                    if (up_scan_counter > min_data_points) {
                        DoubleReal sample_rate_up = (DoubleReal)up_hitting_peak/(DoubleReal)up_scan_counter;

                    if (sample_rate_up < min_sample_rate_)
                    {
                            toggle_up = false;
                        }
                        }

                    }

                }


        if (current_trace.size() >= 2*min_data_points + 1)
            {

            // std::cout << "CURR: " << current_trace.size() << " " << fwhm_counter_up + fwhm_counter_down+1 << std::endl;
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
            MassTrace new_trace(current_trace);
            new_trace.updateWeightedMeanRT();
            new_trace.updateWeightedMeanMZ();

            new_trace.setLabel("T" + tr_num);
            new_trace.setFWHMScansNum(fwhm_counter_down + fwhm_counter_up + 1);

            peaks_detected += new_trace.getSize();
                this->setProgress(peaks_detected);
            found_masstraces.push_back(new_trace);
                ++trace_number;
            }
        }

    std::vector<MassTrace> tmp_mt_vec;

    filterByPeakWidth(found_masstraces, tmp_mt_vec);

    // std::cout << "result: " << found_masstraces.size() << " filt: " << tmp_mt_vec.size() << std::endl;
    found_masstraces = tmp_mt_vec;


        this->endProgress();
        // std::cout << found_masstraces.size() << " traces found" << std::endl;

        return ;
    } // end of MassTraceDetection::run

    void MassTraceDetection::updateMembers_()
    {
        mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
        noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
        chrom_apex_snt_ = (DoubleReal)param_.getValue("chrom_apex_snt");
        chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");

        min_sample_rate_ = (DoubleReal)param_.getValue("min_sample_rate");
    }

}
