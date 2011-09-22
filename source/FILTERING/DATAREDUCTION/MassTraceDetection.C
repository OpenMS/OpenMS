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

#include <vector>
#include <map>
#include <algorithm>
#include <sstream>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
    MassTraceDetection::MassTraceDetection()
        : DefaultParamHandler("MassTraceDetection"), ProgressLogger()
    {
        // defaults_.setValue( "name" , 1 , "descript" );
        defaults_.setValue( "mass_error_ppm" , 20.0 , "Allowed mass deviation (in ppm)");
        defaults_.setValue( "noise_threshold_int" , 10.0 , "Intensity threshold below which peaks are removed as noise");
        defaults_.setValue( "chrom_apex_snt" , 3.0 , "Minimum signal-to-noise a mass trace should have");
        defaults_.setValue( "chrom_fwhm" , 10.0 , "Minimum FWHM (in seconds) a chromatographic peak should have");
        defaults_.setValue( "min_sample_rate" , 0.5 , "Minimum sampling rate of a mass trace");

        defaultsToParam_();

        this->setLogType(CMD);
    }

    MassTraceDetection::~MassTraceDetection()
    {

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

        // this->startProgress(0, input_exp.size(), "Detect potential chromatographic apeces...");
        for (Size scan_idx = 0; scan_idx < input_exp.size(); ++scan_idx)
        {
            // this->setProgress(scan_idx);

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
        }

        // this->endProgress();

        // discard last spectrum's offset
        spec_offsets.pop_back();

        boost::dynamic_bitset<> peak_visited(peak_count);

        //    std::cout << "size work_exp: " << peak_count << std::endl;

        // start extending mass traces beginning with the apex peak

        DoubleReal scan_rt_diff = (input_exp[input_exp.size() - 1].getRT() - input_exp[0].getRT())/(input_exp.size());
        Size min_datapoints = std::floor(chrom_fwhm_/(scan_rt_diff*2));
        Size min_trace_quality = std::floor(2*min_datapoints*min_sample_rate_);

        // std::vector<MassTrace> found_mtraces;

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

            // DoubleReal apex_intensity(m_it->first);

            // std::cout << "int: " << apex_intensity << std::endl;

            Size trace_up_idx(apex_scan_idx);
            Size trace_down_idx(apex_scan_idx);

            MassTrace current_trace;

            current_trace.appendPeak(apex_peak);


            std::vector<std::pair<Size, Size> > gathered_idx;
            gathered_idx.push_back(std::make_pair(apex_scan_idx, apex_peak_idx));

            Size peak_count_downward(0);
            Size peak_count_upward(0);

            Size up_hitting_peak(1), down_hitting_peak(1);
            Size up_scan_counter(1), down_scan_counter(1);

            bool toggle_up = true, toggle_down = true;
            bool is_valid = false;

            while (((trace_down_idx > 0) && toggle_down) || ((trace_up_idx < work_exp.size()-1) && toggle_up)) {

                DoubleReal centroid_mz = current_trace.getCentroidMZ();

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

                            current_trace.prependPeak(next_peak);
                            gathered_idx.push_back(std::make_pair(trace_down_idx - 1, next_down_peak_idx));

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

                    if (down_scan_counter > min_datapoints) {
                        DoubleReal sample_rate_down = (DoubleReal)down_hitting_peak/(DoubleReal)down_scan_counter;

                        if (sample_rate_down < min_sample_rate_) {
                            toggle_down = false;
                        }
                        else {
                            is_valid = true;
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

                            current_trace.appendPeak(next_peak);
                            gathered_idx.push_back(std::make_pair(trace_up_idx + 1, next_up_peak_idx));

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

                    if (up_scan_counter > min_datapoints) {
                        DoubleReal sample_rate_up = (DoubleReal)up_hitting_peak/(DoubleReal)up_scan_counter;

                        if (sample_rate_up < min_sample_rate_) {
                            toggle_up = false;
                        }
                        else {
                            is_valid = true;
                        }

                    }

                }

            }

            if (current_trace.getSize() >= min_trace_quality) {
                // mark all peaks as visited
                for (Size i = 0; i < gathered_idx.size(); ++i)
                {
                    peak_visited[spec_offsets[gathered_idx[i].first] +  gathered_idx[i].second] = true;
                }

                String tr_num;
                std::stringstream read_in;
                read_in << trace_number;
                tr_num = read_in.str();

                current_trace.setLabel("T" + tr_num);
                peaks_detected += current_trace.getSize();
                this->setProgress(peaks_detected);
                found_masstraces.push_back(current_trace);
                ++trace_number;
            }
        }

        this->endProgress();
        // std::cout << found_masstraces.size() << " traces found" << std::endl;


        return ;
    } // end of MassTraceDetection::run






    void MassTraceDetection::updateMembers_()
    {
        // delta_ = (Size)param_.getValue( "delta" );

        mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
        noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
        chrom_apex_snt_ = (DoubleReal)param_.getValue("chrom_apex_snt");
        chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");

        min_sample_rate_ = (DoubleReal)param_.getValue("min_sample_rate");
    }

}
