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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------


#include "OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h"
#include "OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h"
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <sstream>
#include <numeric>
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
ElutionPeakDetection::ElutionPeakDetection() :
    DefaultParamHandler("ElutionPeakDetection"), ProgressLogger()
{
    defaults_.setValue("chrom_fwhm", 5.0, "Expected full-width-at-half-maximum of chromatographic peaks.");
    defaults_.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");

    defaults_.setValue("width_filtering", "fixed", "Enable filtering of unlikely peak widths. The fixed setting filters out mass traces outside the [min_fwhm, max_fwhm] interval (set parameters accordingly!). The auto setting filters with the 5 and 95% quantiles of the peak width distribution.");
    defaults_.setValidStrings("width_filtering", StringList::create(("off,fixed,auto")));
    defaults_.setValue("min_fwhm", 3.0, "Minimum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if paramter width_filtering is off or auto.", StringList::create("advanced"));
    defaults_.setValue("max_fwhm", 60.0, "Maximum full-width-at-half-maximum of chromatographic peaks (in seconds). Ignored if paramter width_filtering is off or auto.", StringList::create("advanced"));

    defaults_.setValue("masstrace_snr_filtering", "false", "Apply post-filtering by signal-to-noise ratio after smoothing.", StringList::create("advanced"));
    defaults_.setValidStrings("masstrace_snr_filtering", StringList::create(("false,true")));

    // defaults_.setValue("min_trace_length", 5.0, "Minimum length of a mass trace (in seconds).", StringList::create("advanced"));
    // defaults_.setValue("max_trace_length", 300.0, "Maximum length of a mass trace (in seconds).", StringList::create("advanced"));


    defaultsToParam_();

    this->setLogType(CMD);
}

ElutionPeakDetection::~ElutionPeakDetection()
{
}


DoubleReal ElutionPeakDetection::computeMassTraceNoise(const MassTrace& tr)
{
    // compute RMSE
    DoubleReal squared_sum(0.0);
    std::vector<DoubleReal> smooth_ints(tr.getSmoothedIntensities());

    for (Size i = 0; i < smooth_ints.size(); ++i)
    {
        squared_sum += (tr[i].getIntensity() - smooth_ints[i])*(tr[i].getIntensity() - smooth_ints[i]);
    }

    DoubleReal rmse(0.0);

    if (smooth_ints.size() > 0)
    {
        rmse = std::sqrt(squared_sum/smooth_ints.size());
    }

    return rmse;
}

DoubleReal ElutionPeakDetection::computeMassTraceSNR(const MassTrace& tr)
{
    DoubleReal noise_area(1.0), signal_area(0.0), snr(0.0);

    if (tr.getSize() > 0)
    {
        noise_area = computeMassTraceNoise(tr) * tr.getTraceLength();
        signal_area = tr.computePeakArea();

        snr = signal_area/noise_area;
    }

    // std::cout << "snr " << snr << " ";

    return snr;
}

DoubleReal ElutionPeakDetection::computeApexSNR(const MassTrace& tr)
{
    DoubleReal snr(0.0);
    DoubleReal noise_level(computeMassTraceNoise(tr));
    DoubleReal smoothed_apex_int(tr.getMaxIntensity(true));

    if (noise_level > 0.0)
    {
        snr = smoothed_apex_int/noise_level;
    }

    // std::cout << "snr " << snr << " ";

    return snr;
}

void ElutionPeakDetection::findLocalExtrema(const MassTrace& tr, const Size & num_neighboring_peaks, std::vector<Size> & chrom_maxes, std::vector<Size> & chrom_mins)
{
    std::vector<DoubleReal> smoothed_ints_vec(tr.getSmoothedIntensities());

    Size mt_length(smoothed_ints_vec.size());

    if (mt_length != tr.getSize())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_ints_vec.size()));
    }

    // first make sure that everything is cleared
    chrom_maxes.clear();
    chrom_mins.clear();

    // Extract RTs from the chromatogram and store them into into vectors for index access

    // std::cout << "neighboring peaks: " << num_neighboring_peaks << std::endl;

    //  Store indices along with smoothed_ints to keep track of the peak order
    std::multimap<DoubleReal, Size> intensity_indices;
    boost::dynamic_bitset<> used_idx(mt_length);

    for (Size i = 0; i < mt_length; ++i)
    {
        intensity_indices.insert(std::make_pair(smoothed_ints_vec[i], i));
    }


    for (std::multimap<DoubleReal, Size>::const_iterator c_it = intensity_indices.begin(); c_it != intensity_indices.end(); ++c_it)
    {
        DoubleReal ref_int = c_it->first;
        Size ref_idx = c_it->second;

        if (!(used_idx[ref_idx]) && ref_int > 0.0)
        {
            bool real_max = true;

            // iterate up the RT
            Size start_idx(0);

            if (ref_idx > num_neighboring_peaks)
            {
                start_idx = ref_idx - num_neighboring_peaks;
            }

            Size end_idx = ref_idx + num_neighboring_peaks;

            if (end_idx > mt_length)
            {
                end_idx = mt_length;
            }

            for (Size j = start_idx; j < end_idx; ++j)
            {
                if (used_idx[j])
                {
                    real_max = false;
                    break;
                }

                if (j == ref_idx)
                {
                    continue;
                }

                if (smoothed_ints_vec[j] > ref_int)
                {
                    real_max = false;
                }
            }

            if (real_max)
            {
                chrom_maxes.push_back(ref_idx);

                for (Size j = start_idx; j < end_idx; ++j)
                {
                    used_idx[j] = true;
                }
            }

        }
    }


    std::sort(chrom_maxes.begin(), chrom_maxes.end());


    if (chrom_maxes.size() > 1)
    {

        Size i(0), j(1);
        //for (Size i = 0; i < chrom_maxes.size() - 1; ++i)

        while (i < j && j < chrom_maxes.size())
        {
            // bisection
            Size left_bound(chrom_maxes[i] + 1);
            Size right_bound(chrom_maxes[j] - 1);

            while ((left_bound + 1) < right_bound)
            {
                DoubleReal mid_dist((right_bound - left_bound) / 2.0);

                Size mid_element_idx(left_bound + std::floor(mid_dist));

                DoubleReal mid_element_int = smoothed_ints_vec[mid_element_idx];

                if (mid_element_int <= smoothed_ints_vec[mid_element_idx + 1])
                {
                    right_bound = mid_element_idx;
                }
                else       // or to the right...
                {
                    left_bound = mid_element_idx;
                }

            }

            Size min_rt((smoothed_ints_vec[left_bound] < smoothed_ints_vec[right_bound]) ? left_bound : right_bound);

            // check for valley depth between chromatographic peaks
            DoubleReal min_int(1.0);
            if (smoothed_ints_vec[min_rt] > min_int)
            {
                min_int = smoothed_ints_vec[min_rt];
            }

            DoubleReal left_max_int(smoothed_ints_vec[chrom_maxes[i]]);
            DoubleReal right_max_int(smoothed_ints_vec[chrom_maxes[j]]);

            DoubleReal left_rt(tr[chrom_maxes[i]].getRT());
            DoubleReal mid_rt(tr[min_rt].getRT());
            DoubleReal right_rt(tr[chrom_maxes[j]].getRT());

            DoubleReal left_dist(std::fabs(mid_rt - left_rt));
            DoubleReal right_dist(std::fabs(right_rt - mid_rt));
            DoubleReal min_dist(min_fwhm_ / 2.0);

            // out debug info
            // std::cout << tr.getLabel() << ": i,j " << i << "," << j << ":" << left_max_int << " min: " << min_int << " " << right_max_int << " l " << left_rt << " r " << right_rt << " m " << mid_rt << std::endl;



            if (left_max_int / min_int >= 2.0
                    && right_max_int / min_int >= 2.0
                    && left_dist >= min_dist
                    && right_dist >= min_dist)
            {
                chrom_mins.push_back(min_rt);

                // std::cout << "min added!" << std::endl;
                i = j;
                ++j;
            }
            else
            {
                // keep one of the chrom_maxes, iterate the other
                if (left_max_int > right_max_int)
                {
                    ++j;
                }
                else
                {
                    i = j;
                    ++j;
                }
            }

            // chrom_mins.push_back(min_rt);
        }
    }

    return;
}


void ElutionPeakDetection::detectPeaks(MassTrace & mt, std::vector<MassTrace> & single_mtraces)
{
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    detectElutionPeaks_(mt, single_mtraces);
    return;
}

void ElutionPeakDetection::detectPeaks(std::vector<MassTrace> & mt_vec, std::vector<MassTrace> & single_mtraces)
{
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    this->startProgress(0, mt_vec.size(), "elution peak detection");
    Size progress(0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < (SignedSize) mt_vec.size(); ++i)
    {
        IF_MASTERTHREAD this->setProgress(progress);
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++progress;

        detectElutionPeaks_(mt_vec[i], single_mtraces);
    }

    this->endProgress();

    return;
}

void ElutionPeakDetection::filterByPeakWidth(std::vector<MassTrace> & mt_vec, std::vector<MassTrace> & filt_mtraces)
{
    filt_mtraces.clear();

    std::multimap<DoubleReal, Size> sorted_by_peakwidth;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        sorted_by_peakwidth.insert(std::make_pair(mt_vec[i].estimateFWHM(true), i));
    }

    DoubleReal mapsize(sorted_by_peakwidth.size());
    Size lower_quartile_idx(std::floor(mapsize * 0.05));
    Size upper_quartile_idx(std::floor(mapsize * 0.95));
    Size count_mt(0);

    // filter out mass traces below lower quartile and above upper quartile
    for (std::multimap<DoubleReal, Size>::const_iterator m_it = sorted_by_peakwidth.begin(); m_it != sorted_by_peakwidth.end(); ++m_it)
    {
        if (count_mt >= lower_quartile_idx && count_mt <= upper_quartile_idx)
        {
            // std::cout << "pw added " << m_it->first << std::endl;
            filt_mtraces.push_back(mt_vec[m_it->second]);
        }
        ++count_mt;
    }

    std::cout << "pw low: " << filt_mtraces[0].estimateFWHM(true) << " " << " pw high: " << filt_mtraces[filt_mtraces.size() - 1].estimateFWHM(true) << std::endl;

    return;
}

void ElutionPeakDetection::detectElutionPeaks_(MassTrace & mt, std::vector<MassTrace> & single_mtraces)
{
    std::vector<DoubleReal> rts, ints;

    for (MassTrace::const_iterator c_it = mt.begin(); c_it != mt.end(); ++c_it)
    {
        rts.push_back(c_it->getRT());
        ints.push_back(c_it->getIntensity());
    }

    std::vector<DoubleReal> smoothed_data;


    LowessSmoothing lowess_smooth;
    Param lowess_params;

    // use dynamically computed window sizes

    // Size win_size = mt.getFWHMScansNum();

    // use one global window size for all mass traces to smooth
    DoubleReal scan_time(mt.getScanTime());
    Size win_size = std::ceil(chrom_fwhm_ / scan_time);

    // std::cout << "win_size elution: " << scan_time << " " << win_size << std::endl;

    // if there is no previous FWHM estimation... do it now
    //    if (win_size == 0)
    //    {
    //        mt.estimateFWHM(false); // estimate FWHM
    //        win_size = mt.getFWHMScansNum();
    //    }

    lowess_params.setValue("window_size", win_size);
    lowess_smooth.setParameters(lowess_params);

    lowess_smooth.smoothData(rts, ints, smoothed_data);

    mt.setSmoothedIntensities(smoothed_data);

    // debug intensities

    // Size i = 0;

    //    std::cout << "*****" << std::endl;
    //    for (MassTrace::const_iterator mt_it = mt.begin(); mt_it != mt.end(); ++mt_it)
    //    {
    //        std::cout << mt_it->getIntensity() << " " << smoothed_data[i] << std::endl;
    //        ++i;
    //    }
    //std::cout << "*****" << std::endl;

    std::vector<Size> maxes, mins;

    // mt.findLocalExtrema(win_size / 2, maxes, mins);

    findLocalExtrema(mt, win_size/2, maxes, mins);

    // if only one maximum exists: finished!
    if (maxes.size() == 1)
    {
        bool pw_ok = true;
        bool snr_ok = true;

        // check mass trace filter criteria (if enabled)
        if (pw_filtering_ == "fixed")
        {
            DoubleReal act_fwhm(mt.estimateFWHM(true));

            // std::cout << "act_fwhm: " << act_fwhm << " ";

            if (act_fwhm < min_fwhm_ || act_fwhm > max_fwhm_)
            {
                pw_ok = false;
            }

            // std::cout << pw_ok << std::endl;
        }

        if (mt_snr_filtering_)
        {
            if (computeApexSNR(mt) < chrom_peak_snr_)
            {
                snr_ok = false;
            }
        }


        if (pw_ok && snr_ok)
        {
            mt.updateSmoothedMaxRT();

            if (pw_filtering_ != "fixed")
            {
                mt.estimateFWHM(true);
            }

            // check for minimum/maximum trace length
            //          DoubleReal mt_length(std::fabs(mt.rbegin()->getRT() - mt.begin()->getRT()));

            //        if ((mt_length >= min_trace_length_) && (mt_length <= max_trace_length_))
            // if (mt_quality >= 1.2)
            //      {
#ifdef _OPENMP
#pragma omp critical
#endif
            single_mtraces.push_back(mt);

        }
    }
    else if (maxes.empty())
    {
        return;
    }
    else // split mt to subtraces
    {
        MassTrace::const_iterator cp_it = mt.begin();
        Size last_idx(0);

        for (Size min_idx = 0; min_idx < mins.size(); ++min_idx)
        {
            // copy subtrace between cp_it and splitpoint
            std::vector<PeakType> tmp_mt;
            std::vector<DoubleReal> smoothed_tmp;

            while (last_idx <= mins[min_idx])
            {
                tmp_mt.push_back(*cp_it);
                smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
                ++cp_it;
                ++last_idx;
            }

            // check if

//            if (tmp_mt.size() >= win_size / 2)
//            {
                DoubleReal scantime(mt.getScanTime());
                MassTrace new_mt(tmp_mt, scantime);

                // copy smoothed ints
                new_mt.setSmoothedIntensities(smoothed_tmp);


                // check filter criteria
                bool pw_ok = true;
                bool snr_ok = true;

                // check mass trace filter criteria (if enabled)
                if (pw_filtering_ == "fixed")
                {
                    DoubleReal act_fwhm(new_mt.estimateFWHM(true));

                    // std::cout << "act_fwhm: " << act_fwhm << " ";

                    if (act_fwhm < min_fwhm_ || act_fwhm > max_fwhm_)
                    {
                        pw_ok = false;
                    }

                    // std::cout << pw_ok << std::endl;
                }

                if (mt_snr_filtering_)
                {
                    if (computeApexSNR(mt) < chrom_peak_snr_)
                    {
                        snr_ok = false;
                    }
                }


                if (pw_ok && snr_ok)
                {

                    // set label of subtrace
                    String tr_num;
                    std::stringstream read_in;
                    read_in << (min_idx + 1);
                    tr_num = "." + read_in.str();

                    new_mt.setLabel(mt.getLabel() + tr_num);
                    //new_mt.updateWeightedMeanRT();
                    new_mt.updateSmoothedMaxRT();
                    //new_mt.updateSmoothedWeightedMeanRT();
                    new_mt.updateWeightedMeanMZ();
                    new_mt.updateWeightedMZsd();

                    if (pw_filtering_ != "fixed")
                    {
                        new_mt.estimateFWHM(true);
                    }
                    // DoubleReal mt_quality(computeApexSNR(new_mt));

                    // DoubleReal new_mt_length(std::fabs(new_mt.rbegin()->getRT() - new_mt.begin()->getRT()));

                    // if ((new_mt_length >= min_trace_length_) && (new_mt_length <= max_trace_length_))
                    //{
#ifdef _OPENMP
#pragma omp critical
#endif
                    single_mtraces.push_back(new_mt);
                }
          //  }
        }

        // don't forget the trailing trace
        std::vector<PeakType> tmp_mt;

        std::vector<DoubleReal> smoothed_tmp;

        while (last_idx < mt.getSize())
        {
            tmp_mt.push_back(*cp_it);
            smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
            ++cp_it;
            ++last_idx;
        }

//        if (tmp_mt.size() >= win_size / 2)
//        {
            DoubleReal scantime(mt.getScanTime());
            MassTrace new_mt(tmp_mt, scantime);

            // copy smoothed ints
            new_mt.setSmoothedIntensities(smoothed_tmp);

            // check filter criteria
            bool pw_ok = true;
            bool snr_ok = true;

            // check mass trace filter criteria (if enabled)
            if (pw_filtering_ == "fixed")
            {
                DoubleReal act_fwhm(new_mt.estimateFWHM(true));

                // std::cout << "act_fwhm: " << act_fwhm << " ";

                if (act_fwhm < min_fwhm_ || act_fwhm > max_fwhm_)
                {
                    pw_ok = false;
                }

                // std::cout << pw_ok << std::endl;
            }

            if (mt_snr_filtering_)
            {
                if (computeApexSNR(mt) < chrom_peak_snr_)
                {
                    snr_ok = false;
                }
            }


            if (pw_ok && snr_ok)
            {
                // set label of subtrace
                String tr_num;
                std::stringstream read_in;
                read_in << (mins.size() + 1);
                tr_num = "." + read_in.str();

                new_mt.setLabel(mt.getLabel() + tr_num);
                new_mt.updateSmoothedMaxRT();
                new_mt.updateWeightedMeanMZ();
                new_mt.updateWeightedMZsd();

                if (pw_filtering_ != "fixed")
                {
                    new_mt.estimateFWHM(true);
                }
                // DoubleReal mt_quality(computeApexSNR(new_mt));

                //                DoubleReal mt_length(std::fabs(new_mt.rbegin()->getRT() - new_mt.begin()->getRT()));

                //                if ((mt_length >= min_trace_length_) && (mt_length <= max_trace_length_))
                //                {
#ifdef _OPENMP
#pragma omp critical
#endif
                single_mtraces.push_back(new_mt);
            }
     //   }
    }
    return;
}

void ElutionPeakDetection::updateMembers_()
{
    chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");
    chrom_peak_snr_ = (DoubleReal)param_.getValue("chrom_peak_snr");
    noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
    // min_trace_length_ = (DoubleReal)param_.getValue("min_trace_length");
    // max_trace_length_ = (DoubleReal)param_.getValue("max_trace_length");
    min_fwhm_ = (DoubleReal)param_.getValue("min_fwhm");
    max_fwhm_ = (DoubleReal)param_.getValue("max_fwhm");

    pw_filtering_ = param_.getValue("width_filtering");
    mt_snr_filtering_ = param_.getValue("masstrace_snr_filtering").toBool();
}

} //namespace OpenMS
