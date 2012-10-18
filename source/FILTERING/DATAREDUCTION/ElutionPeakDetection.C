// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

    defaults_.setValue("width_filtering", "false", "Enable filtering of unlikely peak widths (5 and 95% quantiles of peak width distribution).");
    defaults_.setValidStrings("width_filtering", StringList::create(("false,true")));

    defaults_.setValue("masstrace_snr_filtering", "false", "Apply post-filtering by signal-to-noise ratio after smoothing.", StringList::create("advanced"));
    defaults_.setValidStrings("masstrace_snr_filtering", StringList::create(("false,true")));

    defaultsToParam_();

    this->setLogType(CMD);
}

ElutionPeakDetection::~ElutionPeakDetection()
{
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

void ElutionPeakDetection::estimatePeakWidth(std::vector<MassTrace> & mt_vec)
{
    std::multimap<DoubleReal, Size> histo_map;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        DoubleReal fwhm(mt_vec[i].estimateFWHM(false));

        if (fwhm > 0.0)
        {
            histo_map.insert(std::make_pair(fwhm, i));
        }
    }

    // compute median peak width
    std::vector<DoubleReal> pw_vec;
    std::vector<Size> pw_idx_vec;

    for (std::multimap<DoubleReal, Size>::const_iterator c_it = histo_map.begin(); c_it != histo_map.end(); ++c_it)
    {
        pw_vec.push_back(c_it->first);
        pw_idx_vec.push_back(c_it->second);
    }

    DoubleReal pw_median(Math::median(pw_vec.begin(), pw_vec.end(), true));

    //    // compute median of absolute deviances (MAD)
    //    std::vector<DoubleReal> abs_devs;

    //    for (Size pw_i = 0; pw_i < pw_vec.size(); ++pw_i)
    //    {
    //        abs_devs.push_back(std::fabs(pw_vec[pw_i] - pw_median));
    //    }

    //    // Size abs_devs_size = abs_devs.size();
    //    DoubleReal pw_mad(Math::median(abs_devs.begin(), abs_devs.end(), false));
    // std::cout << "pw_median: " << pw_median << std::endl;
    chrom_fwhm_ = pw_median;


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

    mt.findLocalExtrema(win_size / 2, maxes, mins);

    // if only one maximum exists: finished!
    if (maxes.size() == 1)
    {


        if (!mt_snr_filtering_ || mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
        {
            // mt.updateSmoothedMaxRT();
            // mt.updateWeightedMeanRT();

            mt.updateSmoothedWeightedMeanRT();
            mt.estimateFWHM(true);
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


            if (tmp_mt.size() >= win_size / 2)
            {

                DoubleReal scantime(mt.getScanTime());

                MassTrace new_mt(tmp_mt, scantime);

                // copy smoothed ints
                new_mt.setSmoothedIntensities(smoothed_tmp);

                if (!mt_snr_filtering_ || mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
                {

                    // set label of subtrace
                    String tr_num;
                    std::stringstream read_in;
                    read_in << (min_idx + 1);
                    tr_num = "." + read_in.str();

                    new_mt.setLabel(mt.getLabel() + tr_num);
                    //new_mt.updateWeightedMeanRT();
                    //new_mt.updateSmoothedMaxRT();
                    new_mt.updateSmoothedWeightedMeanRT();
                    new_mt.updateWeightedMeanMZ();
                    new_mt.updateWeightedMZsd();
                    new_mt.estimateFWHM(true);

#ifdef _OPENMP
#pragma omp critical
#endif
                    single_mtraces.push_back(new_mt);
                }
            }
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

        if (tmp_mt.size() >= win_size / 2)
        {
            DoubleReal scantime(mt.getScanTime());

            MassTrace new_mt(tmp_mt, scantime);

            // copy smoothed ints
            new_mt.setSmoothedIntensities(smoothed_tmp);

            if (!mt_snr_filtering_ || mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
            {

                // set label of subtrace
                String tr_num;
                std::stringstream read_in;
                read_in << (mins.size() + 1);
                tr_num = "." + read_in.str();

                new_mt.setLabel(mt.getLabel() + tr_num);
                //new_mt.updateWeightedMeanRT();
                //new_mt.updateSmoothedMaxRT();
                new_mt.updateSmoothedWeightedMeanRT();
                new_mt.updateWeightedMeanMZ();
                new_mt.updateWeightedMZsd();
                new_mt.estimateFWHM(true);


#ifdef _OPENMP
#pragma omp critical
#endif
                single_mtraces.push_back(new_mt);
            }
        }
    }
    return;
}

void ElutionPeakDetection::updateMembers_()
{
    chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");
    chrom_peak_snr_ = (DoubleReal)param_.getValue("chrom_peak_snr");
    noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
    pw_filtering_ = param_.getValue("width_filtering").toBool();
    mt_snr_filtering_ = param_.getValue("masstrace_snr_filtering").toBool();
}

} //namespace OpenMS
