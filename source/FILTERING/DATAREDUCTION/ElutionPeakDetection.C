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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#include "OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h"
#include "OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h"
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <sstream>
#include <numeric>
#include <algorithm>

namespace OpenMS
{
ElutionPeakDetection::ElutionPeakDetection()
    : DefaultParamHandler("ElutionPeakDetection"), ProgressLogger()
{
    defaults_.setValue("chrom_fwhm", 5.0 , "Expected full-width-at-half-maximum of chromatographic peaks.");
    defaults_.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    defaults_.setValue("noise_threshold_int" , 10.0 , "Intensity threshold below which peaks are regarded as noise.");

    defaults_.setValue("width_filtering", "true", "Enable filtering of unlikely peak widths (5 and 95% quantiles of peak width distribution).");
    defaults_.setValidStrings("width_filtering", StringList::create(("true,false")));

    defaultsToParam_();

    this->setLogType(CMD);
}

ElutionPeakDetection::~ElutionPeakDetection()
{
}

void ElutionPeakDetection::detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces)
{
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    detectElutionPeaks_(mt, single_mtraces);
    return ;
}

void ElutionPeakDetection::detectPeaks(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& single_mtraces)
{
    // make sure that single_mtraces is empty
    single_mtraces.clear();

    this->startProgress(0, mt_vec.size(), "elution peak detection");

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        this->setProgress(i);
        detectElutionPeaks_(mt_vec[i], single_mtraces);
    }

    this->endProgress();

    return ;
}


void ElutionPeakDetection::estimatePeakWidth(std::vector<MassTrace>& mt_vec)
{
    std::multimap<DoubleReal, Size> histo_map;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        DoubleReal fwhm(mt_vec[i].estimateFWHM(false));

        if (fwhm > 0.0) {
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


    return ;
}


void ElutionPeakDetection::filterByPeakWidth(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& filt_mtraces)
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

    return ;
}

void ElutionPeakDetection::detectElutionPeaks_(MassTrace& mt, std::vector<MassTrace>& single_mtraces)
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
    Size win_size = std::ceil(chrom_fwhm_/scan_time_);

    // std::cout << "win_size elution: " << win_size << std::endl;

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

    mt.findLocalExtrema(win_size/2, maxes, mins);

    // if only one maximum exists: finished!
    if (maxes.size() == 1)
    {
        if (mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
        {
            mt.updateSmoothedMaxRT();
            single_mtraces.push_back(mt);
        }
    }
    else if (maxes.empty())
    {
        return ;
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


            if (tmp_mt.size() >= win_size/2) {

                DoubleReal scantime(mt.getScanTime());

                MassTrace new_mt(tmp_mt, scantime);

                // copy smoothed ints
                new_mt.setSmoothedIntensities(smoothed_tmp);

                if (mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
                {

                    // set label of subtrace
                    String tr_num;
                    std::stringstream read_in;
                    read_in << (min_idx + 1);
                    tr_num = "." + read_in.str();

                    new_mt.setLabel(mt.getLabel() + tr_num);
                    //new_mt.updateWeightedMeanRT();
                    new_mt.updateSmoothedMaxRT();
                    new_mt.updateWeightedMeanMZ();


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

        if (tmp_mt.size() >= win_size/2) {
            DoubleReal scantime(mt.getScanTime());

            MassTrace new_mt(tmp_mt, scantime);

            // copy smoothed ints
            new_mt.setSmoothedIntensities(smoothed_tmp);

            if (mt.computeSNR(true, noise_threshold_int_) > chrom_peak_snr_)
            {

                // set label of subtrace
                String tr_num;
                std::stringstream read_in;
                read_in << (mins.size() + 1);
                tr_num = "." + read_in.str();

                new_mt.setLabel(mt.getLabel() + tr_num);
                //new_mt.updateWeightedMeanRT();
                new_mt.updateSmoothedMaxRT();
                new_mt.updateWeightedMeanMZ();



                single_mtraces.push_back(new_mt);
            }
        }
    }
    return ;
}


void ElutionPeakDetection::updateMembers_()
{
    chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");
    chrom_peak_snr_ = (DoubleReal)param_.getValue("chrom_peak_snr");
    noise_threshold_int_ = (DoubleReal)param_.getValue("noise_threshold_int");
    pw_filtering_ = param_.getValue("width_filtering").toBool();
}










} //namespace OpenMS
