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
    defaults_.setValue( "chrom_fwhm" , 0.0 , "Allows filtering of mass traces with peak width less than this threshold. Disabled by default (set to 0.0).", StringList::create("advanced"));
    defaults_.setValue("width_filtering", "true", "Enable filtering of unlikely peak widths.");
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

void ElutionPeakDetection::filterByPeakWidth(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& filt_mtraces)
{
    std::multimap<DoubleReal, Size> histo_map;

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        DoubleReal fwhm(mt_vec[i].estimateFWHM(true));

        histo_map.insert(std::make_pair(fwhm, i));
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

    for (Size i = 0; i < mt_vec.size(); ++i)
    {
        // drop any masstrace with lower pw than lower_pw_bound
        if (pw_vec[i] > lower_pw_bound && pw_vec[i] > chrom_fwhm_)
        {
            filt_mtraces.push_back(mt_vec[pw_idx_vec[i]]);
        }
    }

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


    Size win_size = mt.getFWHMScansNum();

    // if there is no previous FWHM estimation... do it now
    if (win_size == 0)
    {
        mt.estimateFWHM(false); // estimate FWHM
        win_size = mt.getFWHMScansNum();
    }

    lowess_params.setValue("window_size", win_size);
    lowess_smooth.setParameters(lowess_params);

    lowess_smooth.smoothData(rts, ints, smoothed_data);

    mt.setSmoothedIntensities(smoothed_data);

    std::vector<Size> maxes, mins;

    mt.findLocalExtrema(win_size/2, maxes, mins);

    // if only one maximum exists: finished!
    if (maxes.size() == 1)
    {
        single_mtraces.push_back(mt);
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

            MassTrace new_mt(tmp_mt);

            // copy smoothed ints
            new_mt.setSmoothedIntensities(smoothed_tmp);

            // set label of subtrace
            String tr_num;
            std::stringstream read_in;
            read_in << (min_idx + 1);
            tr_num = "." + read_in.str();

            new_mt.setLabel(mt.getLabel() + tr_num);
            new_mt.updateWeightedMeanRT();
            new_mt.updateWeightedMeanMZ();

            single_mtraces.push_back(new_mt);
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

        MassTrace new_mt(tmp_mt);

        // copy smoothed ints
        new_mt.setSmoothedIntensities(smoothed_tmp);

        // set label of subtrace
        String tr_num;
        std::stringstream read_in;
        read_in << (mins.size() + 1);
        tr_num = "." + read_in.str();

        new_mt.setLabel(mt.getLabel() + tr_num);
        new_mt.updateWeightedMeanRT();
        new_mt.updateWeightedMeanMZ();


        single_mtraces.push_back(new_mt);
    }
    return ;
}


void ElutionPeakDetection::updateMembers_()
{
    chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");
    pw_filtering_ = param_.getValue("width_filtering").toBool();
}










} //namespace OpenMS
