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
// $Maintainer: $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#include "OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h"
#include "OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h"

#include <sstream>
#include <numeric>
#include <algorithm>

namespace OpenMS
{
    ElutionPeakDetection::ElutionPeakDetection()
        : DefaultParamHandler("ElutionPeakDetection"), ProgressLogger()
    {
        defaults_.setValue("window_size", 5, "Number of neighbouring values in each direction that have to be exceeded to be considered a local maximum");

        defaultsToParam_();

        this->setLogType(CMD);
    }

    ElutionPeakDetection::~ElutionPeakDetection()
    {
    }

    void ElutionPeakDetection::detectPeaks(MassTrace& mt, std::vector<MassTrace>& single_mtraces)
    {
        detectElutionPeaks_(mt, single_mtraces);
        return ;
    }

    void ElutionPeakDetection::detectPeaks(std::vector<MassTrace>& mt_vec, std::vector<MassTrace>& single_mtraces)
    {
        // std::cout << "total mtraces: " << mt_vec.size() << std::endl;
        this->startProgress(0, mt_vec.size(), "elution peak detection");
        for (Size i = 0; i < mt_vec.size(); ++i)
        {
            this->setProgress(i);
            detectElutionPeaks_(mt_vec[i], single_mtraces);
            // std::cout << "trace: " << mt_vec[i].getLabel() << " finished." << std::endl;
        }
        this->endProgress();

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
        lowess_params.setValue("window_size", window_size_);
        lowess_smooth.setParameters(lowess_params);

        lowess_smooth.smoothData(rts, ints, smoothed_data);

        mt.setSmoothedIntensities(smoothed_data);

        std::vector<Size> maxes, mins;

        // std::cout << "winsize: " << window_size_/2 << std::endl;

        mt.findLocalExtrema(window_size_/2, maxes, mins);

        // only one maximum: finished!
        if (maxes.size() == 1)
        {
            single_mtraces.push_back(mt);
        }
        else if (maxes.size() == 0)
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
                MassTrace tmp_mt;
                std::vector<DoubleReal> smoothed_tmp;

                while (last_idx <= mins[min_idx])
                {
                    tmp_mt.appendPeak(*cp_it);
                    smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
                    ++cp_it;
                    ++last_idx;
                }

                // copy smoothed ints
                tmp_mt.setSmoothedIntensities(smoothed_tmp);

                // set label of subtrace
                String tr_num;
                std::stringstream read_in;
                read_in << min_idx;
                tr_num = "_" + read_in.str();

                tmp_mt.setLabel(mt.getLabel() + tr_num);

                single_mtraces.push_back(tmp_mt);
            }

            // don't forget the trailing trace
            MassTrace tmp_mt;

            std::vector<DoubleReal> smoothed_tmp;

            while (last_idx < mt.getSize())
            {
                tmp_mt.appendPeak(*cp_it);
                smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
                ++cp_it;
                ++last_idx;
            }

            // copy smoothed ints
            tmp_mt.setSmoothedIntensities(smoothed_tmp);

            // set label of subtrace
            String tr_num;
            std::stringstream read_in;
            read_in << mins.size();
            tr_num = "_" + read_in.str();

            tmp_mt.setLabel(mt.getLabel() + tr_num);

            single_mtraces.push_back(tmp_mt);
        }
        return ;
    }


    void ElutionPeakDetection::updateMembers_()
    {
        window_size_ = (Size)param_.getValue("window_size");
    }










} //namespace OpenMS
