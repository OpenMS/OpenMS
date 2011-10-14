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


#include "OpenMS/FILTERING/TRANSFORMERS/TICResampling.h"


namespace OpenMS
{
    TICResampling::TICResampling()
        : DefaultParamHandler("TICResampling")
    {
        defaults_.setValue("scan_diff", 0.5, "RT difference between consecutive scans (in seconds)");

        defaultsToParam_();
    }

    TICResampling::~TICResampling()
    {
    }

    void TICResampling::run(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output)
    {



        for (Size scan_idx = 1; scan_idx < input.size(); ++scan_idx)
        {
            MSSpectrum<Peak1D> tmp_scan(input[scan_idx - 1]);
            DoubleReal base_rt(tmp_scan.getRT());

            DoubleReal scan_rt_diff = std::fabs(input[scan_idx].getRT() - base_rt);

            Size num_duplicate(std::floor(scan_rt_diff/scan_diff_));

            if (num_duplicate == 0)
            {
                continue;
            }

            // DoubleReal rt_step(scan_rt_diff/num_duplicate);

            std::cout << "params_ " << num_duplicate << std::endl; // << " " << rt_step << std::endl;

            for (Size dupl_idx = 0; dupl_idx < num_duplicate; ++dupl_idx)
            {
                tmp_scan.setRT(base_rt + dupl_idx*scan_diff_);
                output.push_back(tmp_scan);
            }

        }

        std::cout << "new size: " << output.size() << std::endl;

        return ;

    }



    void TICResampling::updateMembers_()
    {
        scan_diff_ = (DoubleReal)param_.getValue("scan_diff");
    }

} //namespace OpenMS
