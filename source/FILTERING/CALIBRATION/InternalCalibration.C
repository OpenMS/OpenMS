// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>


namespace OpenMS
{

	InternalCalibration::InternalCalibration()
		:DefaultParamHandler("InternalCalibration")
	{
		defaults_.setValue("peak_bound",1400);

		defaultsToParam_();
		updateMembers_();
	}
	
  InternalCalibration::InternalCalibration(InternalCalibration& obj)
		:DefaultParamHandler(obj),exp_peaks_(obj.exp_peaks_),monoiso_peaks_(obj.monoiso_peaks_)
  {
		updateMembers_();
  }
  
  InternalCalibration& InternalCalibration::operator=(const InternalCalibration& obj)
  {
		// take care of self assignments
    if (this == &obj)		return *this;
		DefaultParamHandler::operator=(obj);
		updateMembers_();
		exp_peaks_=obj.exp_peaks_;
		monoiso_peaks_=obj.monoiso_peaks_;
    return *this;

  }

	void InternalCalibration::getMonoisotopicPeaks_()
	{
		
		MSExperiment<PickedPeak1D >::iterator spec_iter = exp_peaks_.begin();
		MSSpectrum<PickedPeak1D >::iterator peak_iter, help_iter;
#ifdef DEBUG_CALIBRATION
		spec_iter = exp_peaks_.begin();
		std::cout << "\n\nbefore---------\n\n";
		// iterate through all spectra
		for(;spec_iter != exp_peaks_.end();++spec_iter)
			{
				peak_iter = spec_iter->begin();
				// go through current scan
				for(;peak_iter != spec_iter->end();++peak_iter)
					{
						std::cout << peak_iter->getMZ() << std::endl;
					}
			}

#endif
		spec_iter = exp_peaks_.begin();
		// iterate through all spectra
		for(;spec_iter != exp_peaks_.end();++spec_iter)
			{
				peak_iter = spec_iter->begin();
				help_iter = peak_iter;
				std::vector<unsigned int> vec;
				// go through current scan
				while(peak_iter < spec_iter->end())
					{
						while(peak_iter+1 < spec_iter->end() && ( (peak_iter+1)->getMZ() - peak_iter->getMZ() < 1.2) )
							{
								++peak_iter;
							}
						
						vec.push_back(distance(spec_iter->begin(),help_iter));
					
						help_iter = peak_iter+1;
						++peak_iter;
						
					}
					monoiso_peaks_.push_back(vec);

			}
		
#ifdef DEBUG_CALIBRATION

		
		std::cout << "\n\nafter---------\n\n";

		for(unsigned int i=0;i<monoiso_peaks_.size();++i)
			{
				for(unsigned int j=0;j<monoiso_peaks_[i].size();++j)
					{
						std::cout << ( (exp_peaks_.begin() + +i)->begin() + (monoiso_peaks_[i])[j])->getMZ() << std::endl;
					}
				std::cout << "--------------\n";
						
			}
		std::cout << "--------------\n\n\n";
#endif
	}

	void InternalCalibration::updateMembers_()
	{
		peak_bound_ = (double)param_.getValue("peak_bound");
	}
}
