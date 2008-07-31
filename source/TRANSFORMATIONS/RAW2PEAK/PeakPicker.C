// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>

namespace OpenMS
{
	PeakPicker::PeakPicker()
    :DefaultParamHandler("PeakPicker")
  {
  	defaults_.setValue("thresholds:signal_to_noise",1.0,"Minimal signal to noise ratio for a peak to be picked.");
		defaults_.setMinFloat("thresholds:signal_to_noise",0.0);
		defaults_.setValue("thresholds:peak_bound",10.0,"Minimal peak intensity.");
		defaults_.setMinFloat("thresholds:peak_bound",0.0);
  	defaults_.setValue("thresholds:peak_bound_ms2_level",10.0,"Minimal peak intensity for MS/MS peaks.");
		defaults_.setMinFloat("thresholds:peak_bound_ms2_level",0.0);
  	defaults_.setValue("thresholds:fwhm_bound",0.2,"Minimal peak width");
		defaults_.setMinFloat("thresholds:fwhm_bound",0.0);
		
		defaultsToParam_();
  }
	
	void PeakPicker::updateMembers_()
	{
    signal_to_noise_ = (float)param_.getValue("thresholds:signal_to_noise");
		peak_bound_ = (float)param_.getValue("thresholds:peak_bound");
		peak_bound_ms2_level_ = (float)param_.getValue("thresholds:peak_bound_ms2_level");
    fwhm_bound_ = (float)param_.getValue("thresholds:fwhm_bound");
	}
} // namespace OpenMS
