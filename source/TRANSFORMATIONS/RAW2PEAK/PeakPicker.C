// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>

namespace OpenMS
{
  PeakPicker::PeakPicker(const String& filename)
  {
    param_.load(filename);

    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
    DataValue dv;
    dv = param_.getValue("thresholds:signal_to_noise");
    if (dv.isEmpty() || dv.toString() == "") signal_to_noise_ = 5;
    else signal_to_noise_ = (float)dv;

    dv = param_.getValue("thresholds:peak_bound");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ = 200;
    else peak_bound_ = (float)dv;

    dv = param_.getValue("thresholds:peak_bound_ms2_level");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ms2_level_ = 30;
    else peak_bound_ms2_level_ = (float)dv;
    	
    dv = param_.getValue("thresholds:fwhm_bound");
    if (dv.isEmpty() || dv.toString() == "") fwhm_bound_ = 0.2;
    else fwhm_bound_ = (float)dv;
  }

  PeakPicker::PeakPicker(const Param& parameters)
  {
    param_ = parameters;

    // if a peak picking parameter is missed in the param object the value should be substituted by a dv value
    DataValue dv;
    dv = param_.getValue("thresholds:signal_to_noise");
    if (dv.isEmpty() || dv.toString() == "") signal_to_noise_ = 3;
    else signal_to_noise_ = (float)dv;

    dv = param_.getValue("thresholds:peak_bound");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ = 200;
    else peak_bound_ = (float)dv;

    dv = param_.getValue("thresholds:peak_bound_ms2_level");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ms2_level_ = 30;
    else peak_bound_ms2_level_ = (float)dv;
    	
    dv = param_.getValue("thresholds:fwhm_bound");
    if (dv.isEmpty() || dv.toString() == "") fwhm_bound_ = 0.2;
    else fwhm_bound_ = (float)dv;
  }

  PeakPicker::PeakPicker(const PeakPicker& pp)
      : param_(pp.param_), 
      peak_bound_(pp.peak_bound_),
      peak_bound_ms2_level_(pp.peak_bound_ms2_level_),
      signal_to_noise_(pp.signal_to_noise_),
      fwhm_bound_(pp.fwhm_bound_)     
{}
} // namespace OpenMS
