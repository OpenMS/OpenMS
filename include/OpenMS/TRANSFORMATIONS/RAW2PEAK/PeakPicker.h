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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iostream>
#include <vector>
#include <math.h>


namespace OpenMS
{
  /**
		@brief This class is the base class for every peak picker.
		 
		@ref PeakPicker_Parameters are explained on a separate page.
		
		@ingroup PeakPicking
  */
  class PeakPicker
  	: public DefaultParamHandler, 
  		public ProgressLogger
  {

  public:
    /// Constructor
    PeakPicker();

    /// Destructor
    virtual ~PeakPicker()
    {	
    }
		
  protected:
    /// Threshold for the peak height in the MS 1 level
    float peak_bound_;

    /// Threshold for the peak height in the MS 2 level
    float peak_bound_ms2_level_;

    /// Signal to noise threshold
    float signal_to_noise_;

    /// The minimal full width at half maximum
    float fwhm_bound_;

		void updateMembers_();
		
  };

}// namespace OpenMS

#endif
