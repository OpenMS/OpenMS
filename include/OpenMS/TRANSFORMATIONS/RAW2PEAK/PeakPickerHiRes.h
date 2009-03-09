// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHPP_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHPP_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
		@brief This class makes nothing.
	
		@htmlinclude OpenMS_PeakPickerCWT.parameters
	  
		@ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerHiRes
		: public DefaultParamHandler,
			public ProgressLogger
  {
	 public:
    /// Constructor
    PeakPickerHiRes();

    /// Destructor
    virtual ~PeakPickerHiRes();

    /** 
			@brief pick
    */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType>& /*input*/, MSSpectrum<PeakType>& /*output*/)
    {
    }

		
    /** 
			@brief 
    */
    template <typename PeakType>
    void pickExperiment(const MSExperiment<PeakType>& /*input*/, MSExperiment<PeakType>& /*output*/)
    {
    }
		
	 protected:
	 	//docu in base class
		void updateMembers_();

  }; // end PeakPickerHiRes


}// namespace OpenMS

#endif
