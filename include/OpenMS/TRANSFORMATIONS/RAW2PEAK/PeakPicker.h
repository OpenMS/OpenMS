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
// $Maintainer: Erhan Kenar$
// --------------------------------------------------------------------------
//
#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
		@brief PeakPicker Factory Class
	
		@htmlinclude OpenMS_PeakPicker.parameters
	  
		@ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPicker
		: public FactoryProduct,
			public ProgressLogger
  {
	 public:

		/// Constructor
    PeakPicker();

    /// Destructor
    virtual ~PeakPicker();

    /** 
			@brief Applies the peak picking algorithm to a single spectrum.
	        
			Picks the peaks in the input spectrum and writes the resulting peaks to the output container.
	
			The ms_level should be one if the spectrum is a normal mass spectrum, or two if it is a tandem mass spectrum.
    */
    virtual void pick(const MSSpectrum<>& /*input*/, MSSpectrum<>& /*output*/);
		
    /** 
				@brief Picks the peaks in an MSExperiment.
			
				Picks the peaks successively in every scan in the spectrum range. The detected peaks are stored in the output MSExperiment. 
    */
    virtual void pickExperiment(const MSExperiment<>& /*input*/, MSExperiment<>& /*output*/);

		
		static void registerChildren();
		
 
	private:
		PeakPicker(const PeakPicker& );
		PeakPicker& operator=(const PeakPicker& );

		void updateMembers_();		
  }; // end PeakPicker


}// namespace OpenMS

#endif
