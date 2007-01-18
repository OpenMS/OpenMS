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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H

#include <OpenMS/FORMAT/Param.h>

#include <iostream>
#include <vector>
#include <math.h>


namespace OpenMS
{
  class String;

  /**
    @defgroup Transformations Transformations

    @brief Classes for the transformation of ms data.

    This module contains all classes that are involved in a data-reduction
    method (e.g. the transformation of raw data into peak data).
  */

  /**
     @brief This class is the base class for every peak picker.
     
     @ingroup PeakPicking
  */
  class PeakPicker
  {

  public:
    /// Constructor
    PeakPicker();

    /// Copy constructor
    PeakPicker(const PeakPicker& pp);

    /// Destructor
    virtual ~PeakPicker()
    {	
    }

    /// Assignment operator
    PeakPicker& operator=(const PeakPicker& pp)
	  {
	    // take care of self assignments
	    if (this == &pp)
	    {
	      return *this;
	    }
			param_ = pp.param_;
	    peak_bound_= pp.peak_bound_;
	    peak_bound_ms2_level_= pp.peak_bound_ms2_level_;
	    signal_to_noise_= pp.signal_to_noise_;
	    fwhm_bound_ = pp.fwhm_bound_;
	
	    return *this;
	  }

    /// Non-mutable access to the threshold of the height
    inline const float& getPeakBound() const 
    { 
    	return peak_bound_; 
    }
    /// Mutable access to the threshold of the height
    inline void setPeakBound(const float& peak_bound) 
    { 
    	peak_bound_ = peak_bound;
    	param_.setValue("thresholds:peak_bound",peak_bound);
    }

    /// Non-mutable access to the threshold of the peak height in the MS 2 level
    inline const float& getPeakBoundMs2Level() const 
    {
    	return peak_bound_ms2_level_;
    }
    /// Mutable access to the threshold of the peak height in the MS 2 level
    inline void setPeakBoundMs2Level(const float& peak_bound_ms2_level) 
    { 
    	peak_bound_ms2_level_ = peak_bound_ms2_level;
    	param_.setValue("thresholds:peak_bound_ms2_level",peak_bound_ms2_level); 
    }

    /// Non-mutable access to the signal to noise threshold
    inline const float& getSignalToNoiseLevel() const 
    { 
    	return signal_to_noise_; 
    }
    /// Mutable access to the signal to noise threshold
    inline void setSignalToNoiseLevel(const float& signal_to_noise) 
    { 
    	signal_to_noise_ = signal_to_noise;
    	param_.setValue("thresholds:signal_to_noise",signal_to_noise);
    }

    /// Non-mutable access to the fwhm threshold
    inline const float& getFwhmBound() const 
    { 
    	return fwhm_bound_; 
    }
    /// Mutable access to the fwhm threshold
    inline void setFwhmBound(const float& fwhm) 
    {
    	fwhm_bound_ = fwhm; 
    	param_.setValue("thresholds:fwhm_bound",fwhm);
    }

		/**
			@brief Checks the give parameters against defaults_ and copies them to member variables.
			
			Call this method in the constructor of derived classes after the defaults_ are set!
			
			@note Overwrite this method in derived classes! You need to copy the subclass-specific parameters to members as well! 
		*/
		virtual void setParam(Param param);
		
		/// Returns the set parameters
		const Param& getParam() const;
		
  protected:
		/// passed parameters
		Param param_;

		/// default parameters
		Param defaults_;
		
    /// Threshold for the peak height in the MS 1 level
    float peak_bound_;

    /// Threshold for the peak height in the MS 2 level
    float peak_bound_ms2_level_;

    /// Signal to noise threshold
    float signal_to_noise_;

    /// The minimal full width at half maximum
    float fwhm_bound_;

  };

}// namespace OpenMS

#endif
