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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKER_H

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>
#include <vector>
#include <math.h>


namespace OpenMS
{
  /**
    @defgroup Transformations Transformations

    @brief Classes for the transformation of ms data.

    This module contains all classes that are involved in a data-reduction
    method (e.g. the transformation of raw data into peak data).
  */

  /**
     @brief This class is the base class for every peak picker.

     @ingroup PeakPicking
     
     @todo Make it work on all classes derived from DRawDataPoint (Eva)

  */
  class PeakPicker
  {

  public:
   /// Constructor
    PeakPicker()
        : peak_bound_(200),
        peak_bound_ms2_level_(50),
        signal_to_noise_(3),
        fwhm_bound_(0.2) {}
        
    /// Constructor given the name of a param file
    PeakPicker(const String& param_filename);
    
    /// Constructor given a param object
    PeakPicker(const Param& parameters);
    
    /// Copy constructor
    PeakPicker(const PeakPicker& pp);
    
    /// Destructor
    virtual ~PeakPicker()
    {   }
   
     /// Assignment operator
    PeakPicker& operator=(const PeakPicker& pp);
    
    
    /// Non-mutable access to the threshold of the height
    inline const float& getPeakBound() const { return peak_bound_; }
    /// Mutable access to the threshold of the height
    inline float& getPeakBound()  { return peak_bound_; }
    /// Mutable access to the threshold of the height
    virtual void setPeakBound(const float& peak_bound) { peak_bound_ = peak_bound;  }

    /// Non-mutable access to the threshold of the peak height in the MS 2 level
    inline const float& getPeakBoundMs2Level() const { return peak_bound_ms2_level_; }
    /// Mutable access to the threshold of the peak height in the MS 2 level
    inline float& getPeakBoundMs2Level() { return peak_bound_ms2_level_; }
    /// Mutable access to the threshold of the peak height in the MS 2 level
    inline void setPeakBoundMs2Level(const float& peak_bound_ms2_level) { peak_bound_ms2_level_ = peak_bound_ms2_level; }

    /// Non-mutable access to the signal to noise threshold
    inline const float& getSignalToNoiseLevel() const { return signal_to_noise_; }
    /// Mutable access to the signal to noise threshold
    inline float& getSignalToNoiseLevel() { return signal_to_noise_; }
    /// Mutable access to the signal to noise threshold
    inline void setSignalToNoiseLevel(const float& signal_to_noise) { signal_to_noise_ = signal_to_noise; }
    
    /// Non-mutable access to the fwhm threshold
    inline const float& getFwhmBound() const { return fwhm_bound_; }
    /// Mutable access to the fwhm threshold
    inline float& getFwhmBound() { return fwhm_bound_; }
    /// Mutable access to the fwhm threshold
    inline void setFwhmBound(const float& fwhm) { fwhm_bound_ = fwhm; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline Param& getParam() { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) { param_ = param; }
    //@}

  protected:
    /// Parameter object
    Param param_;
   
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
