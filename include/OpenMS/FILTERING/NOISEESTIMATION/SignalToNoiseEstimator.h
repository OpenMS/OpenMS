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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATOR_H
#define OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATOR_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

namespace OpenMS
{
  /**
    @brief  Signal to noise estimator. This class represents the abstract base class of a signal to noise estimator.

    A signal to noise estimator should provide the signal to noise ratio of all raw data points
    in a given intervall [first_,last_).
  
  */
class SignalToNoiseEstimator: public DefaultParamHandler
  {
  public:

    /** @name Type definitions
     */
    //@{
    typedef MSSpectrum< >::const_iterator PeakIterator;
    typedef PeakIterator::value_type PeakType;

    //@}
    
    /// Constructor
    inline SignalToNoiseEstimator()
        : DefaultParamHandler("SignalToNoiseEstimator"),
          first_(0),
          last_(0),
          is_result_valid_(false)
    {
    }

    /// Copy constructor
    inline SignalToNoiseEstimator(const SignalToNoiseEstimator&  source)
        : DefaultParamHandler(source),
        stn_estimates_( source.stn_estimates_),
        first_(source.first_),
        last_(source.last_),
        is_result_valid_(source.is_result_valid_)
    {}

    /// Assignment operator
    inline SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator& source)
    {
      if(&source == this) return *this; 
      
      DefaultParamHandler::operator=(source);
      stn_estimates_ = source.stn_estimates_;
      first_ = source.first_;
      last_  = source.last_;
      return *this;
    }
    
    /// Destructor
    virtual ~SignalToNoiseEstimator()
    {}

    /** @name Initialisation of the raw data intervall
        
        Set the start and endpoint of the raw data intervall, for which signal to noise ratios should be estimated
    */
    virtual void init(PeakIterator it_begin, PeakIterator it_end)
    {
      first_=it_begin;
      last_=it_end;
      is_result_valid_=false;
    }

    /// Non-mutable access to the first raw data point
    inline const PeakIterator& getFirstDataPoint() const { return first_; }
    /// Mutable access to the first raw data point
    inline void setFirstDataPoint(const PeakIterator& first) { is_result_valid_=false; first_ = first; }

    /// Non-mutable access to the last raw data point
    inline const PeakIterator& getLastDataPoint() const { return last_; }
    /// Mutable access to the last raw data point
    inline void setLastDataPoint(const PeakIterator& last) { is_result_valid_=false; last_ = last; }

    
    /// Signal To Noise Estimation
    virtual double getSignalToNoise(PeakIterator data_point) = 0;

  protected:

    /** 
      @brief protected struct to store parameters my, sigma for a gaussian distribution
    
      Accessors are : mean and variance      
    */ 
    struct GaussianEstimate
    {
      double mean;   ///mean of estimated gaussian
      double variance; ///variance of estimated gaussian
    };  


    /// calculate mean & stdev of intensities of a DPeakArray
    inline GaussianEstimate estimate(const PeakIterator& scan_first_, const PeakIterator& scan_last_) const
    {
      int size = 0;
      // add up
      double v = 0;
      double m = 0;
      PeakIterator run = scan_first_;
      while (run != scan_last_)
      {
        m += (*run).getIntensity();
        ++size;
        ++run;
      }
      //average
      m = m/size;

      //determine variance
      run = scan_first_;
      while (run != scan_last_)
      {
        v += std::pow(m - (*run).getIntensity(), 2);
        ++run;
      }
      v = v / ((double)size); // divide by n

      GaussianEstimate value = {m, v};
      return value;
    }

    //MEMBERS:

    /// stores the noise estimate for each peak
    std::map< PeakType, double, PeakType::PositionLess > stn_estimates_;
  
    /// points to the first raw data point in the interval
    PeakIterator first_;
    /// points to the right position next to the last raw data point in the interval
    PeakIterator last_;
    /// flag: set to true if SignalToNoise estimates are calculated and none of the params were changed. otherwise false.
    mutable bool is_result_valid_;    
  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATOR_H
