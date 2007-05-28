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
#include <OpenMS/CONCEPT/ProgressLogger.h>

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
  
template < typename Container = MSSpectrum< > >
class SignalToNoiseEstimator: public DefaultParamHandler, public ProgressLogger
  {
  public:

    /** @name Type definitions
     */
    //@{
    typedef typename Container::const_iterator PeakIterator;
    typedef typename PeakIterator::value_type PeakType;
    
    
    //@}
    
    /// Constructor
    inline SignalToNoiseEstimator()
        : DefaultParamHandler("SignalToNoiseEstimator"),
          ProgressLogger(),
          first_(0),
          last_(0),
          is_result_valid_(false)
    {
    }

    /// Copy constructor
    inline SignalToNoiseEstimator(const SignalToNoiseEstimator&  source)
        : DefaultParamHandler(source),
          ProgressLogger(source),
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
      ProgressLogger::operator=(source);
      stn_estimates_ = source.stn_estimates_;
      first_ = source.first_;
      last_  = source.last_;
      return *this;
    }
    
    /// Destructor
    virtual ~SignalToNoiseEstimator()
    {}

    /** @name Accessors
     */

    //@{
    ///
    /// Non-mutable access to the first raw data point
    inline const PeakIterator& getFirstDataPoint() const { return first_; }
    /// Mutable access to the first raw data point
    inline void setFirstDataPoint(const PeakIterator& first) { is_result_valid_=false; first_ = first; }

    /// Non-mutable access to the last raw data point
    inline const PeakIterator& getLastDataPoint() const { return last_; }
    /// Mutable access to the last raw data point
    inline void setLastDataPoint(const PeakIterator& last) { is_result_valid_=false; last_ = last; }
    

    /// Non-mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
    virtual DoubleReal getMaxIntensity() const = 0;
    /// Mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
    virtual void setMaxIntensity(DoubleReal max_intensity) = 0;
//

    /// Non-Mutable access to the AutoMaxStdevFactor-Param, which holds a factor for stddev (only used if autoMode=1)
    virtual DoubleReal getAutoMaxStdevFactor() const = 0;
    /// Mutable access to the AutoMaxStdevFactor-Param, which holds a factor for stddev (only used if autoMode=1)
    virtual void setAutoMaxStdevFactor(DoubleReal value) = 0;
//      

    /// get the AutoMaxPercentile-Param, which holds a percentile (only used if autoMode=2)
    virtual DoubleReal getAutoMaxPercentile() const = 0;
    /// Mutable access to the AutoMaxPercentile-Param, which holds a percentile (only used if autoMode=2)
    virtual void setAutoMaxPercentile(DoubleReal value) = 0;
//

    /// @brief -1 will disable it. 0 is default. 1 is alternative method
    /// Non-mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
    virtual inline Int getAutoMode() const = 0;
    /// @brief -1 will disable it. 0 is default. 1 is alternative method
    /// Mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
    virtual void setAutoMode(Int auto_mode) = 0;
//

    /// Non-mutable access to the window length (in Thomson)
    virtual DoubleReal getWinLen() const = 0;
    /// Mutable access to the window length (in Thomson)
    virtual void setWinLen(DoubleReal win_len) = 0;

//
    /// Non-mutable access to the number of bins used for the histogram (the more bins, the better the approximation, but longer runtime)
    virtual Int getBinCount() const = 0;
    /// Mutable access to the number of bins used for the histogram
    virtual void setBinCount(Int bin_count) = 0;

//
    /// Non-mutable access to the minimum required elements in a window, to be evaluated.
    virtual Int getMinReqElements() const = 0;
    /// Mutable access to the minimum required elements in a window, to be evaluated.
    virtual void setMinReqElements(Int min_required_elements) = 0;

//
    /// Non-mutable access to the noise value that is used if a window contains not enough elements
    virtual DoubleReal getNoiseForEmtpyWindow() const = 0;
    /// Mutable access to the noise value that is used if a window contains not enough elements
    virtual void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window) = 0;
    //@}    
    
        
     
    /// Set the start and endpoint of the raw data intervall, for which signal to noise ratios will be estimated immediately
    virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)
    {
      first_=it_begin;
      last_=it_end;
      computeSTN_(first_, last_);
      is_result_valid_ = true;
    }
      
          
    /// Set the start and endpoint of the raw data intervall, for which signal to noise ratios will be estimated immediately
    virtual void init(const Container& c)
    {
      init(c.begin(), c.end() ); 
    }
    
    /// Return to signal/noise estimate for data point @p data_point
    /// @note the first query to this function will take longer, as
    ///       all SignalToNoise values are calculated
    /// @note you will get a warning to stderr if more than 20% of the
    ///       noise estimates used sparse windows
    virtual double getSignalToNoise(const PeakIterator& data_point)
    {
      if (!is_result_valid_)
      { 
        // recompute ...
        init(first_, last_);
      }

      return stn_estimates_[*data_point];
    }

  protected:

    virtual void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) throw(Exception::InvalidValue) = 0;
        


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
    inline GaussianEstimate estimate_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) const
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
    std::map< PeakType, double, typename PeakType::PositionLess > stn_estimates_;
  
    /// points to the first raw data point in the interval
    PeakIterator first_;
    /// points to the right position next to the last raw data point in the interval
    PeakIterator last_;
    /// flag: set to true if SignalToNoise estimates are calculated and none of the params were changed. otherwise false.
    mutable bool is_result_valid_;    
  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATOR_H
