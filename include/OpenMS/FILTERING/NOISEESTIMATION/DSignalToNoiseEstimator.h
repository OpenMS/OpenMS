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
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATOR_H
#define OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATOR_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

#include <iostream>
#include <vector>

namespace OpenMS
{
  /**
    @brief  Signal to noise estimator. This class represents the abstract base class of a signal to noise estimator.

    A signal to noise estimator should provide the signal to noise ratio of all raw data points
    in a given intervall [first_,last_).
  */
  template <Size D = 1 , typename PeakIterator = MSSpectrum<DRawDataPoint<1> >::const_iterator >
  class DSignalToNoiseEstimator
  {
  public:
    typedef DimensionDescription < LCMS_Tag > DimensionDescription;

    /// Constructor
    inline DSignalToNoiseEstimator()
        : first_(0),
        last_(0)
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }
    }

    /// Constructor
    inline DSignalToNoiseEstimator(const Param& parameters): first_(0), last_(0), param_(parameters)
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }
    }
    
    /// Copy constructor
    inline DSignalToNoiseEstimator(const DSignalToNoiseEstimator&  ne)
        : mz_dim_(ne.mz_dim_),
        rt_dim_(ne.rt_dim_),
        first_(ne.first_),
        last_(ne.last_),
        param_(ne.param_)
    {}
    
    /// Destructor
    virtual ~DSignalToNoiseEstimator()
    {}

    /// Assignment operator
    inline DSignalToNoiseEstimator& operator=(const DSignalToNoiseEstimator& ne)
    {
      mz_dim_=ne.mz_dim_;
      rt_dim_=ne.rt_dim_;
      first_=ne.first_;
      last_=ne.last_;
      param_=ne.param_;

      return *this;
    }

    /** @name Initialisation of the raw data intervall
        
        Set the start and endpoint of the raw data intervall, for which signal to noise ratios should be estimated
    */
    virtual void init(PeakIterator it_begin, PeakIterator it_end)
    {
      first_=it_begin;
      last_=it_end;
    }

    /// Non-mutable access to mz dimension
    inline const int& getMZdim() const { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline int& getMZdim() { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline void setMZdim(const int& mz_dim) { mz_dim_ = mz_dim; }

    /// Non-mutable access to the rt dimension
    inline const int getRTdim() const { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline int& getRTdim() { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline void setRTdim(const int& rt_dim) { rt_dim_ = rt_dim; }

    /// Non-mutable access to the first raw data point
    inline const PeakIterator& getFirstDataPoint() const { return first_; }
    /// Mutable access to the first raw data point
    inline PeakIterator& getFirstDataPoint() { return first_; }
    /// Mutable access to the first raw data point
    inline void setFirstDataPoint(const PeakIterator& first) { first_ = first; }

    /// Non-mutable access to the last raw data point
    inline const PeakIterator& getLastDataPoint() const { return last_; }
    /// Mutable access to the last raw data point
    inline PeakIterator& getLastDataPoint() { return last_; }
    /// Mutable access to the last raw data point
    inline void setLastDataPoint(const PeakIterator& last) { last_ = last; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) { param_ = param; }

    
    /// Signal To Noise Estimation
    virtual double getSignalToNoise(PeakIterator data_point) = 0;

  protected:
    /// m/z dimension
    int mz_dim_;
    /// retention time dimension
    int rt_dim_;
    /// points to the first raw data point in the intervall
    PeakIterator first_;
    /// points to the right position next to the last raw data point in the intervall
    PeakIterator last_;
    /// parameter object
    Param param_;
  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATOR_H
