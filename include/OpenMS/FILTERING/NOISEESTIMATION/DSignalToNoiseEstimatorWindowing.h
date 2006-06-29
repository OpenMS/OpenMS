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
// $Id: DSignalToNoiseEstimatorWindowing.h,v 1.17 2006/05/29 15:50:58 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORWINDOWING_H
#define OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORWINDOWING_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimator.h>

#include <iostream>
#include <vector>
#include <math.h>

#include <fstream>

namespace OpenMS
{
  /**
    @brief Signal to noise estimator, estimates the signal/noise ratio of each data point in a scan.

    This class provides the estimation of the signal to noise ratios in a given raw data
    points intervall using the method of Roegnvaldsson et al. described in "Modular, Scriptable, and Automated
    Analysis Tools for High-Throughput Peptide Mass Fingerprinting".
	
    NOTE: This algorithm works per scan ONLY i.e. you have to call init() with an iterator range
    for each scan, and not for the whole map.

	???? (ost) This class returns occasionally negative (!) s/n ratios ?? I could not figure out why.
	Please use it carefully.
	
    @ingroup PeakPickingCWT
  */
  template <Size D, typename  ContainerType = DPeakArrayNonPolymorphic<D,DRawDataPoint<D> > >
  class DSignalToNoiseEstimatorWindowing : public DSignalToNoiseEstimator<D, ContainerType>
  {

  public:
    /** @name Type definitions
     */
    //@{
    ///
    typedef typename ContainerType::const_iterator ConstIterator;
    //@}

    using DSignalToNoiseEstimator<D, ContainerType>::mz_dim_;
    using DSignalToNoiseEstimator<D, ContainerType>::rt_dim_;
    using DSignalToNoiseEstimator<D, ContainerType>::param_;
    using DSignalToNoiseEstimator<D, ContainerType>::first_;
    using DSignalToNoiseEstimator<D, ContainerType>::last_;


    /** @name Constructors and Destructor
     */
    //@{
    inline DSignalToNoiseEstimatorWindowing()
        : DSignalToNoiseEstimator<D,ContainerType>(),
        bucket_size_(10),
        window_size_(700)
    {}
    ///
    inline DSignalToNoiseEstimatorWindowing(const Param& parameters)
        : DSignalToNoiseEstimator<D,ContainerType>(parameters)
    {
      // if a parameter is missing in the param object, the value is substituted by a default value
      DataValue dv = param_.getValue("SignalToNoiseEstimationParameter:Bucket");
      if (dv.isEmpty() || dv.toString() == "") bucket_size_ = 10;
      else bucket_size_ = (int)dv;

      dv = param_.getValue("SignalToNoiseEstimationParameter:Window");
      if (dv.isEmpty() || dv.toString() == "") window_size_ = 700;
      else window_size_ = (int)dv;
    }

    ///
    inline DSignalToNoiseEstimatorWindowing(const DSignalToNoiseEstimatorWindowing&  ne)
        : DSignalToNoiseEstimator<D,ContainerType>(ne),
        bucket_size_(ne.bucket_size_),
        window_size_(ne.window_size_)
    {
      param_ = ne.param_;

    }
    ///
    virtual ~DSignalToNoiseEstimatorWindowing() {}
    //@}


    /** @name Assignement
     */
    //@{
    ///
    inline DSignalToNoiseEstimatorWindowing& operator=(const DSignalToNoiseEstimatorWindowing& ne)
    {
      bucket_size_= ne.bucket_size_;
      window_size_=ne.window_size_;
      mz_dim_=ne.mz_dim_;
      rt_dim_=ne.rt_dim_;
      param_=ne.param_;
      first_=ne.first_;
      last_=ne.last_;
      return *this;
    }
    //@}


    /** Accessors
     */
    //@{
    /// Non-mutable access to the bucket size
    inline const int& getBucketSize() const { return bucket_size_; }
    /// Mutable access to the bucket size
    inline int& getBucketSize() { return bucket_size_; }
    /// Mutable access to bucket size
    inline void setBucketSize(const int& bucket_size) { bucket_size_ = bucket_size; }

    /// Non-mutable access to the window size
    inline const int& getWindowSize() const { return window_size_; }
    /// Mutable access to the bucket size
    inline int& getWindowSize() { return window_size_; }
    /// Mutable access to the window size
    inline void setWindowSize(const int& window_size) { window_size_ = window_size; }
    //@}


    /** @name Initialisation of the raw data intervall and estimation of noise and baseline levels
     */
    //@{
    void init(ConstIterator it_begin, ConstIterator it_end)
    {
      first_=it_begin;
      last_=it_end;
      float intervall_origin = first_->getPosition()[mz_dim_];
      float intervall_end = (last_-1)->getPosition()[mz_dim_];

      int number_of_buckets = (int)floor(fabs(intervall_origin - intervall_end) / bucket_size_) + 1;
      int buckets_per_win = (int)(window_size_ / bucket_size_);

      int length = distance(first_, last_);
      std::vector<float> Z(number_of_buckets, std::numeric_limits<float>::max());
      std::vector<float> Y(number_of_buckets, -1.*(std::numeric_limits<float>::max() - 10));
      std::vector<float> W(number_of_buckets, 0);
      std::vector<float> w(number_of_buckets, 0);

      y_base_.resize(number_of_buckets, 0);
      y_noise_.resize(number_of_buckets, 0);

      ConstIterator it_help = first_;

      for (int i = 0; i < length; i++)
      {
        int bucket = (int)floor(fabs(it_help->getPosition()[mz_dim_]-intervall_origin) / bucket_size_);
        float value = it_help->getIntensity();

        if (value > Y[bucket])
          Y[bucket] = value;

        if (value < Z[bucket])
          Z[bucket] = value;

        it_help++;
      }

      // Now iterate over all buckets and compute their W values
      for (int i = 0; i < number_of_buckets; i++)
      {
        W[i] = Y[i] - Z[i];
      }

      // Iterate again over all buckets and compute their w - values and their
      // baseline and noise contribution
      for (int i = 0; i < number_of_buckets; i++)
      {
        // starting from this bucket, sum up buckets_per_win to the left and right
        float w_value = 0;

        int start = ((i-buckets_per_win) < 0) ? 0 : (i-buckets_per_win);
        int end = ((i+buckets_per_win) >= number_of_buckets) ? number_of_buckets-1 : (i+buckets_per_win);

        for (; start < end; ++start)
        {
          w_value += (W[start] != 0) ? 1./(W[start]*W[start]) : 0.;
        }


        // now we can compute w_i
        w[i] = ( (W[i] * w_value) != 0.) ? 1. / (W[i] * W[i] * w_value) : 0;

        // and finally we can iterate over the buckets again to build y_base and y_noise
        float y_base_value  = 0.;
        float y_noise_value = 0.;

        start = ((i-buckets_per_win) < 0) ? 0 : (i-buckets_per_win);
        for (; start < end; ++start)
        {
          y_base_value += w[start] * Z[start];
          y_noise_value += w[start] * Y[start];
        }

        y_base_[i]  = y_base_value;
        y_noise_[i] = y_noise_value;
      }
    }
    //@}

    /** @name Signal To Noise Estimation
     */
    //@{
    ///
    double getSignalToNoise(ConstIterator data_point) throw (Exception::OutOfRange)
    {
      if ((data_point->getPosition()[mz_dim_] < first_->getPosition()[mz_dim_]) && (data_point->getPosition()[mz_dim_] <= (last_-1)->getPosition()[mz_dim_]))
      {
        throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      int bucket = (int)floor((data_point->getPosition()[mz_dim_]-first_->getPosition()[mz_dim_])/bucket_size_);

      // TODO: Remove this workaround.
      // if the s/n ratio for the first peak in a scan
      // is request, bucket is set to -1 which results
      // in meaningless results further below. (ost)
      if (bucket < 0) bucket = 0;

      float sn = (fabs(y_noise_[bucket] - y_base_[bucket]) > 0.0001)
                 ? (data_point->getIntensity() - y_base_[bucket]) / (y_noise_[bucket] - y_base_[bucket])
                 : 0.; // something went wrong!

      return sn;
    }
    //@}

  protected:
    /// Baseline levels for every bucket in the intervall
    std::vector<float> y_base_;
    /// Noise level for every bucket in the intervall
    std::vector<float> y_noise_;
    /// Number of data points which belong to one bucket
    int bucket_size_;
    /// Number of data points which belong to the window
    int window_size_;
  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORWINDOWING_H
