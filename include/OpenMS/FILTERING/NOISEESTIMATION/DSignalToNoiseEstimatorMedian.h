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
// $Id: DSignalToNoiseEstimatorMedian.h,v 1.7 2006/04/21 14:54:31 elange Exp $
// $Author: elange $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
#define OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Types.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

namespace OpenMS
{
  /**
  	@brief Simple noise estimator, estimates the signal/noise ratio of each data point in a scan
           based on the median over a small m/z range.
   
    For each datapoint in the map, we collect a range of data points around it (in the same scan).
    We estimate the s/n ratio as the median of the intensities of the points in this range.
    The width of this range is give by window_size_ .
  */
  template <Size D, typename  ContainerType = DPeakArrayNonPolymorphic<D,DRawDataPoint<D> > >
  class DSignalToNoiseEstimatorMedian : public DSignalToNoiseEstimator<D, ContainerType>
  {

  public:
    /** @name Type definitions
     */
    //@{
    ///
    typedef typename ContainerType::const_iterator ConstIterator;
    ///
    typedef typename ContainerType::PeakType PeakType;
    ///
    typedef typename ContainerType::PeakType::CoordinateType CoordinateType;
    ///
    typedef typename ContainerType::PeakType::IntensityType IntensityType;
    //@}

    using DSignalToNoiseEstimator<D, ContainerType>::param_;

    enum DimensionID
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    /** @name Constructors and Destructor
     */
    //@{
    inline DSignalToNoiseEstimatorMedian()
      : DSignalToNoiseEstimator<D,ContainerType>(),
        window_size_(100) {}
    ///
    inline DSignalToNoiseEstimatorMedian(const Param& parameters)
      : DSignalToNoiseEstimator<D,ContainerType>(parameters)
    {
      // set params
      DataValue dv = param_.getValue("SignalToNoiseEstimationParameter:Window");
      if (dv.isEmpty() || dv.toString() == "") window_size_ = 100;
      else window_size_ = (int)dv;

      dv = param_.getValue("SignalToNoiseEstimationParameter:Median_perc");
      if (dv.isEmpty() || dv.toString() == "") median_perc_ = 1.0;
      else median_perc_ = (float)dv;
    }
    ///
    inline DSignalToNoiseEstimatorMedian(const DSignalToNoiseEstimatorMedian&  ne)
      : DSignalToNoiseEstimator<D,ContainerType>(ne),
        window_size_(ne.window_size_),
        median_perc_(ne.median_perc_) {}
    ///
    virtual ~DSignalToNoiseEstimatorMedian() {}
    //@}


    /** @name Assignement
     */
    //@{
    ///
    inline DSignalToNoiseEstimatorMedian& operator=(const DSignalToNoiseEstimatorMedian& ne)
    {
      window_size_ = ne.window_size_;
      median_perc_ = ne.median_perc_;
      param_       = ne.param_;
      return *this;
    }
    //@}


    /** Accessors
     */
    //@{
    /// Non-mutable access to the window size
    inline const int getWindowSize() const { return window_size_; }
    /// Mutable access to the window size
    inline int getWindowSize() { return window_size_; }
    /// Mutable access to the window size
    inline void setWindowSize(const int wsize) { window_size_ = wsize; }
    /// Non-mutable access to the factor
    inline const float getFactor() const { return median_perc_; }
    /// Mutable access to the factor
    inline float getFactor() { return median_perc_; }
    /// Mutable access to the factor
    inline void setFactor(const float factor) { median_perc_ = factor; }
    //@}


    /// Initialisation of the raw data interval and estimation of noise and baseline levels

     **/

    void init(ConstIterator it_begin, ConstIterator it_end)
    {
      CoordinateType current_rt = it_begin->getPosition()[RT];
      std::vector<IntensityType> intensities(window_size_);
      ContainerType scan;
      while (it_begin != it_end)
      {
        CoordinateType next_rt = it_begin->getPosition()[RT];
        if (next_rt != current_rt)
        {
          shiftWindow_(scan);
          scan.clear();
        }
        scan.push_back(*it_begin);
        it_begin++;
      }
    }

    /// Return to signal/noise estimate for data point @p data_point
    double getSignalToNoise(ConstIterator data_point)
    {
      double noise = noise_estimates_[*data_point];

      // if the current noise estimate is zero, we
      // we set the background noise to 2.
      if (noise == 0)
        return  (data_point->getIntensity() / 2);
      else
        return data_point->getIntensity() / noise;
    }

  protected:

    void shiftWindow_(ContainerType current_scan)
    {
      unsigned int left  = (int) floor(window_size_ / 2);

      for (unsigned int i=0; i<current_scan.size(); i++)
      {
        ContainerType window(window_size_);
        // walk to the left and collect
        // and most (window_size_ / 2) peaks
        for (int j=i; j >= 0 && window.size() <= left; j--)
        {
          window.push_back(current_scan[j]);
        }

        // walk to the right and collect
        // at most window_size_ peaks
        for (unsigned int k=(i+1); k < current_scan.size() && window.size() <= window_size_; k++)
        {
          window.push_back(current_scan[k]);
        }

        // compute median of the intensities
        int middle = (int) ceil(window.size() / 2);
        window.sortByIntensity();
        noise_estimates_[current_scan[i]] = window[middle].getIntensity();

      } // end for (i in scan)

    } // end of shiftWindow_

    /// Number of data points which belong to the window
    unsigned int window_size_;

    /// Percentage of the median that we use to set the s/n ratio (default is 1)
    float median_perc_;

    /// stores the noise estimate for each peak
    std::map<PeakType,double, typename PeakType::PositionLess > noise_estimates_;

  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
