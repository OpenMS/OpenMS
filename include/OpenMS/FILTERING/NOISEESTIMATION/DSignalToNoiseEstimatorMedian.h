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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
#define OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H

#include <OpenMS/KERNEL/DPeakArray.h>
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
   template <Size D = 1 , typename PeakIterator = MSSpectrum<DRawDataPoint<1> >::const_iterator >
  class DSignalToNoiseEstimatorMedian : public DSignalToNoiseEstimator<D, PeakIterator>
  {

  public:
    /** @name Type definitions
     */
    //@{
    ///
    typedef typename PeakIterator::value_type PeakType;
    ///
    typedef typename PeakType::CoordinateType CoordinateType;
    ///
    typedef typename PeakType::IntensityType IntensityType;
    //@}

    using DSignalToNoiseEstimator<D, PeakIterator>::param_;
	using DSignalToNoiseEstimator<D, PeakIterator>::first_;
    using DSignalToNoiseEstimator<D, PeakIterator>::last_;

    enum DimensionID
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    /** @name Constructors and Destructor
     */
    //@{
    inline DSignalToNoiseEstimatorMedian()
      : DSignalToNoiseEstimator<D,PeakIterator>(),
        window_size_(100), median_perc_(1) {}
    ///
    inline DSignalToNoiseEstimatorMedian(const Param& parameters)
      : DSignalToNoiseEstimator<D,PeakIterator>(parameters)
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
      : DSignalToNoiseEstimator<D,PeakIterator>(ne),
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
      param_           = ne.param_;
      return *this;
    }
    //@}


    /** Accessors
     */
    //@{
    /// Non-mutable access to the window size
    inline const unsigned int& getWindowSize() const { return window_size_; }
    /// Mutable access to the window size
    inline unsigned int& getWindowSize() { return window_size_; }
    /// Mutable access to the window size
    inline void setWindowSize(const int wsize) { window_size_ = wsize; }
    /// Non-mutable access to the factor
    inline const float& getFactor() const { return median_perc_; }
    /// Mutable access to the factor
    inline float& getFactor() { return median_perc_; }
    /// Mutable access to the factor
    inline void setFactor(const float factor) { median_perc_ = factor; }
    //@}


    /**
    	@brief Initialisation of the raw data interval and estimation of noise and baseline levels
			
			@note Works only scanwise i.e. you have to pass your data to this function scan by scan.
    */
    void init(PeakIterator it_begin, PeakIterator it_end)
    {
	  first_= it_begin;
      last_= it_end;
	
      std::vector<IntensityType> intensities(window_size_);
      DPeakArray<D,PeakType> scan;
      // collect all intensities and compute the background noise as
	  // the median of the intensities in a small window.
	  while (it_begin != it_end)
      {
        scan.push_back(*it_begin);
        it_begin++;
      }
	  
	  shiftWindow_(scan);
      scan.clear();
    }

    /// Return to signal/noise estimate for data point @p data_point
    double getSignalToNoise(PeakIterator data_point)
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

    void shiftWindow_(DPeakArray<D,PeakType> current_scan)
    {
      unsigned int left  = (int) floor(window_size_ / 2);

      for (unsigned int i=0; i<current_scan.size(); i++)
      {
        DPeakArray<D,PeakType> window;
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
