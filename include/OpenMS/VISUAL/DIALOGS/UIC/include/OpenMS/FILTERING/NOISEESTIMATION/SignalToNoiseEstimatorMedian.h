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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEDIAN_H
#define OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEDIAN_H


#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan by using the median (histogram based)
   
    For each datapoint in the given scan, we collect a range of data points around it (param: <i>win_len</i>).
    The noise for a datapoint is estimated to be the median of the intensities of the current window.
    If the number of elements in the current window is not sufficient (param: <i>MinReqElements</i>),
    the noise level is set to a default value (param: <i>noise_for_empty_window</i>).
    The whole computation is histogram based, so the user will need to supply a number of bins (param: <i>bin_count</i>), which determines
    the level of error and runtime. The maximal intensity for a datapoint to be included in the histogram can be either determined 
    automatically (params: <i>AutoMaxIntensity</i>, <i>auto_mode</i>) by two different methods or can be set directly by the user (param: <i>max_intensity</i>).
    If the (estimated) <i>max_intensity</i> value is too low and the median is found to be in the last (&highest) bin, a warning to std:err will be given. In this case you should increase
    <i>max_intensity</i> (and optionally the <i>bin_count</i>).
    
    Changing any of the parameters will invalidate the S/N values (which will invoke a recomputation on the next request).

    @note If more than 20 percent of windows have less than <i>min_required_elements</i> of elements, a warning is issued to <i>stderr</i> and noise estimates in those windows are set to the constant <i>noise_for_empty_window</i>.     
    @note If more than 1 percent of median estimations had to rely on the last(=rightmost) bin (which gives an unreliable result), a warning is issued to <i>stderr</i>.
    
		@htmlinclude OpenMS_SignalToNoiseEstimatorMedian.parameters
    
    @ingroup SignalProcessing
  */
  
  template < typename Container = MSSpectrum< > >
  class SignalToNoiseEstimatorMedian : public SignalToNoiseEstimator< Container >
  {

  public:

    /// method to use for estimating the maximal intensity that is used for histogram calculation
    enum IntensityThresholdCalculation { MANUAL=-1, AUTOMAXBYSTDEV=0, AUTOMAXBYPERCENT=1 };

    using SignalToNoiseEstimator< Container >::stn_estimates_;
    using SignalToNoiseEstimator< Container >::first_;
    using SignalToNoiseEstimator< Container >::last_;
    using SignalToNoiseEstimator< Container >::is_result_valid_;
    using SignalToNoiseEstimator< Container >::defaults_;
    using SignalToNoiseEstimator< Container >::param_;
    
    typedef typename SignalToNoiseEstimator< Container >::PeakIterator PeakIterator;
    typedef typename SignalToNoiseEstimator< Container >::PeakType PeakType;
    
    typedef typename SignalToNoiseEstimator< Container >::GaussianEstimate GaussianEstimate;

    /// default constructor
    inline SignalToNoiseEstimatorMedian()
    {
    	//set the name for DefaultParamHandler error messages
    	this->setName("SignalToNoiseEstimatorMedian");	

      defaults_.setValue("max_intensity", -1, "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode)."\
" Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')!"\
" All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin."\
" If you choose 'max_intensity' too small, the noise estimate might be too small as well. "\
" If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime)."\
" In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.", StringList::create("advanced")); 
			defaults_.setMinInt ("max_intensity", -1);

      defaults_.setValue("auto_max_stdev_factor", 3.0, "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev", StringList::create("advanced"));
			defaults_.setMinFloat ("auto_max_stdev_factor", 0.0);
			defaults_.setMaxFloat ("auto_max_stdev_factor", 999.0);
							
      defaults_.setValue("auto_max_percentile", 95, "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile", StringList::create("advanced"));
			defaults_.setMinInt ("auto_max_percentile", 0);
			defaults_.setMaxInt ("auto_max_percentile", 100);
							
      defaults_.setValue("auto_mode", 0, "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method", StringList::create("advanced"));
			defaults_.setMinInt ("auto_mode", -1);
			defaults_.setMaxInt ("auto_mode", 1);  
							
      defaults_.setValue("win_len", 200.0, "window length in Thomson");
			defaults_.setMinFloat ("win_len", 1.0);
							
      defaults_.setValue("bin_count", 30, "number of bins for intensity values");
			defaults_.setMinInt ("bin_count", 3);
							
      defaults_.setValue("min_required_elements", 10, "minimum number of elements required in a window (otherwise it is considered sparse)");
			defaults_.setMinInt ("min_required_elements", 1);
			
      defaults_.setValue("noise_for_empty_window", std::pow(10.0,20), "noise value used for sparse windows", StringList::create("advanced"));
			

      SignalToNoiseEstimator< Container >::defaultsToParam_();
    }


    /// Copy Constructor
    inline SignalToNoiseEstimatorMedian(const SignalToNoiseEstimatorMedian&  source)
        : SignalToNoiseEstimator< Container >(source)
    {
      updateMembers_();
    }


    /** @name Assignment
     */
    //@{
    ///
    inline SignalToNoiseEstimatorMedian& operator=(const SignalToNoiseEstimatorMedian& source)
    {
      if(&source == this) return *this; 
      SignalToNoiseEstimator< Container >::operator=(source);
      updateMembers_();
      return *this;
    }
    //@}


    /// Destructor
    virtual ~SignalToNoiseEstimatorMedian()
    {}


  protected:


		/** calculate StN values for all datapoints given, by using a sliding window approach
				@param scan_first_ first element in the scan
				@param scan_last_ last element in the scan (disregarded)
				@exception Throws Exception::InvalidValue
	  */
    void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_)
    {
      // reset counter for sparse windows
      double sparse_window_percent = 0;
      // reset counter for histogram overflow
      double histogram_oob_percent = 0;
      
      // reset the results
      stn_estimates_.clear();
      
      // maximal range of histogram needs to be calculated first
      if (auto_mode_ == AUTOMAXBYSTDEV)
      {
        // use MEAN+auto_max_intensity_*STDEV as threshold
        GaussianEstimate gauss_global = SignalToNoiseEstimator< Container >::estimate_(scan_first_, scan_last_);
        max_intensity_ = gauss_global.mean + std::sqrt(gauss_global.variance)*auto_max_stdev_Factor_;
      }
      else if (auto_mode_ == AUTOMAXBYPERCENT)
      {
        // get value at "auto_max_percentile_"th percentile
        // we use a histogram approach here as well.
        if ((auto_max_percentile_ < 0) || (auto_max_percentile_ > 100))
        {
          String s = auto_max_percentile_;
          throw Exception::InvalidValue(__FILE__, 
                                         __LINE__, 
                                         __PRETTY_FUNCTION__, 
                                         "auto_mode is on AUTOMAXBYPERCENT! auto_max_percentile is not in [0,100]. Use setAutoMaxPercentile(<value>) to change it!", 
                                         s);
        }

        std::vector <int> histogram_auto(100, 0);

        // find maximum of current scan
        int size = 0;
        typename PeakType::IntensityType maxInt = 0;
        PeakIterator run = scan_first_;
        while (run != scan_last_)
        {
          maxInt = std::max(maxInt, (*run).getIntensity());
          ++size;
          ++run;
        }

        double bin_size = maxInt / 100;

        // fill histogram
        run = scan_first_;
        while (run != scan_last_)
        {
          ++histogram_auto[(int) (((*run).getIntensity()-1) / bin_size)];
          ++run;
        }

        // add up element counts in histogram until ?th percentile is reached
        int elements_below_percentile = (int) (auto_max_percentile_ * size / 100);
        int elements_seen = 0;
        int i = -1;
        run = scan_first_;

        while (run != scan_last_ && elements_seen < elements_below_percentile)
        {
          ++i;
          elements_seen += histogram_auto[i];
          ++run;
        }

        max_intensity_ = (((double)i) + 0.5) * bin_size;
      }
      else //if (auto_mode_ == MANUAL)
      {
        if (max_intensity_<=0) 
        {
          String s = max_intensity_;
          throw Exception::InvalidValue(__FILE__, 
                                         __LINE__, 
                                         __PRETTY_FUNCTION__, 
                                         "auto_mode is on MANUAL! max_intensity is <=0. Needs to be positive! Use setMaxIntensity(<value>) or enable auto_mode!", 
                                         s);
        }
      }
      
			if (max_intensity_ < 0)
			{
				std::cerr << "TODO SignalToNoiseEstimatorMedian: the max_intensity_ value should be positive! " << max_intensity_ << std::endl;
				return;
			}

      PeakIterator window_pos_center  = scan_first_;
      PeakIterator window_pos_borderleft = scan_first_;
      PeakIterator window_pos_borderright = scan_first_;
      
      double window_half_size = win_len_ / 2;
      double bin_size = std::max(1.0, max_intensity_ / bin_count_); // at least size of 1 for intensity bins
      int bin_count_minus_1 = bin_count_ - 1;
      
      std::vector <int> histogram(bin_count_, 0);
      std::vector <double> bin_value(bin_count_, 0);
      // calculate average intensity that is represented by a bin
      for (int bin=0; bin<bin_count_; bin++)
      {
         histogram[bin] = 0;
         bin_value[bin] = (bin + 0.5) * bin_size;           
      }
      // bin in which a datapoint would fall
      int to_bin = 0;

      // index of bin where the median is located
      int median_bin = 0;
      // additive number of elements from left to x in histogram
      int element_inc_count = 0;
      
      // tracks elements in current window, which may vary because of unevenly spaced data
      int elements_in_window = 0;
      // number of windows
      int window_count = 0;
      
      // number of elements where we find the median
      int element_in_window_half = 0;
      
      double noise;    // noise value of a datapoint      

      // determine how many elements we need to estimate (for progress estimation)
      int windows_overall = 0;
      PeakIterator run = scan_first_;
      while (run != scan_last_)
      {
        ++windows_overall;
        ++run;
      }
      SignalToNoiseEstimator< Container >::startProgress(0,windows_overall,"noise estimation of data");

      // MAIN LOOP
      while (window_pos_center != scan_last_)
      {
        
        // erase all elements from histogram that will leave the window on the LEFT side
        while ( (*window_pos_borderleft).getMZ() <  (*window_pos_center).getMZ() - window_half_size )
        {
					to_bin = std::max( std::min<int>( (int)((*window_pos_borderleft).getIntensity() / bin_size), bin_count_minus_1), 0);
          --histogram[to_bin];
          --elements_in_window;
          ++window_pos_borderleft;
        }
        
        // add all elements to histogram that will enter the window on the RIGHT side
        while (    (window_pos_borderright != scan_last_)
                &&((*window_pos_borderright).getMZ() <= (*window_pos_center).getMZ() + window_half_size ) )
        {
					//std::cerr << (*window_pos_borderright).getIntensity() << " " << bin_size << " " << bin_count_minus_1 << std::endl;
          to_bin = std::max( std::min<int>( (int)((*window_pos_borderright).getIntensity() / bin_size), bin_count_minus_1), 0);
          ++histogram[to_bin];
          ++elements_in_window;
          ++window_pos_borderright;
        }

        if (elements_in_window < min_required_elements_)
        {
          noise = noise_for_empty_window_;
          ++sparse_window_percent;
        }
        else
        {
          // find bin i where ceil[elements_in_window/2] <= sum_c(0..i){ histogram[c] }
          median_bin = -1;
          element_inc_count = 0;
          element_in_window_half = (elements_in_window+1) / 2;
          while (median_bin < bin_count_minus_1 && element_inc_count < element_in_window_half) {
            ++median_bin;
            element_inc_count += histogram[median_bin];
          }

          // increase the error count
          if (median_bin == bin_count_minus_1) {++histogram_oob_percent;}
          
          // just avoid division by 0
          noise = std::max(1.0, bin_value[median_bin]);
        }
        
        // store result
        stn_estimates_[*window_pos_center] = (*window_pos_center).getIntensity() / noise;
        
        
        // advance the window center by one datapoint
        ++window_pos_center;
        ++window_count;  
        // update progress 
        SignalToNoiseEstimator< Container >::setProgress(window_count);
                  
      } // end while

      SignalToNoiseEstimator< Container >::endProgress();
        
      sparse_window_percent = sparse_window_percent *100 / window_count;
      histogram_oob_percent = histogram_oob_percent *100 / window_count;
      
      // warn if percentage of sparse windows is above 20%
      if (sparse_window_percent > 20) 
      {
        std::cerr << "WARNING in SignalToNoiseEstimatorMedian: " 
                 << sparse_window_percent 
                 << "% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'" 
                 << std::endl;
      }
      
      // warn if percentage of possibly wrong median estimates is above 1%
      if (histogram_oob_percent > 1) 
      {
        std::cerr << "WARNING in SignalToNoiseEstimatorMedian: " 
                 << histogram_oob_percent 
                 << "% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. " 
                 << "You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)" 
                 << std::endl;
      }      
      
    } // end of shiftWindow_

    /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_()
    {
      max_intensity_         = (double)param_.getValue("max_intensity"); 
      auto_max_stdev_Factor_ = (double)param_.getValue("auto_max_stdev_factor"); 
      auto_max_percentile_   = param_.getValue("auto_max_percentile"); 
      auto_mode_             = param_.getValue("auto_mode"); 
      win_len_               = (double)param_.getValue("win_len"); 
      bin_count_             = param_.getValue("bin_count"); 
      min_required_elements_ = param_.getValue("min_required_elements"); 
      noise_for_empty_window_= (double)param_.getValue("noise_for_empty_window"); 
      is_result_valid_ = false;
    }

    /// maximal intensity considered during binning (values above get discarded)
    double max_intensity_;
    /// parameter for initial automatic estimation of "max_intensity_": a stdev multiplier
    double auto_max_stdev_Factor_;
    /// parameter for initial automatic estimation of "max_intensity_" percentile or a stdev
    double auto_max_percentile_;
    /// determines which method shall be used for estimating "max_intensity_". valid are MANUAL=-1, AUTOMAXBYSTDEV=0 or AUTOMAXBYPERCENT=1
    int    auto_mode_;
    /// range of data points which belong to a window in Thomson
    double win_len_;
    /// number of bins in the histogram
    int    bin_count_;
    /// minimal number of elements a window needs to cover to be used
    int min_required_elements_;
    /// used as noise value for windows which cover less than "min_required_elements_" 
    /// use a very high value if you want to get a low S/N result
    double noise_for_empty_window_;



  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
