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

#ifndef OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
#define OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H


#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <cmath>
#include <map>

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise ratio of each data point in a scan by using the median (histogram based)
   
    For each datapoint in the given scan, we collect a range of data points around it (param: WindowLength).
    The noise for a datapoint is estimated to be the median of the intensities of the current window.
    If the number of elements in the current window is not sufficient (param: MinReqElementsInWindow),
    the noise level is set to a default value (param: NoiseEmptyWindow).
    The whole computation is histogram based, so the user will need to supply a number of bins (param: BinCount), which determines
    the level of error and runtime. The maximal intensity for a datapoint to be included in the histogram can be either determined 
    automatically (params: AutoMaxIntensity, AutoMode) by two different methods or can be set directly by the user (param: MaxIntensity).
    If provided, the MaxIntensity param is always favoured over AutoMaxIntensity&AutoMode. If the (estimated) MaxIntensity value is too
    low and the median is found to be in the last (&highest) bin, a warning to std:err will be given. In this case you should increase
    MaxIntensity (and optionally the BinCount).
    
    The class uses a lazy evaluation approach and will compute ALL StN values when the first request is made. 
    Changing any of the parameters will invalidate the StN values (which will invoke a recomputation on the next request).

    Accepted Ini-File Parameters:
    <table>
    <tr><th>Parameter                                               </th><th>  Description                                </th></tr>
    <tr><td> WindowLength          </td><td>  window length in Thomson                   </td></tr>
    <tr><td> BinCount              </td><td>  number of bins used                        </td></tr>
    <tr><td> MinReqElementsInWindow</td><td>  minimum number of elements required in a window   </td></tr>
    <tr><td> NoiseEmptyWindow      </td><td>  noise value used for sparse windows        </td></tr>   
      
    <tr><td> MaxIntensity          </td><td> maximal intensity used for histogram construction
                                                                   If this parameter is given, the "AutoMode/AutoMaxIntensity" is ignored 
                                                                   Be aware though that choosing an adverse "MaxIntensity" might lead to bad results
                                                                   All intensities EQUAL/ABOVE "MaxIntensity" will not be added to the histogram. If you 
                                                                   choose "MaxIntensity" too small, the noise estimate might be too small as well.
                                                                   If chosen too big, the bins become quite large (which you could counter by increasing
                                                                   "BinCount", which increases runtime). 
                                                                                                                          </td></tr>
    <tr><td colspan=2> [The following 2 parameters belong together and should both be provided if <b>MaxIntensity</b> was NOT given]</td></tr>

    <tr><td> AutoMode                 </td><td> method to use to determine "MaxIntensity": 
                                                             <pre>0 [default]  use mean + "AutoMaxIntensity" * stdev</pre>  
                                                             <pre>1            "AutoMaxIntensity"th percentile</pre>
                                                                                                                          </td></tr>
    <tr><td> AutoMaxIntensity         </td><td> parameter for "MaxIntensity" estimation; 
                                                             <pre>if "AutoMode" == 0, default is "AutoMaxIntensity"=3</pre> 
                                                             <pre>else            the default is "AutoMaxIntensity"=95 (do NOT exceed 100)</pre>
                                                                                                                          </td></tr>
    </table>  
    
        
   @note 
    Warning to *stderr* if sparse_window_percent > 20
            - percent of windows that have less than "MinReqElementsInWindow" of elements
              (noise estimates in those windows are simply a constant "NoiseEmptyWindow")
            .             
    Warning to *stderr* if histogram_oob_percent (oob=_out_of_bounds) > 1
            - percentage of median estimations that had to rely on the last(=rightmost) bin
              which gives an unreliable result
            .  
    
    
    @todo use DefaultParamHandler (Chris)
    
    @ingroup SignalToNoiseEstimators
    
  */
  
  
  template <Size D = 1 , typename PeakIterator = MSSpectrum<DRawDataPoint<1> >::const_iterator >
     class DSignalToNoiseEstimatorMedian : public DSignalToNoiseEstimator<D, PeakIterator>
  {

  public:
  
    /// method to use for estimating the maximal intensity that is used for histogram calculation
    enum IntensityThresholdCalculation { AUTOMAXBYSTDEV, AUTOMAXBYPERCENT };
    
    /// define default values
    enum DEFAULTS 
    {
       DEFAULT_WINLEN   = 200,  // window length in Thomson
       DEFAULT_BINCOUNT = 30,  // number of bins for histogramm 
       DEFAULT_MAXINTENSITY_BYSTDEV   = 3,  // "mean+DEFAULT_MAXINTENSITY_BYSTDEV*stdev" will give the treshold during preprocessing on whole dataset
       DEFAULT_MAXINTENSITY_BYPERCENT = 95, // 95th percentile value will give the treshold during preprocessing on whole dataset
       DEFAULT_MIN_REQUIRED_ELEMENTS  = 10, // number of elements required for a window to be considered NON-empty
       DEFAULT_NOISE_ON_EMTPY_WINDOW  = 2   // noise value if number of elements in a window is not sufficient 
    };
    
      
    /** @name Type definitions
     */
    //@{
    ///
    typedef typename PeakIterator::value_type PeakType;
    //@}

    using DSignalToNoiseEstimator<D, PeakIterator>::param_;
    using DSignalToNoiseEstimator<D, PeakIterator>::first_;
    using DSignalToNoiseEstimator<D, PeakIterator>::last_;
    using DSignalToNoiseEstimator<D, PeakIterator>::rt_dim_;
    using DSignalToNoiseEstimator<D, PeakIterator>::mz_dim_;

    enum DimensionID
    {
      RT = DimensionDescription < LCMS_Tag >::RT,
      MZ = DimensionDescription < LCMS_Tag >::MZ
    };

    
    /// (default) constructor
    /// accepts up to 4 parameters (for a description see above)
    /// This constructor does not accept a "max_intensity" parameter, so the algorithm will
    /// use its default heuristic to estimate a "max_intensity"
    inline DSignalToNoiseEstimatorMedian(
            double win_len               = (double) DEFAULT_WINLEN,
            int    bin_count             = DEFAULT_BINCOUNT,
            int    min_required_elements = DEFAULT_MIN_REQUIRED_ELEMENTS,
            double noise_for_empty_window= (double) DEFAULT_NOISE_ON_EMTPY_WINDOW)
       : DSignalToNoiseEstimator<D,PeakIterator>(), 
         win_len_(win_len), 
         bin_count_(bin_count),
         min_required_elements_(min_required_elements),
         noise_for_empty_window_(noise_for_empty_window)
    {
        // we were not given a max_intensity_, so we will need to calculate one later using AUTOMAXBYSTDEV
        max_intensity_ = -1;
        auto_max_intensity_ = DEFAULT_MAXINTENSITY_BYSTDEV;
        auto_mode_ = AUTOMAXBYSTDEV;
        
        updateParam();
    }
   
    /// constructor
    /// accepts 5 Parameters (for a description see above)
    /// As "max_intensity" is provided by the user, no heuristic will be used
    /// Be aware though, that a misjudged "max_intensity" (either too small or too big) 
    /// can lead to unreliable results
    inline DSignalToNoiseEstimatorMedian(
            double win_len,
            int    bin_count,
            int    min_required_elements,
            double noise_for_empty_window,
            int    max_intensity)
       : DSignalToNoiseEstimator<D,PeakIterator>(), 
         win_len_(win_len), 
         bin_count_(bin_count),
         min_required_elements_(min_required_elements),
         noise_for_empty_window_(noise_for_empty_window),
         max_intensity_((double) max_intensity)
    {
        // we were given a max_intensity_ by the user, so we can directly go about constructing the histogram later on  
        // -nevertheless we set the other options         
        auto_max_intensity_ = DEFAULT_MAXINTENSITY_BYSTDEV;
        auto_mode_ = AUTOMAXBYSTDEV;
        
        updateParam();
    }


    /// constructor
    /// accepts 6 Parameters (for a description see above)
    /// "max_intensity" will be estimated using a heuristic,
    /// which parameters you can specify using
    /// "auto_max_intensity" and "auto_mode".
    inline DSignalToNoiseEstimatorMedian(
            double win_len,
            int    bin_count,
            int    min_required_elements,
            double noise_for_empty_window,
            double auto_max_intensity,
            int    auto_mode = AUTOMAXBYSTDEV)
       : DSignalToNoiseEstimator<D,PeakIterator>(), 
         win_len_(win_len), 
         bin_count_(bin_count),
         min_required_elements_(min_required_elements),
         noise_for_empty_window_(noise_for_empty_window),
         auto_max_intensity_(auto_max_intensity), 
         auto_mode_(auto_mode)
    {
        // the user did not decide on a max_intensity_ but gave parameters on how to calculate one
        max_intensity_ = -1.0;

        updateParam();
    }
    
    /// constructor (using a param object)
    /// @see setParam(const Param& param) 
    inline DSignalToNoiseEstimatorMedian(const Param& parameters)
      : DSignalToNoiseEstimator<D,PeakIterator>(parameters)
    {
      setParam(parameters);

      is_result_valid_ = false;
    }


    /// Copy Constructor
    inline DSignalToNoiseEstimatorMedian(const DSignalToNoiseEstimatorMedian&  ne)
      : DSignalToNoiseEstimator<D,PeakIterator>(ne),
        stn_estimates_         (ne.stn_estimates_),
        win_len_               (ne.win_len_),
        bin_count_             (ne.bin_count_),
        min_required_elements_ (ne.min_required_elements_),
        noise_for_empty_window_(ne.noise_for_empty_window_),
        max_intensity_         (ne.max_intensity_),
        auto_max_intensity_    (ne.auto_max_intensity_), 
        auto_mode_             (ne.auto_mode_),
        is_result_valid_       (ne.is_result_valid_)
    {
      param_  = ne.getParam();
    }


    /// Destructor
    virtual ~DSignalToNoiseEstimatorMedian() {}


    /** @name Assignment
     * WARNING: does not assign the StNs
     */
    //@{
    ///
    inline DSignalToNoiseEstimatorMedian& operator=(const DSignalToNoiseEstimatorMedian& ne)
    {
      mz_dim_                 = ne.mz_dim_;
      rt_dim_                 = ne.rt_dim_;
      first_                  = ne.first_;
      last_                   = ne.last_;
      
      win_len_                = ne.win_len_; 
      bin_count_              = ne.bin_count_;
      min_required_elements_  = ne.min_required_elements_;
      noise_for_empty_window_ = ne.noise_for_empty_window_;
      max_intensity_          = ne.max_intensity_;
      auto_max_intensity_     = ne.auto_max_intensity_;  
      auto_mode_              = ne.auto_mode_;
      param_                  = ne.getParam();
      is_result_valid_        = ne.is_result_valid_;
      stn_estimates_          = ne.stn_estimates_;
      
      return *this;
    }
    //@}

    /** Accessors
     */
    
    /* new accessors */

    /// Non-mutable access to the window length (in Thomson)
    inline const double& getWinLen() const { return win_len_; }
    /// Mutable access to the window length (in Thomson)
    inline double& getWinLen() { return win_len_; }
    /// Mutable access to the window length (in Thomson)
    inline void setWinLen(const double& win_len) { win_len_ = win_len; updateParam();}
     
    /// Non-mutable access to the number of bins used for the histogram (the more bins, the better the approximation, but longer runtime)
    inline const int& getBinCount() const { return bin_count_; }
    /// Mutable access to the number of bins used for the histogram
    inline int& getBinCount() { return bin_count_; }
    /// Mutable access to the number of bins used for the histogram
    inline void setBinCount(const int& bin_count) { bin_count_ = bin_count; updateParam();}

    /// Non-mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
    inline const double& getMaxIntensity() const { return max_intensity_; }
    /// Mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
    inline double& getMaxIntensity() { return max_intensity_; }
    /// Mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
    inline void setMaxIntensity(const double& max_intensity) { max_intensity_ = max_intensity; updateParam();}
    
    /// Non-mutable access to the AutoMaxIntensity-Param, which holds either a percentile or a stdev-multiplier (depending on AutoMode)
    inline const double& getAutoMaxIntensity() const { return auto_max_intensity_; }
    /// Mutable access to the AutoMaxIntensity-Param, which holds either a percentile or a stdev-multiplier (depending on AutoMode)
    inline double& getAutoMaxIntensity() { return auto_max_intensity_; }
    /// Mutable access to the AutoMaxIntensity-Param, which holds either a percentile or a stdev-multiplier (depending on AutoMode)
    inline void setAutoMaxIntensity(const double& auto_max_intensity) { auto_max_intensity_ = auto_max_intensity; updateParam();}
    
    /// Non-mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
    inline const int& getAutoMode() const { return auto_mode_; }
    /// Mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
    inline int& getAutoMode() { return auto_mode_; }
    /// Mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
    inline void setAutoMode(const int& auto_mode) { auto_mode_ = auto_mode; updateParam();}
    
    /// Non-mutable access to the minimum required elements in a window, to be evaluated.
    inline const int& getMinReqElements() const { return min_required_elements_; }
    /// Mutable access to the minimum required elements in a window, to be evaluated.
    inline int& getMinReqElements() { return min_required_elements_; }
    /// Mutable access to the minimum required elements in a window, to be evaluated.
    inline void setMinReqElements(const int& min_required_elements) { min_required_elements_ = min_required_elements; updateParam();}
    
    /// Non-mutable access to the noise value that is used if a window contains not enough elements 
    inline const double& getNoiseForEmtpyWindow() const { return noise_for_empty_window_; }
    /// Mutable access to the noise value that is used if a window contains not enough elements
    inline double& getNoiseForEmtpyWindow() { return noise_for_empty_window_; }
    /// Mutable access to the noise value that is used if a window contains not enough elements
    inline void setNoiseForEmtpyWindow(const double& noise_for_empty_window) { noise_for_empty_window_ = noise_for_empty_window; updateParam();}


    /* overridden members */

    /// Mutable access to the mz dimension
    inline void setMZdim(const int& mz_dim) { DSignalToNoiseEstimator<D, PeakIterator>::setMZdim(mz_dim); is_result_valid_ = false;}
    /// Mutable access to the rt dimension
    inline void setRTdim(const int& rt_dim) { DSignalToNoiseEstimator<D, PeakIterator>::setRTdim(rt_dim); is_result_valid_ = false;}
    /// Mutable access to the first raw data point
    inline void setFirstDataPoint(const PeakIterator& first) { DSignalToNoiseEstimator<D, PeakIterator>::setFirstDataPoint(first); is_result_valid_ = false;}
    /// Mutable access to the last raw data point
    inline void setLastDataPoint(const PeakIterator& last) { DSignalToNoiseEstimator<D, PeakIterator>::setLastDataPoint(last); is_result_valid_ = false;}    
    
    /// Mutable access to the parameter object
    /// For valid parameter names see class description
    /// Left out parameters will be set to default
    /// @note a waring will be printed if unknown parameters are found
    inline void setParam(const Param& param) 
    { 
      param_ = param; 

      // print Warning message if there are unknown parameters (typos?)
      // ... define a default object
      Param default_p;
      default_p.setValue("WindowLength", win_len_);
      default_p.setValue("BinCount", (double) bin_count_);
      default_p.setValue("MinReqElementsInWindow", (double) min_required_elements_);
      default_p.setValue("NoiseEmptyWindow", noise_for_empty_window_);
      default_p.setValue("MaxIntensity", max_intensity_);
      default_p.setValue("AutoMode", (double) auto_mode_);
      default_p.setValue("AutoMaxIntensity", auto_max_intensity_);

      param_.setDefaults(default_p);

      // ... and check it against current param object:
      param_.checkDefaults("DSignalToNoiseEstimatorMedian",default_p);
			
      // set local members
      win_len_ = (double) param_.getValue("WindowLength");

      bin_count_ = (int) param_.getValue("BinCount");
      
      min_required_elements_ = (int) param_.getValue("MinReqElementsInWindow");
      
      noise_for_empty_window_ = (double) param_.getValue("NoiseEmptyWindow");
      
      max_intensity_ = (double) param_.getValue("MaxIntensity");
      
      auto_mode_ = (int) param_.getValue("AutoMode");
      
      auto_max_intensity_ = (double) param_.getValue("AutoMaxIntensity");
      
      is_result_valid_ = false;

    }

    /// Initialisation of the raw data interval and estimation of noise and baseline levels
    /// @note you can provide a "flat" 2D DPeakArray iterator. The class will evaluate the datapoints 
    /// scan by scan. A scan is found by collecting all items with the same (*iterator).getPosition()[rt_dim_]
    void init(PeakIterator it_begin, PeakIterator it_end)
    {
      first_= it_begin;
      last_= it_end;
      
      is_result_valid_ = false;
    }

    /// Return to signal/noise estimate for data point @p data_point
    /// @note the first query to this function will take longer, as
    ///        all SignalToNoise values are calculated
    /// @note you will get a warning to stderr if more than 20% of the 
    ///        noise estimates used sparse windows
    /// @note you will get a warning to stderr if more than 1% of the 
    ///        median estimates was in the rightmost histogram bin
    ///        (you should consider increasing MaxIntensity)
    virtual double getSignalToNoise(PeakIterator data_point)
    {
      double sparse_window = 0;
      double histogram_oob = 0;
      
      if (!is_result_valid_)
      {
        if (rt_dim_ == -1) {
            // D==1
            shiftWindow_(first_, last_, sparse_window, histogram_oob);
        }
        else
        {
            // D==2
            PeakIterator scan_start = first_;
            PeakIterator scan_end = first_;
            double scan_count = 0;
            double rt;
            double sparse_window_scan;
            double histogram_oob_scan;
            while (scan_end != last_)
            {
              scan_count++;
              rt = (*scan_end).getPosition()[rt_dim_];
              // find end of scan
              while (scan_end != last_ && rt == (*scan_end).getPosition()[rt_dim_])
              {
                ++scan_end;
              }
              shiftWindow_(scan_start, scan_end, sparse_window_scan, histogram_oob_scan);
              
              scan_start = scan_end;
              
              // add up percentages for sparse windows & histogram overflow in scans
              sparse_window += sparse_window_scan;
              histogram_oob += histogram_oob_scan;
            }
            // average result over scans
            sparse_window = sparse_window / scan_count;
            histogram_oob = histogram_oob / scan_count;
        }
        
        is_result_valid_ = true;

        // warn if percentage of sparse windows is above 20%
        if (sparse_window > 20) 
        {
          std::cerr << "WARNING in DSignalToNoiseEstimatorMedian: " 
                   << sparse_window 
                   << "% of all windows were sparse. You should consider decreasing WindowLength and/or MinReqElementsInWindow" 
                   << std::endl;
        }
        
        // warn if percentage of possibly wrong median estimates is above 1%
        if (histogram_oob > 1) 
        {
          std::cerr << "WARNING in DSignalToNoiseEstimatorMedian: " 
                   << histogram_oob 
                   << "% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. " 
                   << "You should consider increasing MaxIntensity (and maybe BinCount with it, to keep bin width reasonable)" 
                   << std::endl;
        }
        
      }  
      
      return stn_estimates_[*data_point];
    }

  protected:

    /// this function should be called whenever a parameter was changed to keep the param object up to date
    void updateParam()  
    {
      // update param-object
      param_.setValue("WindowLength", win_len_);
      param_.setValue("BinCount", (double) bin_count_);
      param_.setValue("MinReqElementsInWindow", (double) min_required_elements_);
      param_.setValue("NoiseEmptyWindow", noise_for_empty_window_);
      param_.setValue("MaxIntensity", max_intensity_);
      param_.setValue("AutoMode", (double) auto_mode_);
      param_.setValue("AutoMaxIntensity", auto_max_intensity_);
      
      is_result_valid_ = false;
    }
  
    /** 
      @brief protected struct to store parameters my, sigma for a gaussian distribution
      
      Accessors are : mean and variance      
    */ 
    struct GaussianEstimate 
    {
      double mean;
      double variance;
    };

    /// calculate StN values for all datapoints given, by using a sliding window approach
    /// @param scan_first_ first element in the scan
    /// @param scan_last_ last element in the scan (disregarded)
    /// @param sparse_window_percent percent of windows that have less than "min_required_elements_" of elements
    ///          (noise estimates in those windows are simply a constant "noise_for_empty_window_")           
    /// @param histogram_oob_percent (oob=_out_of_bounds) percentage of median estimations that had to rely on the last(=rightmost) bin
    ///          which gives an unreliable result
    void shiftWindow_(const PeakIterator& scan_first_, const PeakIterator& scan_last_, double &sparse_window_percent, double &histogram_oob_percent)
    {
      // reset counter for sparse windows
      sparse_window_percent = 0;
      // reset counter for histogram overflow
      histogram_oob_percent = 0;
      
      // reset the results
      stn_estimates_.clear();
      
      // local maximal intensity for the scan. We do not change "max_intensity_", because we might have
      // more scans to work on and the result of auto_mode would be different
      double max_intensity_local;
      
      // maximal range of histogram needs to be calculated first
      if (max_intensity_ == -1.0)
      {
        if (auto_mode_ == AUTOMAXBYSTDEV)
        {
          // use MEAN+auto_max_intensity_*STDEV as threshold
          GaussianEstimate gauss_global = estimate(scan_first_, scan_last_);
          max_intensity_local = gauss_global.mean + std::sqrt(gauss_global.variance)*auto_max_intensity_;
        } else
        {
          // get value at "auto_max_intensity_"th percentile
          // we use a histogram approach here as well.
          if (auto_max_intensity_ > 100) {auto_max_intensity_ = DEFAULT_MAXINTENSITY_BYPERCENT;}
          
          std::vector <int> histogram_auto(100, 0);
          
          // find maximum of current scan
          int size = 0;
          double maxInt = 0;
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
          int elements_below_percentile = (int) (auto_max_intensity_ * size / 100);
          int elements_seen = 0;
          int i = -1;
          run = scan_first_;
          
          while (run != scan_last_ && elements_seen < elements_below_percentile) {
            ++i;
            elements_seen += histogram_auto[i];
            ++run;
          }
          
          max_intensity_local = (((double)i) + 0.5) * bin_size;
        }
      } else 
      {
        max_intensity_local = max_intensity_;
      }
      
      PeakIterator window_pos_center  = scan_first_;
      PeakIterator window_pos_borderleft = scan_first_;
      PeakIterator window_pos_borderright = scan_first_;
      
      double window_half_size = win_len_ / 2;
      double bin_size = max_intensity_local / bin_count_;
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
      int to_bin;

      // index of bin where the median is located
      int median_bin;
      // additive number of elements from left to x in histogram
      int element_inc_count;
      
      // tracks elements in current window, which may vary because of uneven spaced data
      int elements_in_window = 0;
      // number of windows
      int window_count = 0;
      
      // number of elements where we find the median
      int element_in_window_half;
      
      double noise;    // noise value of a datapoint      

      while (window_pos_center != scan_last_)
      {
        
        // erase all elements from histogram that will leave the window on the LEFT side
        while ( (*window_pos_borderleft).getPosition()[mz_dim_] <  (*window_pos_center).getPosition()[mz_dim_] - window_half_size )
        {
          to_bin = std::min((int) (((*window_pos_borderleft).getIntensity()) / bin_size), bin_count_minus_1);
          --histogram[to_bin];
          --elements_in_window;
          ++window_pos_borderleft;
        }
        
        // add all elements to histogram that will enter the window on the RIGHT side
        while (    (window_pos_borderright != scan_last_)
                &&((*window_pos_borderright).getPosition()[mz_dim_] <= (*window_pos_center).getPosition()[mz_dim_] + window_half_size ) )
        {
          to_bin = std::min((int) (((*window_pos_borderright).getIntensity()) / bin_size), bin_count_minus_1);
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
          
          noise = std::max(noise_for_empty_window_, bin_value[median_bin]);
        }
        
        // store result
        stn_estimates_[*window_pos_center] = (*window_pos_center).getIntensity() / noise;
        
        
        // advance the window center by one datapoint
        ++window_pos_center;
        ++window_count;  
      } // end while

      sparse_window_percent = sparse_window_percent *100 / window_count;
      histogram_oob_percent = histogram_oob_percent *100 / window_count;
      
    } // end of shiftWindow_

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

    /// stores the noise estimate for each peak
    ///TODO change to std::hash_map<const char*, int, hash<const char*>, eqstr>  == the current one is SORTED (WHY??)!
    std::map< PeakType, double, typename PeakType::PositionLess > stn_estimates_;
    //HashMap <PeakType, double> stn_estimates_; // --> error cast to int

    /// range of data points which belong to a window in Thomson
    double win_len_;    
    /// number of bins in the histogram
    int    bin_count_;
    /// minimal number of elements a window need to cover to be used
    int min_required_elements_;
    /// used as noise value for windows which cover less than "min_required_elements_" AND for windows where the mean value got iterated down to 0
    double noise_for_empty_window_;
    /// maximal intensity considered during binning (values above get discarded)
    double max_intensity_;
    /// parameter for initial automatic estimation of "max_intensity_": either a percentile or a stdev, interpretation depending on "auto_mode_" 
    double auto_max_intensity_;
    /// determines which method shall be used for estimating "max_intensity_". valid are AUTOMAXBYSTDEV or AUTOMAXBYPERCENT
    int    auto_mode_;

    /// flag: set to true if SignalToNoise estimates are calculated and none of the params were changed. otherwise false.
    mutable bool is_result_valid_;


  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
