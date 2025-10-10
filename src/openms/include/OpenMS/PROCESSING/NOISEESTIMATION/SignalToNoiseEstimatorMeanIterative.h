// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <vector>
#include <algorithm> //for std::max_element

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan
           based on an iterative scheme which discards high intensities

    For each datapoint in the given scan, we collect a range of data points around it (param: <i>win_len</i>).
    The noise for a datapoint is estimated iteratively by discarding peaks which are more than
    (<i>stdev_mp</i> * StDev) above the mean value. After three iterations, the mean value is
    considered to be the noise level. If the number of elements in the current window is not sufficient (param: <i>min_required_elements</i>),
    the noise level is set to a default value (param: <i>noise_for_empty_window</i>).

    The whole computation is histogram based, so the user will need to supply a number of bins (param: <i>bin_count</i>), which determines
    the level of error and runtime. The maximal intensity for a datapoint to be included in the histogram can be either determined
    automatically (param: <i>auto_mode</i>) by two different methods or can be set directly by the user (param: <i>max_intensity</i>).

    Changing any of the parameters will invalidate the S/N values (which will invoke a recomputation on the next request).


    @note If more than 20 percent of windows have less than <i>min_required_elements</i> of elements, a warning is issued to <i>stderr</i> and noise estimates in those windows are set to the constant <i>noise_for_empty_window</i>.

    @htmlinclude OpenMS_SignalToNoiseEstimatorMeanIterative.parameters

    @ingroup SignalProcessing
  */
  template <typename Container = MSSpectrum>
  class SignalToNoiseEstimatorMeanIterative :
    public SignalToNoiseEstimator<Container>
  {

public:

    /// method to use for estimating the maximal intensity that is used for histogram calculation
    enum IntensityThresholdCalculation {MANUAL = -1, AUTOMAXBYSTDEV = 0, AUTOMAXBYPERCENT = 1};

    using SignalToNoiseEstimator<Container>::stn_estimates_;
    using SignalToNoiseEstimator<Container>::defaults_;
    using SignalToNoiseEstimator<Container>::param_;

    typedef typename SignalToNoiseEstimator<Container>::PeakIterator PeakIterator;
    typedef typename SignalToNoiseEstimator<Container>::PeakType PeakType;

    typedef typename SignalToNoiseEstimator<Container>::GaussianEstimate GaussianEstimate;


    /// default constructor
    inline SignalToNoiseEstimatorMeanIterative()
    {
      //set the name for DefaultParamHandler error messages
      this->setName("SignalToNoiseEstimatorMeanIterative");

      defaults_.setValue("max_intensity", -1, "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode)." \
                                              " Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')!" \
                                              " All intensities EQUAL/ABOVE 'max_intensity' will not be added to the histogram." \
                                              " If you choose 'max_intensity' too small, the noise estimate might be too small as well." \
                                              " If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime).", {"advanced"});
      defaults_.setMinInt("max_intensity", -1);

      defaults_.setValue("auto_max_stdev_factor", 3.0, "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev", {"advanced"});
      defaults_.setMinFloat("auto_max_stdev_factor", 0.0);
      defaults_.setMaxFloat("auto_max_stdev_factor", 999.0);


      defaults_.setValue("auto_max_percentile", 95, "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile", {"advanced"});
      defaults_.setMinInt("auto_max_percentile", 0);
      defaults_.setMaxInt("auto_max_percentile", 100);

      defaults_.setValue("auto_mode", 0, "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method", {"advanced"});
      defaults_.setMinInt("auto_mode", -1);
      defaults_.setMaxInt("auto_mode", 1);

      defaults_.setValue("win_len", 200.0, "window length in Thomson");
      defaults_.setMinFloat("win_len", 1.0);

      defaults_.setValue("bin_count", 30, "number of bins for intensity values");
      defaults_.setMinInt("bin_count", 3);

      defaults_.setValue("stdev_mp", 3.0, "multiplier for stdev", {"advanced"});
      defaults_.setMinFloat("stdev_mp", 0.01);
      defaults_.setMaxFloat("stdev_mp", 999.0);

      defaults_.setValue("min_required_elements", 10, "minimum number of elements required in a window (otherwise it is considered sparse)");
      defaults_.setMinInt("min_required_elements", 1);

      defaults_.setValue("noise_for_empty_window", std::pow(10.0, 20), "noise value used for sparse windows", {"advanced"});

      SignalToNoiseEstimator<Container>::defaultsToParam_();
    }

    /// Copy Constructor
    inline SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative & source) :
      SignalToNoiseEstimator<Container>(source)
    {
      updateMembers_();
    }

    /** @name Assignment
     */
    //@{
    ///
    inline SignalToNoiseEstimatorMeanIterative & operator=(const SignalToNoiseEstimatorMeanIterative & source)
    {
      if (&source == this) return *this;

      SignalToNoiseEstimator<Container>::operator=(source);
      updateMembers_();
      return *this;
    }

    //@}


    /// Destructor
    ~SignalToNoiseEstimatorMeanIterative() override
    {}


protected:


    /** calculate StN values for all datapoints given, by using a sliding window approach
                  @param c raw data
                  @exception Throws Exception::InvalidValue
           */
    void computeSTN_(const Container& c) override
    {
      //first element in the scan
      PeakIterator scan_first_ = c.begin();
      //last element in the scan
      PeakIterator scan_last_ = c.end();

      // reset counter for sparse windows
      double sparse_window_percent = 0;

      // reset the results
      stn_estimates_.clear();
      stn_estimates_.resize(c.size());

      // maximal range of histogram needs to be calculated first
      if (auto_mode_ == AUTOMAXBYSTDEV)
      {
        // use MEAN+auto_max_intensity_*STDEV as threshold
        GaussianEstimate gauss_global = SignalToNoiseEstimator<Container>::estimate_(scan_first_, scan_last_);
        max_intensity_ = gauss_global.mean + std::sqrt(gauss_global.variance) * auto_max_stdev_Factor_;
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
                                        OPENMS_PRETTY_FUNCTION,
                                        "auto_mode is on AUTOMAXBYPERCENT! auto_max_percentile is not in [0,100]. Use setAutoMaxPercentile(<value>) to change it!",
                                        s);
        }

        std::vector<int> histogram_auto(100, 0);

        // find maximum of current scan
        auto maxIt = std::max_element(c.begin(), c.end() ,[](const PeakType& a, const PeakType& b){ return a.getIntensity() > b.getIntensity();});
        typename PeakType::IntensityType maxInt = maxIt->getIntensity();

        double bin_size = maxInt / 100;

        // fill histogram
        for(auto& run : c)
        {
          ++histogram_auto[(int) (((run).getIntensity() - 1) / bin_size)];
        }

        // add up element counts in histogram until ?th percentile is reached
        int elements_below_percentile = (int) (auto_max_percentile_ * c.size() / 100);
        int elements_seen = 0;
        int i = -1;
        PeakIterator run = scan_first_;

        while (run != scan_last_ && elements_seen < elements_below_percentile)
        {
          ++i;
          elements_seen += histogram_auto[i];
          ++run;
        }

        max_intensity_ = (((double)i) + 0.5) * bin_size;
      }
      else   //if (auto_mode_ == MANUAL)
      {
        if (max_intensity_ <= 0)
        {
          String s = max_intensity_;
          throw Exception::InvalidValue(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
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
      double bin_size = std::max(1.0, max_intensity_ / bin_count_);   // at least size of 1 for intensity bins

      std::vector<int> histogram(bin_count_, 0);
      std::vector<double> bin_value(bin_count_, 0);
      // calculate average intensity that is represented by a bin
      for (int bin = 0; bin < bin_count_; bin++)
      {
        histogram[bin] = 0;
        bin_value[bin] = (bin + 0.5) * bin_size;
      }
      // index of last valid bin during iteration
      int hist_rightmost_bin;
      // bin in which a datapoint would fall
      int to_bin;
      // mean & stdev of the histogram
      double hist_mean;
      double hist_stdev;

      // tracks elements in current window, which may vary because of unevenly spaced data
      int elements_in_window = 0;
      int window_count = 0;

      double noise;      // noise value of a datapoint

      ///start progress estimation
      SignalToNoiseEstimator<Container>::startProgress(0, c.size(), "noise estimation of data");

      // MAIN LOOP
      while (window_pos_center != scan_last_)
      {
        // erase all elements from histogram that will leave the window on the LEFT side
        while ((*window_pos_borderleft).getMZ() <  (*window_pos_center).getMZ() - window_half_size)
        {
          //std::cout << "S: " << (*window_pos_borderleft).getMZ()  <<  " " << ( (*window_pos_center).getMZ() - window_half_size ) << "\n";
          to_bin = (int) ((std::max((*window_pos_borderleft).getIntensity(), 0.0f)) / bin_size);
          if (to_bin < bin_count_)
          {
            --histogram[to_bin];
            --elements_in_window;
          }
          ++window_pos_borderleft;
        }

        //std::printf("S1: %E %E\n", (*window_pos_borderright).getMZ(), (*window_pos_center).getMZ() + window_half_size);


        // add all elements to histogram that will enter the window on the RIGHT side
        while ((window_pos_borderright != scan_last_)
              && ((*window_pos_borderright).getMZ() < (*window_pos_center).getMZ() + window_half_size))
        {
          //std::printf("Sb: %E %E %E\n", (*window_pos_borderright).getMZ(), (*window_pos_center).getMZ() + window_half_size, (*window_pos_borderright).getMZ() - ((*window_pos_center).getMZ() + window_half_size));

          to_bin = (int) ((std::max((*window_pos_borderright).getIntensity(), 0.0f)) / bin_size);
          if (to_bin < bin_count_)
          {
            ++histogram[to_bin];
            ++elements_in_window;
          }
          ++window_pos_borderright;
        }

        if (elements_in_window < min_required_elements_)
        {
          noise = noise_for_empty_window_;
          ++sparse_window_percent;
        }
        else
        {

          hist_rightmost_bin = bin_count_;

          // do iteration on histogram and find threshold
          for (int i = 0; i < 3; ++i)
          {
            // mean
            hist_mean = 0;
            for (int bin = 0; bin < hist_rightmost_bin; ++bin)
            {
              //std::cout << "V: " << bin << " " << hist_mean << " " << histogram[bin] << " " << elements_in_window << " " << bin_value[bin] << "\n";
              // immediate division is numerically more stable
              hist_mean += histogram[bin] / (double) elements_in_window * bin_value[bin];
            }
            //hist_mean = hist_mean / elements_in_window;

            // stdev
            hist_stdev = 0;
            for (int bin = 0; bin < hist_rightmost_bin; ++bin)
            {
              double tmp(bin_value[bin] - hist_mean);
              hist_stdev += histogram[bin] / (double) elements_in_window * tmp * tmp;
            }
            hist_stdev = std::sqrt(hist_stdev);

            //determine new threshold (i.e. the rightmost bin we consider)
            int estimate = (int) ((hist_mean + hist_stdev * stdev_ - 1) / bin_size + 1);
            //std::cout << "E: " << hist_mean << " " << hist_stdev << " " << stdev_ << " " << bin_size<< " " << estimate << "\n";
            hist_rightmost_bin = std::min(estimate, bin_count_);
          }

          // just avoid division by 0
          noise = std::max(1.0, hist_mean);
        }

        // store result
        stn_estimates_[window_count] = (*window_pos_center).getIntensity() / noise;



        // advance the window center by one datapoint
        ++window_pos_center;
        ++window_count;
        // update progress
        SignalToNoiseEstimator<Container>::setProgress(window_count);

      }   // end while

      SignalToNoiseEstimator<Container>::endProgress();

      sparse_window_percent = sparse_window_percent * 100 / window_count;
      // warn if percentage of sparse windows is above 20%
      if (sparse_window_percent > 20)
      {
        std::cerr << "WARNING in SignalToNoiseEstimatorMeanIterative: "
                  << sparse_window_percent
                  << "% of all windows were sparse. You should consider increasing 'win_len' or increasing 'min_required_elements'"
                  << " You should also check the MaximalIntensity value (or the parameters for its heuristic estimation)"
                  << " If it is too low, then too many high intensity peaks will be discarded, which leads to a sparse window!"
                  << std::endl;
      }

      return;

    }   // end of shiftWindow_

    /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_() override
    {
      max_intensity_         = (double)param_.getValue("max_intensity");
      auto_max_stdev_Factor_ = (double)param_.getValue("auto_max_stdev_factor");
      auto_max_percentile_   = param_.getValue("auto_max_percentile");
      auto_mode_             = param_.getValue("auto_mode");
      win_len_               = (double)param_.getValue("win_len");
      bin_count_             = param_.getValue("bin_count");
      stdev_                 = (double)param_.getValue("stdev_mp");
      min_required_elements_ = param_.getValue("min_required_elements");
      noise_for_empty_window_ = (double)param_.getValue("noise_for_empty_window");
      stn_estimates_.clear();
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
    /// multiplier for the stdev of intensities
    double stdev_;
    /// minimal number of elements a window needs to cover to be used
    int min_required_elements_;
    /// used as noise value for windows which cover less than "min_required_elements_"
    /// use a very high value if you want to get a low S/N result
    double noise_for_empty_window_;




  };

} // namespace OpenMS

