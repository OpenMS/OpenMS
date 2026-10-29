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
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianSettings.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <algorithm> //for std::max_element
#include <cmath>
#include <vector>

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan by using the median (histogram based)

    For each datapoint in the given scan, we collect a range of data points
    around it (param: <i>win_len</i>).  The noise for a datapoint is estimated
    to be the median of the intensities of the current window.  If the number
    of elements in the current window is not sufficient (param:
    <i>MinReqElements</i>), the noise level is set to a default value (param:
    <i>noise_for_empty_window</i>).  The whole computation is histogram based,
    so the user will need to supply a number of bins (param: <i>bin_count</i>),
    which determines the level of error and runtime. The maximal intensity for
    a datapoint to be included in the histogram can be either determined
    automatically (params: <i>AutoMaxIntensity</i>, <i>auto_mode</i>) by two
    different methods or can be set directly by the user (param:
    <i>max_intensity</i>).
    Within the histogram bin that contains the median, the noise value is
    linearly interpolated between that bin's edges (assuming a uniform
    distribution of points within the bin), rather than using the bin's
    center. This keeps the estimate close to continuous as the underlying data
    changes, instead of only ever taking one of <i>bin_count</i> discrete
    values.
    If the (estimated) <i>max_intensity</i> value is too low and the median is
    found to be in the last (&highest) bin, a warning will be given.  In this
    case you should increase <i>max_intensity</i> (and optionally the
    <i>bin_count</i>).

    Changing any of the parameters will invalidate the S/N values (which will invoke a recomputation on the next request).

    @note If more than 20 percent of windows have less than <i>min_required_elements</i> of elements, a warning is issued to <i>OPENMS_LOG_WARN</i> and noise estimates in those windows are set to the constant <i>noise_for_empty_window</i>.
    @note If more than 1 percent of median estimations had to rely on the last(=rightmost) bin (which gives an unreliable result), a warning is issued to <i>OPENMS_LOG_WARN</i>.  In this case you should increase <i>max_intensity</i> (and optionally the <i>bin_count</i>). 
    @note You can disable logging this error by setting <i>write_log_messages</i> and read out the values 


        @htmlinclude OpenMS_SignalToNoiseEstimatorMedian.parameters

    @ingroup SignalProcessing
  */

  template <typename Container = MSSpectrum>
  class SignalToNoiseEstimatorMedian :
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
    inline SignalToNoiseEstimatorMedian()
    {
      //set the name for DefaultParamHandler error messages
      this->setName("SignalToNoiseEstimatorMedian");

      defaults_.setValue("max_intensity", -1, "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode)." \
                                              " Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')!" \
                                              " All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin." \
                                              " If you choose 'max_intensity' too small, the noise estimate might be too small as well. " \
                                              " If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime)." \
                                              " In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.", {"advanced"});
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

      defaults_.setValue("min_required_elements", 10, "minimum number of elements required in a window (otherwise it is considered sparse)");
      defaults_.setMinInt("min_required_elements", 1);

      defaults_.setValue("noise_for_empty_window", std::pow(10.0, 20), "noise value used for sparse windows", {"advanced"});

      defaults_.setValue("write_log_messages", "true", "Write out log messages in case of sparse windows or median in rightmost histogram bin");
      defaults_.setValidStrings("write_log_messages", {"true","false"});

      SignalToNoiseEstimator<Container>::defaultsToParam_();
    }

    /// Copy Constructor
    inline SignalToNoiseEstimatorMedian(const SignalToNoiseEstimatorMedian & source) :
      SignalToNoiseEstimator<Container>(source)
    {
      updateMembers_();
    }

    /** @name Assignment
     */
    //@{
    ///
    inline SignalToNoiseEstimatorMedian & operator=(const SignalToNoiseEstimatorMedian & source)
    {
      if (&source == this) return *this;

      SignalToNoiseEstimator<Container>::operator=(source);
      updateMembers_();
      return *this;
    }

    //@}

    /// Destructor
    ~SignalToNoiseEstimatorMedian() override
    {}

    /// Returns how many percent of the windows were sparse
    double getSparseWindowPercent() const
    {
      return sparse_window_percent_;
    }

    /// Returns the percentage where the median was found in the rightmost bin
    double getHistogramRightmostPercent() const
    {
      return histogram_oob_percent_;
    }

    /// The parameters of an estimation as plain values (see the parameters of the class); the defaults are theirs
    using Settings = SignalToNoiseEstimatorMedianSettings;

    /// What an estimation determined besides the S/N values
    struct Statistics
    {
      double max_intensity = 0;          ///< maximal intensity of the histogram (estimated in the auto modes)
      double sparse_window_percent = 0;  ///< percentage of the windows that were sparse
      double histogram_oob_percent = 0;  ///< percentage of the medians found in the rightmost bin
    };

    /**
      @brief Calculates the signal-to-noise values of all data points of @p c with a sliding window, without an estimator object

      Gives the values that init() and getSignalToNoise() give for an estimator with the same parameters, without
      building and checking a parameter set: for callers that estimate many containers with the same settings.

      @param[in] c Raw data, usually an MSSpectrum
      @param[in] settings The parameters of the estimation
      @param[out] sn The S/N of every data point of @p c
      @param[in] progress Logger to report the progress to, or nullptr
      @return The maximal intensity that was used and the percentages of sparse windows and of medians in the rightmost bin
      @exception Exception::InvalidParameter if a setting is outside the range of its parameter
      @exception Exception::InvalidValue if the maximal intensity is not valid for the auto mode
    */
    static Statistics estimate(const Container& c, const Settings& settings, std::vector<double>& sn,
                               const ProgressLogger* progress = nullptr)
    {
      // the ranges of the parameters, which Settings do not go through
      if (!(settings.win_len >= 1.0) || settings.bin_count < 3 || settings.min_required_elements < 1 ||
          settings.auto_mode < MANUAL || settings.auto_mode > AUTOMAXBYPERCENT || settings.max_intensity < -1 ||
          !(settings.auto_max_stdev_factor >= 0.0 && settings.auto_max_stdev_factor <= 999.0))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "SignalToNoiseEstimatorMedian: a setting is outside the range of its parameter (win_len >= 1, "
                                          "bin_count >= 3, min_required_elements >= 1, auto_mode in [-1, 1], max_intensity >= -1, "
                                          "auto_max_stdev_factor in [0, 999])");
      }

      Statistics stats;
      stats.max_intensity = settings.max_intensity;

      //first element in the scan
      PeakIterator scan_first_ = c.begin();
      //last element in the scan
      PeakIterator scan_last_ = c.end();

      // reset the results
      sn.clear();
      sn.resize(c.size());

      // maximal range of histogram needs to be calculated first
      if (settings.auto_mode == AUTOMAXBYSTDEV)
      {
        // use MEAN+auto_max_intensity_*STDEV as threshold
        GaussianEstimate gauss_global = SignalToNoiseEstimator<Container>::estimate_(scan_first_, scan_last_);
        stats.max_intensity = gauss_global.mean + std::sqrt(gauss_global.variance) * settings.auto_max_stdev_factor;
      }
      else if (settings.auto_mode == AUTOMAXBYPERCENT)
      {
        // get value at "auto_max_percentile"th percentile
        // we use a histogram approach here as well.
        if ((settings.auto_max_percentile < 0) || (settings.auto_max_percentile > 100))
        {
          std::string s = StringUtils::toStr(settings.auto_max_percentile);
          throw Exception::InvalidValue(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "auto_mode is on AUTOMAXBYPERCENT! auto_max_percentile is not in [0,100]. Use setAutoMaxPercentile(<value>) to change it!",
                                        s);
        }

        if (c.empty())
        {
          return stats;
        }
        std::vector<int> histogram_auto(100, 0);

        // find maximum of current scan
        auto maxIt = std::max_element(c.begin(), c.end() ,[](const PeakType& a, const PeakType& b){ return a.getIntensity() < b.getIntensity();});
        typename PeakType::IntensityType maxInt = maxIt->getIntensity();

        double bin_size = maxInt > 0 ? static_cast<double>(maxInt) / 100.0 : 1.0; // no division by zero if no intensity is positive

        // fill histogram
        for(const auto& peak : c)
        {
            ++histogram_auto[(int) std::max(0.0, std::min(peak.getIntensity() / bin_size, 99.0))];
        }

        // add up element counts in histogram until ?th percentile is reached
        int elements_below_percentile = std::max(1, (int) (settings.auto_max_percentile * c.size() / 100));
        int elements_seen = 0;
        int i = -1;

        while (i < 99 && elements_seen < elements_below_percentile)
        {
          ++i;
          elements_seen += histogram_auto[i];
        }

        stats.max_intensity = (((double)i) + 0.5) * bin_size;
      }
      else //if (settings.auto_mode == MANUAL)
      {
        if (stats.max_intensity <= 0)
        {
          std::string s = StringUtils::toStr(stats.max_intensity);
          throw Exception::InvalidValue(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "auto_mode is on MANUAL! max_intensity is <=0. Needs to be positive! Use setMaxIntensity(<value>) or enable auto_mode!",
                                        s);
        }
      }

      if (stats.max_intensity < 0)
      {
        OPENMS_LOG_WARN << "SignalToNoiseEstimatorMedian: the max_intensity value should be positive! " << stats.max_intensity << std::endl;
        return stats;
      }

      PeakIterator window_pos_center  = scan_first_;
      PeakIterator window_pos_borderleft = scan_first_;
      PeakIterator window_pos_borderright = scan_first_;

      double window_half_size = settings.win_len / 2;
      double bin_size = std::max(1.0, stats.max_intensity / settings.bin_count); // at least size of 1 for intensity bins
      int bin_count_minus_1 = settings.bin_count - 1;

      std::vector<int> histogram(settings.bin_count, 0);
      std::vector<double> bin_value(settings.bin_count, 0);
      // calculate average intensity that is represented by a bin
      for (int bin = 0; bin < settings.bin_count; bin++)
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

      ///start progress estimation
      if (progress != nullptr) progress->startProgress(0, c.size(), "noise estimation of data");

      // MAIN LOOP
      while (window_pos_center != scan_last_)
      {

        // erase all elements from histogram that will leave the window on the LEFT side
        while ((*window_pos_borderleft).getPos() < (*window_pos_center).getPos() - window_half_size)
        {
          to_bin = (int) std::max(std::min((*window_pos_borderleft).getIntensity() / bin_size, (double) bin_count_minus_1), 0.0); // clamp before converting to int
          --histogram[to_bin];
          --elements_in_window;
          ++window_pos_borderleft;
        }

        // add all elements to histogram that will enter the window on the RIGHT side
        while ((window_pos_borderright != scan_last_)
              && ((*window_pos_borderright).getPos() <= (*window_pos_center).getPos() + window_half_size))
        {
          //std::cerr << (*window_pos_borderright).getIntensity() << " " << bin_size << " " << bin_count_minus_1 << std::endl;
          to_bin = (int) std::max(std::min((*window_pos_borderright).getIntensity() / bin_size, (double) bin_count_minus_1), 0.0); // clamp before converting to int
          ++histogram[to_bin];
          ++elements_in_window;
          ++window_pos_borderright;
        }

        if (elements_in_window < settings.min_required_elements)
        {
          noise = settings.noise_for_empty_window;
          ++stats.sparse_window_percent;
        }
        else
        {
          // find bin i where ceil[elements_in_window/2] <= sum_c(0..i){ histogram[c] }
          median_bin = -1;
          element_inc_count = 0;
          element_in_window_half = (elements_in_window + 1) / 2;
          while (median_bin < bin_count_minus_1 && element_inc_count < element_in_window_half)
          {
            ++median_bin;
            element_inc_count += histogram[median_bin];
          }

          // increase the error count
          if (median_bin == bin_count_minus_1) {++stats.histogram_oob_percent; }

          // Interpolate within the median bin instead of just reporting its center.
          // Using the bin center alone means the noise estimate can only take one of
          // 'bin_count' discrete values; a tiny, insignificant shift in which side of
          // a bin boundary the cumulative-median crossing falls on then causes a full
          // bin-width jump in the reported noise, even though the underlying data barely
          // changed. Assuming points are uniformly distributed within the bin, we can
          // instead estimate where inside the bin the crossing actually occurs.
          const int elements_in_median_bin = histogram[median_bin];
          double noise_estimate;
          if (elements_in_median_bin > 0)
          {
            const int elements_before_median_bin = element_inc_count - elements_in_median_bin;
            const double median_bin_lower_edge = median_bin * bin_size;
            noise_estimate = median_bin_lower_edge
              + ((double)(element_in_window_half - elements_before_median_bin) / elements_in_median_bin) * bin_size;
          }
          else // only possible if the rightmost bin was hit while empty (already flagged above)
          {
            noise_estimate = bin_value[median_bin];
          }

          // just avoid division by 0
          noise = std::max(1.0, noise_estimate);
        }

        // store result
        sn[window_count] = (*window_pos_center).getIntensity() / noise;


        // advance the window center by one datapoint
        ++window_pos_center;
        ++window_count;
        // update progress
        if (progress != nullptr) progress->setProgress(window_count);

      } // end while

      if (progress != nullptr) progress->endProgress();

      stats.sparse_window_percent = stats.sparse_window_percent * 100 / window_count;
      stats.histogram_oob_percent = stats.histogram_oob_percent * 100 / window_count;

      // warn if percentage of sparse windows is above 20%
      if (stats.sparse_window_percent > 20 && settings.write_log_messages)
      {
        OPENMS_LOG_WARN << "WARNING in SignalToNoiseEstimatorMedian: "
                 << stats.sparse_window_percent
                 << "% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'"
                 << std::endl;
      }

      // warn if percentage of possibly wrong median estimates is above 1%
      if (stats.histogram_oob_percent > 1 && settings.write_log_messages)
      {
        OPENMS_LOG_WARN << "WARNING in SignalToNoiseEstimatorMedian: "
                 << stats.histogram_oob_percent
                 << "% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. "
                 << "You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)"
                 << std::endl;
      }

      return stats;
    }

protected:

    /** Calculate signal-to-noise values for all data points given, by using a sliding window approach

        @param[in] c Raw data, usually an MSSpectrum
        @exception Throws Exception::InvalidValue
    */
    void computeSTN_(const Container& c) override
    {
      Settings settings;
      settings.max_intensity = max_intensity_;
      settings.auto_max_stdev_factor = auto_max_stdev_Factor_;
      settings.auto_max_percentile = auto_max_percentile_;
      settings.auto_mode = auto_mode_;
      settings.win_len = win_len_;
      settings.bin_count = bin_count_;
      settings.min_required_elements = min_required_elements_;
      settings.noise_for_empty_window = noise_for_empty_window_;
      settings.write_log_messages = write_log_messages_;
      // reset the counters, also for an estimation that throws
      sparse_window_percent_ = 0;
      histogram_oob_percent_ = 0;
      const Statistics stats = estimate(c, settings, stn_estimates_, this);
      max_intensity_ = stats.max_intensity;
      sparse_window_percent_ = stats.sparse_window_percent;
      histogram_oob_percent_ = stats.histogram_oob_percent;
    }

    /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_() override
    {
      max_intensity_           = (double)param_.getValue("max_intensity");
      auto_max_stdev_Factor_   = (double)param_.getValue("auto_max_stdev_factor");
      auto_max_percentile_     = param_.getValue("auto_max_percentile");
      auto_mode_               = param_.getValue("auto_mode");
      win_len_                 = (double)param_.getValue("win_len");
      bin_count_               = param_.getValue("bin_count");
      min_required_elements_   = param_.getValue("min_required_elements");
      noise_for_empty_window_  = (double)param_.getValue("noise_for_empty_window");
      write_log_messages_      = (bool)param_.getValue("write_log_messages").toBool();
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
    /// minimal number of elements a window needs to cover to be used
    int min_required_elements_;
    /// used as noise value for windows which cover less than "min_required_elements_"
    /// use a very high value if you want to get a low S/N result
    double noise_for_empty_window_;

    // whether to write out log messages in the case of failure
    bool write_log_messages_;

    // counter for sparse windows
    double sparse_window_percent_;
    // counter for histogram overflow
    double histogram_oob_percent_;


  };

} // namespace OpenMS

