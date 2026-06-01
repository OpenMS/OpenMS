// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>
#include <cmath>

namespace OpenMS
{
  class MSExperiment;

  /**
    @brief This class represents the abstract base class of a signal to noise estimator.

    A signal to noise estimator should provide the signal to noise ratio of all raw data points
    in a given interval [first_,last_).
  */

  template <typename Container = MSSpectrum>
  class SignalToNoiseEstimator :
    public DefaultParamHandler, public ProgressLogger
  {
public:

    /** @name Type definitions
     */
    //@{
    typedef typename Container::const_iterator PeakIterator;
    typedef typename PeakIterator::value_type PeakType;


    //@}

    /// Constructor
    inline SignalToNoiseEstimator() :
      DefaultParamHandler("SignalToNoiseEstimator"),
      ProgressLogger()
    {
    }

    /// Copy constructor
    inline SignalToNoiseEstimator(const SignalToNoiseEstimator & source) :
      DefaultParamHandler(source),
      ProgressLogger(source),
      stn_estimates_(source.stn_estimates_)
    {}

    /// Assignment operator
    inline SignalToNoiseEstimator & operator=(const SignalToNoiseEstimator & source)
    {
      if (&source == this) return *this;

      DefaultParamHandler::operator=(source);
      ProgressLogger::operator=(source);
      stn_estimates_ = source.stn_estimates_;
      return *this;
    }

    /// Destructor
    ~SignalToNoiseEstimator() override
    {}

    /// Set the start and endpoint of the raw data interval, for which signal to noise ratios will be estimated immediately
    virtual void init(const Container& c)
    {
      computeSTN_(c);
    }

    ///Return to signal/noise estimate for date point @p index
    ///@note you will get a warning to stderr if more than 20% of the
    ///      noise estimates used sparse windows
    virtual double getSignalToNoise(const Size index) const
    {
      OPENMS_POSTCONDITION(index < stn_estimates_.size(),"SignalToNoiseEstimator estimates beyond container size was requested.");
      return stn_estimates_[index];
    }

protected:

    /**
         * @brief computes the S/N values when init() is called
         *
         * @exception Throws Exception::InvalidValue
         */
    virtual void computeSTN_(const Container& c) = 0;



    /**
      @brief protected struct to store parameters my, sigma for a Gaussian distribution

      Accessors are : mean and variance
    */
    struct GaussianEstimate
    {
      double mean;   ///< mean of estimated Gaussian
      double variance; ///< variance of estimated Gaussian
    };


    /// calculate mean & stdev of intensities of a spectrum
    inline GaussianEstimate estimate_(const PeakIterator & scan_first_, const PeakIterator & scan_last_) const
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
      m = m / size;

      //determine variance
      run = scan_first_;
      while (run != scan_last_)
      {
        double tmp(m - (*run).getIntensity());
        v += tmp * tmp;
        ++run;
      }
      v = v / ((double)size); // divide by n

      GaussianEstimate value = {m, v};
      return value;
    }

    //MEMBERS:

    /// stores the noise estimate for each peak
    std::vector<double> stn_estimates_;
  };

  /// Picks @p n_scans from the given @p ms_level randomly and returns either average intensity at a certain @p percentile.
  /// If no scans with the required level are present, 0.0 is returned
  OPENMS_DLLAPI float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans = 10, const double percentile = 80);

} // namespace OpenMS

