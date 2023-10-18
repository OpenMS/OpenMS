// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <map>

//#define DEBUG_CALIBRATION
namespace OpenMS
{
  /**
   @brief This class implements an external calibration for TOF data using external calibrant spectra.

   The procedure is very similar to the one described in Gobom et al. (Anal Chem. 2002, 74 (15) pp 3915-23).
   The input experiment data need to be flight times. They are converted into m/z-values using the
       calibrant spectra. The calibrant spectra and their expected masses are used to determine the
       quadratic dependence of TOF and m/z values.

       @note The input spectra need to contain flight times.

       @note The peaks must be sorted according to ascending m/z!

     @ingroup SignalProcessing
*/
  class OPENMS_DLLAPI TOFCalibration :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:

    /// Default constructor
    TOFCalibration();

    /// Destructor
    ~TOFCalibration() override;


    /*
        @ brief Apply the external calibration using raw calibrant spectra.

        @exception Exception::UnableToCalibrate is thrown if not enough reference masses are observed.

    */
    void pickAndCalibrate(PeakMap & calib_spectra, PeakMap & exp, std::vector<double> & exp_masses);

    /*
        @ brief Apply the external calibration using picked calibrant spectra.

        @exception Exception::UnableToCalibrate is thrown if not enough reference masses are observed.

    */
    void calibrate(PeakMap & calib_spectra, PeakMap & exp, std::vector<double> & exp_masses);

    /// Non-mutable access to the first calibration constant
    inline const std::vector<double> & getML1s() const {return ml1s_; }
    ///mutable access to the first calibration constant
    inline void setML1s(const std::vector<double> & ml1s)
    {
      ml1s_ = ml1s;
    }

    /// Non-mutable access to the second calibration constant
    inline const std::vector<double> & getML2s() const {return ml2s_; }
    /// mutable access to the second calibration constant
    inline void setML2s(const std::vector<double> & ml2s)
    {
      ml2s_ = ml2s;
    }

    /// Non-mutable access to the third calibration constant
    inline const std::vector<double> & getML3s() const {return ml3s_; }
    /// mutable access to the third calibration constant
    inline void setML3s(const std::vector<double> & ml3s)
    {
      ml3s_ = ml3s;
    }

private:
    ///the calibrant spectra still using flight times instead of m/z-values
    PeakMap calib_peaks_ft_;


    /// the expected calibrant masses
    std::vector<double> exp_masses_;

    /// error in ppm after quadratic fit
    std::map<double, std::vector<double> > errors_;

    /// median errors
    std::vector<double> error_medians_;

    ///
    std::vector<double> calib_masses_;

    ///calibration constants from the instrument needed for the conversion of the calibrant spectra
    std::vector<double> ml1s_;
    std::vector<double> ml2s_;
    std::vector<double> ml3s_;

    /// all coefficients of the quadratic fit
    std::vector<double> coeff_quad_fit_;

    /// mean coefficients
    double a_, b_, c_;


    /// Calculates the coefficients of the quadratic fit used for external calibration.
    void calculateCalibCoeffs_(PeakMap & calib_peaks_ft);


    /// determines the monoisotopic peaks
    void getMonoisotopicPeaks_(PeakMap & calib_peaks, std::vector<std::vector<unsigned int> > & monoiso_peaks);

    /**
             @brief Applies the conversion from TOF to m/z-values to all peaks

             Either a 2-point or a 3-point time of flight conversion can be used, as well as
             different constants for each calibrant spectra or one set for all of them.

             The 2-point equation is mass = ml1/10^12 * (TOF * 1000 - ml2).
             The 3-point equation is time = ml2 + sqrt(10^12/ml1 * mass) +  ml3*mass.

        */
    void applyTOFConversion_(PeakMap & calib_spectra);

    /// determine the monoisotopic masses that have matching expected masses
    void matchMasses_(PeakMap & calib_peaks, 
                      std::vector<std::vector<unsigned int> > & monoiso_peaks, 
                      std::vector<unsigned int> & obs_masses, 
                      std::vector<double> & exp_masses, 
                      unsigned int idx);

    /// Calculate the mass value for a given flight time using the coefficients of the quadratic fit in a specific spectrum.
    inline double mQ_(double ft, unsigned int spec)
    {
      return coeff_quad_fit_[3 * spec] + ft * coeff_quad_fit_[3 * spec + 1] + ft * ft * coeff_quad_fit_[3 * spec + 2];
    }

    /// Calculate the mass value for a given flight time using the averaged coefficients of the quadratic fit.
    inline double mQAv_(double ft) const
    {
      return a_ + ft * b_ + ft * ft * c_;
    }

    /// Calculate the average errors of the reference masses over all scans
    void averageErrors_();

    /// Average the coefficients of the quadratic fit
    void averageCoefficients_();
  };

} // namespace OpenMS

