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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_CALIBRATION_TOFCALIBRATION_H
#define OPENMS_FILTERING_CALIBRATION_TOFCALIBRATION_H


#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>

//#define DEBUG_CALIBRATION
namespace OpenMS
{
	/**
     @brief This class implements an external calibration for TOF data using external calibrant spectra.
     
     The procedure is very similar to the one described in Gobom et al. (Anal Chem. 2002, 74 (15) pp 3915-23).
     The input experiment data need to be flight times. They are converted into m/z-values using the
		 calibrant spectra. The calibrant spectra and their expected masses are used to determine the
		 quadratic dependance of tof and m/z values.

		 @note The input spectra need to contain flight times.

		 @note The peaks must be sorted according to ascending m/z!

	   @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI TOFCalibration 
    : public DefaultParamHandler,
    	public ProgressLogger
  {
  public:
    
    /// Default constructor
    TOFCalibration();
    
    /// Destructor
    ~TOFCalibration();
    

		/*
			@ brief Apply the external calibration using raw calibrant spectra.
			
			@exception Exception::UnableToCalibrate is thrown if not enough reference masses are observed.
			
		*/
    template<typename PeakType>
    void pickAndCalibrate(MSExperiment<Peak1D>& calib_spectra,MSExperiment<PeakType >& exp, std::vector<double>& exp_masses);
		
		/*
			@ brief Apply the external calibration using picked calibrant spectra.
			
			@exception Exception::UnableToCalibrate is thrown if not enough reference masses are observed.

		*/
		template<typename PeakType>
    void calibrate(MSExperiment<Peak1D>& calib_spectra,MSExperiment<PeakType >& exp, std::vector<double>& exp_masses);

		/// Non-mutable access to the first calibration constant 
		inline const std::vector<double>& getML1s() const {return ml1s_;}
    ///mutable access to the first calibration constant
    inline void setML1s(const std::vector<double>& ml1s) 
    {
      ml1s_ = ml1s;
    }
		
		/// Non-mutable access to the second calibration constant
		inline const std::vector<double>& getML2s() const {return ml2s_;}		
    /// mutable access to the second calibration constant
    inline void setML2s(const std::vector<double>& ml2s) 
    {
      ml2s_ = ml2s;
    }
		
		/// Non-mutable access to the third calibration constant
		inline const std::vector<double>& getML3s() const {return ml3s_;}
    /// mutable access to the third calibration constant
    inline void setML3s(const std::vector<double>& ml3s) 
    {
      ml3s_ = ml3s;
    }
		
  private:
		///the calibrant spectra still using flight times instead of m/z-values
		MSExperiment<> calib_peaks_ft_;

		
    /// the expected calibrant masses
    std::vector<double> exp_masses_;

    /// error in ppm after quadratic fit
    std::map<double,std::vector<double> > errors_;

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
    double a_,b_,c_;
		

    gsl_interp_accel* acc_;

    gsl_spline* spline_;

    /// Calculates the coefficients of the quadratic fit used for external calibration.
    void calculateCalibCoeffs_(MSExperiment<>& calib_peaks_ft) ;

		
    /// determines the monoisotopic peaks
    void getMonoisotopicPeaks_(MSExperiment<>& calib_peaks, std::vector<std::vector<unsigned int> >& monoiso_peaks);

    /**
			 @brief Applies the conversion from TOF to m/z-values to all peaks

			 Either a 2-point or a 3-point time of flight conversion can be used, as well as
			 different constants for each calibrant spectra or one set for all of them.

			 The 2-point equation is mass = ml1/10^12 * (tof * 1000 - ml2).
			 The 3-point equation is time = ml2 + sqrt(10^12/ml1 * mass) +  ml3*mass.

		*/
    void applyTOFConversion_(MSExperiment<>& calib_spectra);
    
    /// determine the monoisotopic masses that have matching expected masses
    void matchMasses_(MSExperiment<>& calib_peaks,std::vector<std::vector<unsigned int> >& monoiso_peaks, std::vector<unsigned int>& obs_masses,std::vector<double>& exp_masses,unsigned int idx);
		
    /// Calculate the mass value for a given flight time using the coefficients of the quadratic fit in a specific spectrum.
    inline double mQ_(double ft, unsigned int spec)
    {
      return coeff_quad_fit_[3*spec] + ft*coeff_quad_fit_[3*spec+1] + ft*ft*coeff_quad_fit_[3*spec+2]; 
    }
		
    /// Calculate the mass value for a given flight time using the averaged coefficients of the quadratic fit.		
    inline double mQAv_(double ft)
    {
      return a_ + ft*b_ + ft*ft*c_; 
    }
	
    /// Calculate the average errors of the reference masses over all scans
    void averageErrors_();

    /// Average the coefficients of the quadratic fit
    void averageCoefficients_();
  };

	template<typename PeakType>
  void TOFCalibration::pickAndCalibrate(MSExperiment<Peak1D>& calib_spectra,MSExperiment<PeakType >& exp, std::vector<double>& exp_masses)
	{
		MSExperiment<Peak1D> p_calib_spectra;
		
		// pick peaks
		PeakPickerCWT pp;
		pp.setParameters(param_.copy("PeakPicker:",true));
		pp.pickExperiment(calib_spectra,p_calib_spectra);
		
		//calibrate
		calibrate(p_calib_spectra,exp,exp_masses);
	}
	
  template<typename PeakType>
  void TOFCalibration::calibrate(MSExperiment<Peak1D>& calib_spectra,MSExperiment<PeakType >& exp, std::vector<double>& exp_masses)
  {
		exp_masses_ = exp_masses;
		calculateCalibCoeffs_(calib_spectra);
		double m;
    for(unsigned int spec=0;spec <  exp.size(); ++spec)
    {
			for(unsigned int peak=0;peak <  exp[spec].size(); ++peak)
			{
				m = mQAv_(exp[spec][peak].getMZ());
				exp[spec][peak].setPos(m - gsl_spline_eval(spline_,m,acc_));
			}
    }
  }


	
} // namespace OpenMS

#endif // OPENMS_FILTERING_CALIBRATION_TOFCALIBRATION_H

