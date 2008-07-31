// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H
#define OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <gsl/gsl_spline.h>

namespace OpenMS
{
	
  /**
     @brief A simple calibration method using linear interpolation of given reference masses.

     This class implements a simle calibration method: given a list of reference masses,
     the relative errors of the peaks in the data are approximated by linear interpolation and
     subtracted from the data. If the input data is raw data peak picking is done first.

	   @ref InternalCalibration_Parameters are explained on a separate page.
	   
	   @ingroup SignalProcessing
  */
  class InternalCalibration 
  	: public DefaultParamHandler, 
  		public ProgressLogger
  {
  public:
    /// Default constructor
    InternalCalibration();

		/// Destructor
    ~InternalCalibration(){}

		/// Copy c'tor
    InternalCalibration(InternalCalibration& obj);

		/// Assignment operator
    InternalCalibration& operator=(const InternalCalibration& obj);


		/**
			 Calibrate a map using given reference masses. Set the flag peak data if you want to calibrate
			 already picked data.

		*/		
    template<typename InputPeakType>
    void calibrate(MSExperiment<InputPeakType>& exp, std::vector<double>& ref_masses,bool peak_data=false);

    /// Non-mutable access to the picked peaks
		inline DoubleReal getWindowLength() const {return window_length_;}
		/// Mutable access to the peaks
		inline void setWindowLength(const DoubleReal window_length) 
		{
			window_length_ = window_length;
			param_.setValue("window_length",window_length);
		}

		/// Non-mutable access to the picked peaks
		inline const MSExperiment<>& getPeaks() const {return exp_peaks_;}
		/// Mutable access to the peaks
		inline void setPeaks(const MSExperiment<>& exp_peaks) {exp_peaks_ = exp_peaks;}

		/// Non-mutable access to the monoisotopic peaks
		inline const std::vector<std::vector<UInt> >& getMonoisotopicPeaks() const {return monoiso_peaks_;}
		/// Mutable access to the calibrant spectra
		inline void setMonoisotopicPeaks(const std::vector<std::vector<UInt> >& monoiso_peaks) {monoiso_peaks_ = monoiso_peaks;}

  protected:

		DoubleReal window_length_;

		MSExperiment<> exp_peaks_;

		std::vector<std::vector<UInt> > monoiso_peaks_;
		
    /// Determines the monoisotopic peaks
    void getMonoisotopicPeaks_();

		// The actual calibration function
		template<typename InputPeakType>
		void calibrate_(MSExperiment<InputPeakType>& exp, std::vector<double>& ref_masses);
		
		void updateMembers_();	

  };// class InternalCalibration


	template<typename InputPeakType>
  void InternalCalibration::calibrate(MSExperiment<InputPeakType>& exp, std::vector<double>& ref_masses,bool peak_data)
  {
#ifdef DEBUG_CALIBRATION
		std::cout.precision(12);
#endif
	
		if(peak_data)
			{
				exp_peaks_ = exp;
			}
		else
			{
				exp_peaks_.clear();
				monoiso_peaks_.clear();
				
				// sort data
				exp.sortSpectra(true);
				
				// pick peaks (only in a certain distance to the reference masses)
				PeakPickerCWT pp;
				pp.setParameters(param_.copy("PeakPicker:",true));
				typename MSExperiment<InputPeakType>::ConstIterator exp_iter = exp.begin();
				typename MSExperiment<InputPeakType>::SpectrumType::ConstIterator spec_iter_l,spec_iter_r; 
				for(;exp_iter != exp.end();++exp_iter)
				{
						MSSpectrum<> spec;
						// pick region around each reference mass
						std::vector<double>::iterator vec_iter = ref_masses.begin();
						for(;vec_iter != ref_masses.end();++vec_iter)
						{
								MSSpectrum<> tmp_spec;		
								// determine region
								spec_iter_l =  (exp_iter->MZBegin(*vec_iter-window_length_));
								// check borders (avoid )
								spec_iter_r =  (exp_iter->MZBegin(*vec_iter+window_length_));
								if((spec_iter_l >= exp_iter->end()) || (spec_iter_r >= exp_iter->end())) continue;

								// pick region
								pp.pick(spec_iter_l,spec_iter_r,tmp_spec);
								typename MSSpectrum<>::Iterator spec_iter = tmp_spec.begin();
								// store
								for(;spec_iter != tmp_spec.end(); ++spec_iter)
								{
										spec.push_back(*spec_iter);
								}
						}
						if(!spec.empty()) exp_peaks_.push_back(spec);
				}
			}
		
		calibrate_(exp,ref_masses);
	}

	template<typename InputPeakType>
  void InternalCalibration::calibrate_(MSExperiment<InputPeakType>& exp, std::vector<double>& ref_masses)
  {
    // get monoisotopic peaks
    getMonoisotopicPeaks_();

    
    size_t num_ref_peaks = ref_masses.size();
    std::vector<double> corr_masses,rel_errors;
    corr_masses.resize(num_ref_peaks,0.);
    rel_errors.resize(num_ref_peaks,0.);
    startProgress(0,monoiso_peaks_.size(),"calibrate spectra");    
    // for each spectrum
    for(size_t spec=0;spec <  monoiso_peaks_.size(); ++spec)
      {
				UInt corr_peaks=0;
				for(size_t peak=0;peak <  monoiso_peaks_[spec].size(); ++peak)
					{
						for(size_t ref_peak=0; ref_peak < num_ref_peaks;++ref_peak)
							{
								if( fabs(exp_peaks_[spec][monoiso_peaks_[spec][peak]].getMZ() - ref_masses[ref_peak]) < 1 )
									{
										corr_masses[ref_peak] = exp_peaks_[spec][monoiso_peaks_[spec][peak]].getMZ();
										++corr_peaks;
										break;
									}
								
							}
					}
				if(corr_peaks < 2)
					{
						std::cout << "spec: "<<spec
											<< " less than 2 reference masses were detected within a reasonable error range\n";
						std::cout << "This spectrum cannot be calibrated!\n";
						continue;
					}
				
				double* x = new double[corr_peaks];
				double* y = new double[corr_peaks];
				UInt p =0;
				// determine rel error in ppm for the two reference masses
				for(size_t ref_peak=0; ref_peak < num_ref_peaks;++ref_peak)
					{
						if(corr_masses[ref_peak] != 0.)
							{
								rel_errors[ref_peak] = (ref_masses[ref_peak]-corr_masses[ref_peak])/corr_masses[ref_peak] * 1e6;
								x[p] =corr_masses[ref_peak];
								y[p] =rel_errors[ref_peak];
								
								++p;
								
							}
					}
				
				
				// linear interpolation
				gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear,corr_peaks);
				gsl_interp_init(interp, x, y, corr_peaks);
				gsl_interp_accel* acc = gsl_interp_accel_alloc();

				
				// use interp to internally calibrate the whole spectrum
				for(unsigned int peak=0;peak <  exp[spec].size(); ++peak)
					{
						exp[spec][peak].setMZ(exp[spec][peak].getMZ() + gsl_interp_eval(interp,x,y,
																																						exp[spec][peak].getMZ(),
																																						acc)/1e6*exp[spec][peak].getMZ());
#ifdef DEBUG_CALIBRATION
						std::cout << exp[spec][peak].getMZ()<< "\t"
											<< exp[spec][peak].getMZ() + gsl_interp_eval(interp,x,y,
																																	 exp[spec][peak].getMZ(),
																																	 acc)/1e6*exp[spec][peak].getMZ()
											<< std::endl;
#endif
					}
				delete[] x;
				delete[] y;
				setProgress(spec);
      }// for(size_t spec=0;spec <  monoiso_peaks.size(); ++spec)
		endProgress();
		
  }// calibrate(MSExperiment<InputPeakType> exp, std::vector<Real> ref_masses)

  
} // namespace OpenMS

#endif // OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H

