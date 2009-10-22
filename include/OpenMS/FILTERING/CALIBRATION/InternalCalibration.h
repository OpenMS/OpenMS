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


#ifndef OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H
#define OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/MzMLFile.h>

 #include <gsl/gsl_fit.h>

namespace OpenMS
{
	
  /**
     @brief A simple calibration method using linear interpolation of given reference masses.

     This class implements a simle calibration method: given a list of reference masses,
     the relative errors of the peaks in the data are approximated by linear interpolation and
     subtracted from the data. If the input data is raw data peak picking is done first.
	
	   @htmlinclude OpenMS_InternalCalibration.parameters

	   @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI InternalCalibration 
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
			 Calibrate a map using given reference masses. The calibration function is calculated for each spectrum
			 separately. If not enough reference masses are found for a spectrum it is left uncalibrated.

		*/		
    template<typename InputPeakType>
    void calibrateMapSpectrumwise(const MSExperiment<InputPeakType>& exp,MSExperiment<InputPeakType>& calibrated_exp, std::vector<DoubleReal>& ref_masses);

		/**
			 Calibrate a map using given reference masses. The calibration function is calculated for the whole map.

		*/		
    template<typename InputPeakType>
    void calibrateMapGlobally(const MSExperiment<InputPeakType>& exp,MSExperiment<InputPeakType>& calibrated_exp, std::vector<DoubleReal>& ref_masses);

		/**
			 Calibrate a map using given identifications. The calibration function is calculated for the whole map.

		*/
    template<typename InputPeakType>
    void calibrateMapGlobally(const MSExperiment<InputPeakType>& exp, MSExperiment<InputPeakType>& calibrated_exp,std::vector<PeptideIdentification>& ref_ids);

    template<typename InputPeakType>
		void calibrateMapList(std::vector<MSExperiment<InputPeakType> >& exp_list,std::vector<MSExperiment<InputPeakType> >& calibrated_exp_list, std::vector<DoubleReal>& ref_masses, std::vector<DoubleReal>& detected_background_masses);


  protected:

		// the actual calibration function
		void makeLinearRegression_(std::vector<DoubleReal>& observed_masses, std::vector<DoubleReal>& theoretical_masses);

		void checkReferenceIds_(std::vector<PeptideIdentification>& pep_ids);

		// here the transformation is stored
		TransformationDescription trafo_;
  };// class InternalCalibration


	template<typename InputPeakType>
  void InternalCalibration::calibrateMapSpectrumwise(const MSExperiment<InputPeakType>& exp, MSExperiment<InputPeakType>& calibrated_exp,std::vector<DoubleReal>& ref_masses)
  {
#ifdef DEBUG_CALIBRATION
		std::cout.precision(writtenDigits<DoubleReal>());
#endif
		if(exp.empty())
			{
				std::cout << "Input is empty."<<std::endl;
				return;
			}
		
		if(exp[0].getType() != SpectrumSettings::PEAKS)
			{
				std::cout << "Attention: this function is assuming peak data."<<std::endl;
			}
		calibrated_exp = exp;
		
    Size num_ref_peaks = ref_masses.size();
    bool use_ppm = (param_.getValue("mz_tolerance_unit") == "ppm" ) ? true : false;
		DoubleReal mz_tol = param_.getValue("mz_tolerance");
    startProgress(0,exp.size(),"calibrate spectra");    
    // for each spectrum
    for(Size spec=0;spec <  exp.size(); ++spec)
      {
				// calibrate only MS1 spectra
				if(exp[spec].getMSLevel() != 1)
					{
						continue;
					}
				
				
				std::vector<DoubleReal> corr_masses,rel_errors,found_ref_masses;
				UInt corr_peaks=0;
				for(Size peak=0;peak <  exp[spec].size(); ++peak)
					{
						for(Size ref_peak=0; ref_peak < num_ref_peaks;++ref_peak)
							{
									if(!use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) <  mz_tol)
									{
										found_ref_masses.push_back(ref_masses[ref_peak]);
										corr_masses.push_back(exp[spec][peak].getMZ());
										++corr_peaks;
										break;
									}
									else if(use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) / ref_masses[ref_peak] * 1e6<  mz_tol)
									{
										found_ref_masses.push_back(ref_masses[ref_peak]);
										corr_masses.push_back(exp[spec][peak].getMZ());
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
				
				// determine rel error in ppm for the two reference masses
				for(Size ref_peak=0; ref_peak < found_ref_masses.size();++ref_peak)
					{
							rel_errors.push_back((found_ref_masses[ref_peak]-corr_masses[ref_peak])/corr_masses[ref_peak] * 1e6);
					}

				makeLinearRegression_(corr_masses,found_ref_masses);
				
				// now calibrate the whole spectrum
				for(unsigned int peak=0;peak <  calibrated_exp[spec].size(); ++peak)
					{
#ifdef DEBUG_CALIBRATION
							std::cout << calibrated_exp[spec][peak].getMZ()<< "\t";
#endif
							DoubleReal mz = calibrated_exp[spec][peak].getMZ();
							trafo_.apply(mz);
							calibrated_exp[spec][peak].setMZ(mz);
#ifdef DEBUG_CALIBRATION
						std::cout	<< calibrated_exp[spec][peak].getMZ()<< std::endl;
#endif

					}
				setProgress(spec);
      }// for(Size spec=0;spec <  exp.size(); ++spec)
		endProgress();
	}

	 
  template<typename InputPeakType>
  void InternalCalibration::calibrateMapGlobally(const MSExperiment<InputPeakType>& exp, MSExperiment<InputPeakType>& calibrated_exp,
																								 std::vector<PeptideIdentification>& ref_ids)
	{
		bool use_ppm = param_.getValue("mz_tolerance_unit") == "ppm" ? true : false;
		DoubleReal mz_tolerance = param_.getValue("mz_tolerance");
		if(exp.empty())
			{
				std::cout << "Input is empty."<<std::endl;
				return;
			}
		
		if(exp[0].getType() != SpectrumSettings::PEAKS)
			{
				std::cout << "Attention: this function is assuming peak data."<<std::endl;
			}
		// check if the ids contain meta information about the peak positions
		checkReferenceIds_(ref_ids);

		std::vector<DoubleReal> theoretical_masses,observed_masses;
		for(Size p_id = 0; p_id < ref_ids.size();++p_id)
			{
				for(Size p_h = 0; p_h < ref_ids[p_id].getHits().size();++p_h)
					{
						Int charge = ref_ids[p_id].getHits()[p_h].getCharge();
						DoubleReal theo_mass = ref_ids[p_id].getHits()[p_h].getSequence().getMonoWeight(Residue::Full,charge)/(DoubleReal)charge;
						// first find corresponding ms1-spectrum
						typename MSExperiment<InputPeakType>::ConstIterator rt_iter = exp.RTBegin(ref_ids[p_id].getMetaValue("RT"));
						while(rt_iter != exp.begin() && rt_iter->getMSLevel() != 1) 
							{
								--rt_iter;
							}
						// now find closest peak
						typename MSSpectrum<InputPeakType>::ConstIterator mz_iter = rt_iter->MZBegin(ref_ids[p_id].getMetaValue("MZ"));
						//					std::cout << mz_iter->getMZ() <<" "<<(DoubleReal)ref_ids[p_id].getMetaValue("MZ")<<"\t";
						DoubleReal dist = (DoubleReal)ref_ids[p_id].getMetaValue("MZ") - mz_iter->getMZ();
						//					std::cout << dist << "\t";
						if((mz_iter +1) != rt_iter->end()
							 && fabs((mz_iter +1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < fabs(dist)
							 && mz_iter != rt_iter->begin()
							 && fabs((mz_iter -1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < fabs((mz_iter +1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ"))) // if mz_iter +1 has smaller dist than mz_iter and mz_iter-1
							{
								if((use_ppm &&
									 fabs((mz_iter +1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) / (DoubleReal)ref_ids[p_id].getMetaValue("MZ") *1e06< mz_tolerance) ||
									 (!use_ppm && fabs((mz_iter+1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < mz_tolerance))
									{
										//		std::cout <<(mz_iter +1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")<<"\t";
										observed_masses.push_back((mz_iter +1)->getMZ());
										theoretical_masses.push_back(theo_mass);
										//									std::cout << (mz_iter +1)->getMZ() << " ~ "<<theo_mass << " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
										//			<< "\tplus 1"<< std::endl;
									}
							}
						else if(mz_iter != rt_iter->begin()
										&& fabs((mz_iter -1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < fabs(dist)) // if mz_iter-1 has smaller dist than mz_iter
							{
								if((use_ppm &&
									 fabs((mz_iter -1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) / (DoubleReal)ref_ids[p_id].getMetaValue("MZ") *1e06< mz_tolerance) ||
									 (!use_ppm && fabs((mz_iter-1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < mz_tolerance))
									{
										//									std::cout <<(mz_iter -1)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")<<"\t";
										observed_masses.push_back((mz_iter -1)->getMZ());
										theoretical_masses.push_back(theo_mass);
										//									std::cout << (mz_iter -1)->getMZ() << " ~ "<<theo_mass << " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
										//			<< "\tminus 1"<< std::endl;
									}
							}
						else
							{
								if((use_ppm &&
										fabs((mz_iter)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) / (DoubleReal)ref_ids[p_id].getMetaValue("MZ") *1e06< mz_tolerance) ||
									 (!use_ppm && fabs((mz_iter)->getMZ() - (DoubleReal)ref_ids[p_id].getMetaValue("MZ")) < mz_tolerance))
									{
										
										observed_masses.push_back(mz_iter->getMZ());
										theoretical_masses.push_back(theo_mass);
// 										std::cout <<"\t"<< mz_iter->getMZ() << " ~ "<<theo_mass<< " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
// 															<< "\tat mz_iter"<< std::endl;
									}
							}
					}
			}

		makeLinearRegression_(observed_masses,theoretical_masses);
		static_cast<ExperimentalSettings&>(calibrated_exp) = exp;
		calibrated_exp.resize(exp.size());

		// for each spectrum
		for(Size spec=0;spec <  exp.size(); ++spec)
      {
				// calibrate only MS1 spectra
				if(exp[spec].getMSLevel() != 1)
					{
						calibrated_exp[spec] = exp[spec];
						continue;
					}
				// copy the spectrum meta data
				calibrated_exp[spec] = exp[spec];

				for(unsigned int peak=0;peak <  exp[spec].size(); ++peak)
					{
#ifdef DEBUG_CALIBRATION
						std::cout << exp[spec][peak].getMZ()<< "\t";
#endif
						DoubleReal mz = exp[spec][peak].getMZ();
						trafo_.apply(mz);
						calibrated_exp[spec][peak].setMZ(mz);
#ifdef DEBUG_CALIBRATION
						std::cout << calibrated_exp[spec][peak].getMZ()<< std::endl;
#endif

					}
      }// for(Size spec=0;spec <  exp.size(); ++spec)
	}


	template<typename InputPeakType>
  void InternalCalibration::calibrateMapGlobally(const MSExperiment<InputPeakType>& exp, MSExperiment<InputPeakType>& calibrated_exp,std::vector<DoubleReal>& ref_masses)
	{
		if(exp.empty())
			{
				std::cout << "Input is empty."<<std::endl;
				return;
			}
			
		if(exp[0].getType() != SpectrumSettings::PEAKS)
			{
				std::cout << "Attention: this function is assuming peak data."<<std::endl;
			}


    Size num_ref_peaks = ref_masses.size();
    bool use_ppm = (param_.getValue("mz_tolerance_unit") == "ppm" ) ? true : false;
		DoubleReal mz_tol = param_.getValue("mz_tolerance");
    startProgress(0,exp.size(),"calibrate spectra");    
		std::vector<DoubleReal> corr_masses,rel_errors,found_ref_masses;
		UInt corr_peaks=0;
    // for each spectrum
    for(Size spec=0;spec <  exp.size(); ++spec)
      {
        // calibrate only MS1 spectra
				if(exp[spec].getMSLevel() != 1) continue;
				for(Size peak=0;peak <  exp[spec].size(); ++peak)
					{
						for(Size ref_peak=0; ref_peak < num_ref_peaks;++ref_peak)
							{
								if(!use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) <  mz_tol)
									{
										found_ref_masses.push_back(ref_masses[ref_peak]);
										corr_masses.push_back(exp[spec][peak].getMZ());
										++corr_peaks;
										break;
									}
								else if(use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) / ref_masses[ref_peak] * 1e6<  mz_tol)
									{
										found_ref_masses.push_back(ref_masses[ref_peak]);
										corr_masses.push_back(exp[spec][peak].getMZ());
										++corr_peaks;
										break;
									}
							}
					}
			}
		if(corr_peaks < 2)
			{
				std::cout << "Less than 2 reference masses were detected within a reasonable error range\n";
				std::cout << "This spectrum cannot be calibrated!\n";
				return;
			}
			
		// calculate the (linear) calibration function
		makeLinearRegression_(corr_masses,found_ref_masses);
		static_cast<ExperimentalSettings&>(calibrated_exp) = exp;
		calibrated_exp.resize(exp.size());
    
		// apply the calibration function to each peak
		for(Size spec=0;spec <  exp.size(); ++spec)
      {
				// calibrate only MS1 spectra
				if(exp[spec].getMSLevel() != 1)
					{
						calibrated_exp[spec] = exp[spec];
						continue;
					}

				// copy the spectrum data
				calibrated_exp[spec] = exp[spec];

				for(unsigned int peak=0;peak <  exp[spec].size(); ++peak)
					{
#ifdef DEBUG_CALIBRATION
							std::cout << exp[spec][peak].getMZ()<< "\t";											
#endif
							DoubleReal mz = exp[spec][peak].getMZ();
							trafo_.apply(mz);
							calibrated_exp[spec][peak].setMZ(mz);

#ifdef DEBUG_CALIBRATION
						std::cout << calibrated_exp[spec][peak].getMZ()	<< std::endl;
#endif

					}
				setProgress(spec);
      }// for(Size spec=0;spec <  exp.size(); ++spec)
		endProgress();
	}


} // namespace OpenMS

#endif // OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H

