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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMISOTOPEWAVELET_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMISOTOPEWAVELET_H

#ifndef NULL
#define NULL 0
#endif

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>


namespace OpenMS
{
	/** @brief Implements the isotope wavelet feature finder.
	 *
	 * 	The FeatureFinderAlgorithmIsotopeWavelet class has been designed for finding features in 1D or 2D MS data sets using the isotope wavelet.
	 * 	In the case of two dimensional data, the class provides additionally the sweep line algorithm. Please note that in
	 * 	its current implementation the isotope wavelet feature finder is only applicable to raw data (not to picked data). 
	 *
	 * 	@ref FeatureFinderAlgorithmIsotopeWavelet_Parameters are explained on a separate page.
	 *
	 *	@ingroup FeatureFinder */
	template <typename PeakType, typename FeatureType>
	class FeatureFinderAlgorithmIsotopeWavelet : public FeatureFinderAlgorithm<PeakType, FeatureType> 
	{

		public:

			typedef FeatureFinderAlgorithm<PeakType, FeatureType> Base;

			/** @brief Default Constructor */
			FeatureFinderAlgorithmIsotopeWavelet() 
			{ 
				Base::defaults_.setValue ("max_charge", 1, "The maximal charge state to be considered.", false);
				Base::defaults_.setValue ("intensity_threshold", 1, "The final threshold t' is build upon the formula: t' = av+t*sd\n" 
														"where t is the intensity_threshold, av the average intensity within the wavelet transformed signal\n" 
														"and sd the standard deviation of the transform.\n"
														"If you set intensity_threshold=-1, t' will be zero.\n"
														"For single scan analysis (e.g. MALDI peptide fingerprints) you should start with an intensity_threshold\n"
														"around 0..1 and increase if necessary.", false);
				Base::defaults_.setValue ("rt_votes_cutoff", 5, "A parameter of the sweep line algorithm. It determines the minimum number of\n"
														"subsequent scans a pattern must occur to be considered as a feature.", false);
				Base::defaults_.setValue ("rt_interleave", 2, "A parameter of the sweep line algorithm. It determines the maximum number of\n"
														"scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing.", true);
				Base::defaults_.setValue ("recording_mode", 1, "Determines if the spectra have been recorded in positive ion (1) or\n" 
																	"negative ion (-1) mode.", true);
				Base::defaultsToParam_();
			}


			/** @brief Destructor. */		
			virtual ~FeatureFinderAlgorithmIsotopeWavelet() 
			{
			}	


			/** @brief The working horse of this class. */
			void run ()
			{
				DoubleReal max_mz = Base::map_->getMax()[1];
				DoubleReal min_mz = Base::map_->getMin()[1];

				IsotopeWaveletTransform<PeakType> iwt (min_mz, max_mz, max_charge_);
		
				Base::ff_->setLogType (ProgressLogger::CMD);
				Base::ff_->startProgress (0, 3*Base::map_->size(), "analyzing spectra");  

				UInt RT_votes_cutoff = RT_votes_cutoff_;
				//Check for useless RT_votes_cutoff_ parameter
				if (RT_votes_cutoff_ > Base::map_->size())
				{
					RT_votes_cutoff = 0;
				};
				
				for (UInt i=0, j=0; i<Base::map_->size(); ++i)
				{	
					std::vector<MSSpectrum<PeakType> > pwts (max_charge_, Base::map_->at(i));
					std::cout << "Spectrum " << i+1 << " (" << Base::map_->at(i).getRT() << ") of " << Base::map_->size() << "\t" ; 
					std::cout.flush();
				
					iwt.getTransforms (Base::map_->at(i), pwts, max_charge_, mode_);
					Base::ff_->setProgress (++j);

					#ifdef DEBUG_FEATUREFINDER
						std::cout << "transform O.K. ... "; std::cout.flush();
					#endif
					
					iwt.identifyCharges (pwts,  Base::map_->at(i), i, ampl_cutoff_);
					Base::ff_->setProgress (++j);

					#ifdef DEBUG_FEATUREFINDER
						std::cout << "charge recognition O.K. ... "; std::cout.flush();
					#endif

					iwt.updateBoxStates(*Base::map_, i, RT_interleave_, RT_votes_cutoff);
					Base::ff_->setProgress (++j);

					#ifdef DEBUG_FEATUREFINDER
						std::cout << "updated box states." << std::endl;
					#endif

					std::cout.flush();
				};

				Base::ff_->endProgress();
				
				//Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
				iwt.updateBoxStates(*Base::map_, INT_MAX, RT_interleave_, RT_votes_cutoff); 

				#ifdef DEBUG_FEATUREFINDER
					std::cout << "Final mapping."; std::cout.flush();
				#endif
	
				*Base::features_ = iwt.mapSeeds2Features (*Base::map_, max_charge_, RT_votes_cutoff_);

				#ifdef DEBUG_FEATUREFINDER
					std::vector<DoubleReal> error_prone_scans = iwt.getErrorProneScans();
					if (!error_prone_scans.empty())
					{
						std::cerr << "Warning: some of your scans triggered errors while passing the isotope wavelet transform (IWT)." << std::endl;
						std::cerr << "Please remember that the IWT is only suited for MS and not for MS/MS scans. Hence you should always exclude tandem MS signals from the IWT." << std::endl;
						std::cerr << "Another reason might be a very bad resolution of your scan, s.t. the wavelet is unable to adapt its own spacing in a still reasonable manner." << std::endl;
						std::cerr << "The problematic scans are: " << std::endl;
						for (UInt i=0; i<error_prone_scans.size(); ++i)
						{
							std::cerr << error_prone_scans[i] << "\t"; 
						};
						std::cerr << std::endl;
					};
				#endif

			}

			static const String getProductName()
			{ 
				return ("isotope_wavelet_nofit"); 
			}
					
			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmIsotopeWavelet();
			}


		protected:

			/** @brief Internally used data structure for the sweep line algorithm. */
			struct BoxElement
			{			
				DoubleReal mz;
				UInt c; ///<Note, this is not the charge (it is charge-1!!!)
				DoubleReal score;
				DoubleReal intens;
				DoubleReal RT; ///<The elution time (not the scan index)
			};				

			typedef std::map<UInt, BoxElement> Box; ///<Key: RT (index), value: BoxElement
			typedef DPeak<2> Peak2D; 

			UInt max_charge_; ///<The maximal charge state we will consider
			DoubleReal ampl_cutoff_; ///<The only parameter of the isotope wavelet
			UInt RT_votes_cutoff_; ///<The number of subsequent scans a pattern must cover in order to be considered as signal 
			UInt RT_interleave_; ///<The number of scans we allow to be missed within RT_votes_cutoff_
			Int mode_; ///<Negative or positive charged 

			void updateMembers_() 
			{
				max_charge_ = Base::param_.getValue ("max_charge"); 
				ampl_cutoff_ = Base::param_.getValue ("intensity_threshold");
				RT_votes_cutoff_ = Base::param_.getValue ("rt_votes_cutoff");
				RT_interleave_ = Base::param_.getValue ("rt_interleave");
				mode_ = Base::param_.getValue ("recording_mode");
				IsotopeWavelet::setMaxCharge(max_charge_);
			}
	};

} //namespace

#endif 
