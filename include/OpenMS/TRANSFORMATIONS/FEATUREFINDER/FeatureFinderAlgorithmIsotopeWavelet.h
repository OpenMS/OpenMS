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
// $Maintainer: Rene Hussong $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMISOTOPEWAVELET_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMISOTOPEWAVELET_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/CoupledMarrWavelet.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>
#include <time.h>


namespace OpenMS
{
	/**
		@brief Implements the isotope wavelet feature finder.

		The FeatureFinderAlgorithmIsotopeWavelet class has been designed for finding features in 1D or 2D MS data sets using the isotope wavelet.
		In the case of two dimensional data, the class provides additionally the sweep line algorithm. Please note that in
		its current implementation the isotope wavelet feature finder is only applicable to raw data (not to picked data). 
		
		@htmlinclude OpenMS_FeatureFinderAlgorithmIsotopeWavelet.parameters
		
		@ingroup FeatureFinder
	*/
	template <typename PeakType, typename FeatureType>
	class FeatureFinderAlgorithmIsotopeWavelet : public FeatureFinderAlgorithm<PeakType, FeatureType> 
	{

	 public:

		typedef FeatureFinderAlgorithm<PeakType, FeatureType> Base;

		/** @brief Default Constructor */
		FeatureFinderAlgorithmIsotopeWavelet() 
		{ 
				#ifdef OPENMS_HAS_CUDA
					this->defaults_.setValue ("use_cuda", -1, "Negative, if the computations should be performed on the CPU.\n For evaluation on GPUs," 
																		"enter the corresponding device ID.", false);
				#endif
				this->defaults_.setValue ("max_charge", 1, "The maximal charge state to be considered.", false);
				this->defaults_.setValue ("intensity_threshold", 2., "The final threshold t' is build upon the formula: t' = av+t*sd\n" 
																	"where t is the intensity_threshold, av the average intensity within the wavelet transformed signal\n" 
																	"and sd the standard deviation of the transform.\n"
																	"If you set intensity_threshold=-1, t' will be zero.\n"
																	"For single scan analysis (e.g. MALDI peptide fingerprints) you should start with an intensity_threshold\n"
																	"around 0..1 and increase if necessary.", false);
				this->defaults_.setValue ("rt_votes_cutoff", 5, "A parameter of the sweep line algorithm. It determines the minimum number of\n"
																	"subsequent scans a pattern must occur to be considered as a feature.", true);
				this->defaults_.setValue ("rt_interleave", 2, "A parameter of the sweep line algorithm. It determines the maximum number of\n"
																	"scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing.\n There is usually not reason to change the default value", true);
				this->defaults_.setValue ("use_cmarr", 0, "Experimental, do not enable this feature at the moment!", true);
				this->defaultsToParam_();
		}


		/** @brief Destructor. */		
		virtual ~FeatureFinderAlgorithmIsotopeWavelet() 
		{
		}	


		/** @brief The working horse of this class. */
		void run ()
		{
				clock_t start=clock(), end;
				DoubleReal max_mz = this->map_->getMax()[1];
				DoubleReal min_mz = this->map_->getMin()[1];
					
				UInt max_size=0;
				#ifdef OPENMS_HAS_CUDA 
					if (use_cuda_ >=0) //some preprocessing necessary for the GPU computation
					{
						for (UInt i=0; i<this->map_->size(); ++i)
						{
							max_size = std::max (max_size, (*this->map_)[i].size());
						};
					};
				#endif
				
				IsotopeWaveletTransform<PeakType> iwt (min_mz, max_mz, max_charge_, 0.2, max_size);
	
				this->ff_->setLogType (ProgressLogger::CMD);
				this->ff_->startProgress (0, 3*this->map_->size(), "analyzing spectra");  

			UInt RT_votes_cutoff = RT_votes_cutoff_;
			//Check for useless RT_votes_cutoff_ parameter
				if (RT_votes_cutoff_ > this->map_->size())
			{
				RT_votes_cutoff = 0;
			};
				
				for (UInt i=0, j=0; i<this->map_->size(); ++i)
				{			
					
					if (use_cmarr_ > 0)
					{
						if (iwt.estimateCMarrWidth (this->map_->at(i)))
						{
							std::cout << "Sigma estimation for coupled Marr wavelet successful: " << iwt.getSigma() << std::endl; 
						}
						else
						{
							std::cout << "Sigma estimation for coupled Marr wavelet failed.\n";
							std::cout << "Estimating sigma via the average sampling rate: " << iwt.getSigma() << std::endl; 
						};
					};

					#ifdef OPENMS_DEBUG
						std::cout << "Spectrum " << i+1 << " (" << this->map_->at(i).getRT() << ") of " << this->map_->size() << " ... " ; 
				std::cout.flush();
					#endif
					
					if (this->map_->at(i).size() <= 1) //unable to do transform anything
					{					
						this->ff_->setProgress (j+=3);
						continue;
					};

 
					if (use_cuda_ < 0)
					{	
						std::vector<MSSpectrum<PeakType> > pwts (max_charge_, this->map_->at(i));
						iwt.getTransforms (this->map_->at(i), pwts, max_charge_);
		
						#ifdef DEBUG_FEATUREFINDER
							for (UInt c=0; c<max_charge_; ++c)
							{
								std::stringstream stream;
								stream << "org_" << this->map_->at(i).getRT() << "_" << c+1 << ".trans\0"; 
								std::ofstream ofile (stream.str().c_str());
								for (UInt k=0; k < pwts[c].size(); ++k)
								{
									ofile << pwts[c][k] << std::endl;
								};
								ofile.close();
							};
						#endif
					
						this->ff_->setProgress (++j);

						#ifdef OPENMS_DEBUG
							std::cout << "transform O.K. ... "; std::cout.flush();
						#endif
					
						iwt.identifyCharges (pwts,  this->map_->at(i), i, intensity_threshold_, (use_cmarr_>0) ? true : false);
						this->ff_->setProgress (++j);

						#ifdef OPENMS_DEBUG
							std::cout << "charge recognition O.K. ... "; std::cout.flush();
						#endif
					}	
					else
					{
						#ifdef OPENMS_HAS_CUDA
							MSSpectrum<PeakType> c_trans (this->map_->at(i));
							if (iwt.initializeCudaScan (this->map_->at(i), use_cuda_) == Constants::CUDA_INIT_SUCCESS)
							{
								for (UInt c=0; c<max_charge_; ++c)
								{	
									iwt.getCudaTransforms (c_trans, c);


									#ifdef DEBUG_FEATUREFINDER
										std::stringstream stream;
										stream << "cuda_" << this->map_->at(i).getRT() << "_" << c+1 << ".trans\0"; 
										std::ofstream ofile (stream.str().c_str());
										for (UInt k=0; k < c_trans.size(); ++k)
										{
											ofile << c_trans[k].getMZ() << "\t" <<  c_trans[k].getIntensity() << std::endl;
										};
										ofile.close();
									#endif					

									#ifdef OPENMS_DEBUG
										std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
									#endif
							
										
									iwt.identifyCudaCharges (c_trans, this->map_->at(i), i, c, intensity_threshold_, (use_cmarr_>0) ? true : false);

									#ifdef OPENMS_DEBUG
										std::cout << "cuda charge recognition for charge " << c+1 << " O.K. ... "; std::cout.flush();
									#endif						

								};
								iwt.finalizeCudaScan();
								this->ff_->setProgress (j+=2);
							#else
								std::cerr << "You requested computation on GPU, but OpenMS has not been configured for CUDA usage." << std::endl;
								std::cerr << "You need to rebuild OpenMS using the configure flag \"--enable-cuda\"." << std::endl; 
							#endif
						};
					};
	
					iwt.updateBoxStates(*this->map_, i, RT_interleave_, RT_votes_cutoff);
					this->ff_->setProgress (++j);
				std::cout << "updated box states." << std::endl;
#endif

				std::cout.flush();
			};

				this->ff_->endProgress();

			//Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
				iwt.updateBoxStates(*this->map_, INT_MAX, RT_interleave_, RT_votes_cutoff); 

#ifdef DEBUG_FEATUREFINDER
			std::cout << "Final mapping."; std::cout.flush();
#endif
	
				*this->features_ = iwt.mapSeeds2Features (*this->map_, max_charge_, RT_votes_cutoff_);

				/*std::vector<DoubleReal> error_prone_scans = iwt.getErrorProneScans();
			if (!error_prone_scans.empty())
			{
				std::cerr << "Warning: some of your scans triggered errors while passing the isotope wavelet transform (IWT)." << std::endl;
				std::cerr << "Please remember that the IWT is only suited for MS and not for MS/MS scans. Hence you should always exclude tandem MS signals from the IWT." << std::endl;
				std::cerr << "Another reason might be a very bad resolution of your scan, s.t. the wavelet is unable to adapt its own spacing in a still reasonable manner." << std::endl;
				std::cerr << "The problematic scans are: " << std::endl;
				for (Size i=0; i<error_prone_scans.size(); ++i)
				{
					std::cerr << error_prone_scans[i] << "\t"; 
				};
				std::cerr << std::endl;
				};*/
#endif

				end=clock();

				std::cout << "Running time in seconds: " << (end-start)/(float)(CLOCKS_PER_SEC) << std::endl; 
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

		UInt max_charge_; ///<The maximal charge state we will consider
			DoubleReal intensity_threshold_; ///<The only parameter of the isotope wavelet
		UInt RT_votes_cutoff_; ///<The number of subsequent scans a pattern must cover in order to be considered as signal 
		UInt RT_interleave_; ///<The number of scans we allow to be missed within RT_votes_cutoff_
			Int use_cmarr_, use_cuda_;

		void updateMembers_() 
		{
				max_charge_ = this->param_.getValue ("max_charge");
				intensity_threshold_ = this->param_.getValue ("intensity_threshold");
				RT_votes_cutoff_ = this->param_.getValue ("rt_votes_cutoff");
				RT_interleave_ = this->param_.getValue ("rt_interleave");
			IsotopeWavelet::setMaxCharge(max_charge_);
				use_cmarr_ = this->param_.getValue ("use_cmarr");
				#ifdef OPENMS_HAS_CUDA 
					use_cuda_ = this->param_.getValue ("use_cuda");
				#else
					use_cuda_ = -1;
				#endif
		}
	};

} //namespace

#endif 
