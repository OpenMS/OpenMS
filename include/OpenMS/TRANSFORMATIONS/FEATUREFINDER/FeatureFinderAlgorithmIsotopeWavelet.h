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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletParallelFor.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/CoupledMarrWavelet.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>
#include <time.h>
				
#ifdef OPENMS_HAS_TBB_H
#include <tbb/task_scheduler_init.h>
#include <tbb/pipeline.h>
#include <tbb/parallel_for.h>
#endif


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
		friend class IsotopeWaveletParallelFor<PeakType, FeatureType>;

	 public:
			
		typedef FeatureFinderAlgorithm<PeakType, FeatureType> Base;

		/** @brief Default Constructor */
		FeatureFinderAlgorithmIsotopeWavelet() 
		{ 
				this->defaults_.setValue ("max_charge", 1, "The maximal charge state to be considered.", false);
				this->defaults_.setMinInt ("max_charge", 1);

				this->defaults_.setValue ("intensity_threshold", 2., "The final threshold t' is build upon the formula: t' = av+t*sd," 
																	"where t is the intensity_threshold, av the average intensity within the wavelet transformed signal" 
																	"and sd the standard deviation of the transform. "
																	"If you set intensity_threshold=-1, t' will be zero.\n");
				
				this->defaults_.setValue ("check_ppm", "false", "Enables/disable a ppm test vs. the averagine model.", true); 
				this->defaults_.setValidStrings("check_ppm",StringList::create("true,false"));

				#ifdef OPENMS_HAS_CUDA
					device_num_=0; cudaGetDeviceCount(&device_num_);
					this->defaults_.setValue ("cuda:use_cuda", 1, "Enables/disables computation on CUDA GPUs.\n"
																		"The number indicates the ID of the GPU to use.'-1' disables CUDA at all.\n"
																		"If you plan to use several GPUs simultaneously, just use a non-negative number for this flag.", false);
					this->defaults_.setMaxInt ("cuda:use_cuda", device_num_); 
					#ifdef OPENMS_HAS_TBB_H
						this->defaults_.setValue ("cuda:use_tbb", "true", "Enables/disables computation on several GPUs via Intel TBBs.", true);
						this->defaults_.setValidStrings("cuda:use_tbb",StringList::create("true,false"));
						this->defaults_.setValue ("cuda:tbb:exclude_id", -1, "Disables computation on GPU device with corresponding ID."
																			"Hence, you are able to exclude your 'normal' graphic card from the computation, while all other GPU devices"
																			"will be in use.\nIf you want to include all devices set this entry to -1." , true);
						this->defaults_.setMinInt ("cuda:tbb:exclude_id", -1);
						this->defaults_.setMaxInt ("cuda:tbb:exclude_id", device_num_);	
					#endif						
				#endif

				this->defaults_.setValue ("rt_votes_cutoff", 5, "A parameter of the sweep line algorithm. It determines the minimum number of"
																	"subsequent scans a pattern must occur to be considered as a feature.", true);
				this->defaults_.setMinInt("rt_votes_cutoff", 0);
				this->defaults_.setValue ("rt_interleave", 2, "A parameter of the sweep line algorithm. It determines the maximum number of"
																	"scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing.\nThere is usually no reason to change the default value.", true);
				this->defaults_.setMinInt ("rt_interleave", 0);
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
				time_t start=time(NULL), end;
				DoubleReal max_mz = this->map_->getMax()[1];
				DoubleReal min_mz = this->map_->getMin()[1];
					
				UInt max_size=0;
				#ifdef OPENMS_HAS_CUDA 
					if (use_cuda_ >= 0) //some preprocessing necessary for the GPU computation
					{
						for (UInt i=0; i<this->map_->size(); ++i)
						{
							max_size = std::max (max_size, (*this->map_)[i].size());
						};
					};
				#endif
			
			//Check for useless RT_votes_cutoff_ parameter
				if (RT_votes_cutoff_ > this->map_->size())
				{
					real_RT_votes_cutoff_ = 0;
				}
				else
				{
					real_RT_votes_cutoff_ = RT_votes_cutoff_;
				};
					
				this->ff_->setLogType (ProgressLogger::CMD);
				progress_counter_ = 0;
				this->ff_->startProgress (0, 2*this->map_->size()*max_charge_, "analyzing spectra");
				
				#ifdef OPENMS_HAS_TBB_H
					if (use_tbb_)
					{
						tbb::task_scheduler_init init (gpu_ids_.size());
						std::vector<IsotopeWaveletTransform<PeakType>*> iwts (gpu_ids_.size()); 
						for (UInt t=0; t<gpu_ids_.size(); ++t)
						{
							iwts[t] = new IsotopeWaveletTransform<PeakType> (min_mz, max_mz, max_charge_, 0.2, max_size);
						};

						static tbb::affinity_partitioner ap;
						//The parallel execution over all available GPU devices
						tbb::parallel_for(tbb::blocked_range<size_t>(0, gpu_ids_.size(), 1), IsotopeWaveletParallelFor<PeakType, FeatureType>(iwts, this), ap);
					
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "Final mapping."; std::cout.flush();
						#endif
						*this->features_ = iwts[gpu_ids_.size()-1]->mapSeeds2Features (*this->map_, max_charge_, real_RT_votes_cutoff_); 
					};
				#else
					if (use_tbb_)
					{
						std::cerr << "Error: You requested computation via TBB, but OpenMS has not been configured for TBB usage." << std::endl;
						std::cerr << "Error: You need to rebuild OpenMS using the configure flag \"--enable-tbb-release\" or \"--enable-tbb-debug\"." << std::endl;
						std::cerr << "Error: Please note that the multithreaded FeatureFinder needs necessarily the CUDA libary, which must be enabled with \"--enable-cuda\"." << std::endl; 
					};
				#endif

				if (!use_tbb_)
				{
					IsotopeWaveletTransform<PeakType> iwt (min_mz, max_mz, max_charge_, 0.2, max_size);
				
					for (UInt i=0; i<this->map_->size(); ++i)
					{			
						const MSSpectrum<PeakType>& c_ref ((*this->map_)[i]);

						if (c_ref.size() <= 1) //unable to do transform anything
						{					
							this->ff_->setProgress (progress_counter_+=2);
							continue;
						};

						if (use_cmarr_ > 0)
						{
							if (iwt.estimateCMarrWidth (c_ref))
							{
								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "Sigma estimation for coupled Marr wavelet successful: " << iwt.getSigma() << std::endl; 	
								#endif
							}
							else
							{
								std::cout << "Note: Sigma estimation for coupled Marr wavelet failed.\n";
								std::cout << "Note: Estimating sigma via the average sampling rate: " << iwt.getSigma() << std::endl; 
							};
						};

						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "Spectrum " << i+1 << " (" << (*this->map_)[i].getRT() << ") of " << this->map_->size() << " ... " ; 
							std::cout.flush();
						#endif
						
						if (use_cuda_ < 0)
						{	
							for (UInt c=0; c<max_charge_; ++c)
							{
								MSSpectrum<PeakType> c_trans (c_ref);
								
								iwt.getTransform (c_trans, c_ref, c);
								
								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::stringstream stream;
									stream << "cpu_" << c_ref.getRT() << "_" << c+1 << ".trans\0"; 
									std::ofstream ofile (stream.str().c_str());
									for (UInt k=0; k < c_ref.size(); ++k)
									{
										ofile << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << c_ref[k].getMZ() << "\t" << c_ref[k].getIntensity() << std::endl;
									};
									ofile.close();
								#endif

								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "transform O.K. ... "; std::cout.flush();
								#endif							
								this->ff_->setProgress (++progress_counter_);

								iwt.identifyCharge (c_trans, c_ref, i, c, intensity_threshold_, check_PPMs_, use_cmarr_);
						
								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "charge recognition O.K. ... "; std::cout.flush();
								#endif
								this->ff_->setProgress (++progress_counter_);
							};
						}	
						else
						{
							#ifdef OPENMS_HAS_CUDA
								cudaSetDevice(use_cuda_);
								typename IsotopeWaveletTransform<PeakType>::TransSpectrum c_trans (&(*this->map_)[i]);
								if (iwt.initializeScanCuda ((*this->map_)[i]) == Constants::CUDA_INIT_SUCCESS)
								{
									for (UInt c=0; c<max_charge_; ++c)
									{	
										iwt.getTransformCuda (c_trans, c);

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::stringstream stream;
											stream << "gpu_" << (*this->map_)[i].getRT() << "_" << c+1 << ".trans\0"; 
											std::ofstream ofile (stream.str().c_str());
											for (UInt k=0; k < c_trans.size(); ++k)
											{
												ofile << c_trans.getMZ(k) << "\t" <<  c_trans.getTransIntensity(k) << "\t" << c_ref[k].getIntensity() << std::endl;
											};
											ofile.close();
										#endif					

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
										#endif
										this->ff_->setProgress (++progress_counter_);

										iwt.identifyChargeCuda (c_trans, i, c, intensity_threshold_, check_PPMs_, use_cmarr_);

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda charge recognition for charge " << c+1 << " O.K. ... "; std::cout.flush();
										#endif					
										this->ff_->setProgress (++progress_counter_);	
									};
									iwt.finalizeScanCuda();
								};
							#else
								std::cerr << "Error: You requested computation on GPU, but OpenMS has not been configured for CUDA usage." << std::endl;
								std::cerr << "Error: You need to rebuild OpenMS using the configure flag \"--enable-cuda\"." << std::endl; 
							#endif
						};
		
						iwt.updateBoxStates(*this->map_, i, RT_interleave_, real_RT_votes_cutoff_);
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "updated box states." << std::endl;
						#endif

						std::cout.flush();
					};

					this->ff_->endProgress();

					//Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
					iwt.updateBoxStates(*this->map_, INT_MAX, RT_interleave_, real_RT_votes_cutoff_);
									
					#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			std::cout << "Final mapping."; std::cout.flush();
					#endif
					*this->features_ = iwt.mapSeeds2Features (*this->map_, max_charge_, real_RT_votes_cutoff_); 
				}
				else
				{
					#ifndef OPENMS_HAS_TBB_H
						std::cerr << "Error: You requested multi-threaded computation via threading building blocks, but OpenMS has not been configured for TBB usage." << std::endl;
						std::cerr << "Error: You need to rebuild OpenMS using the configure flag \"--enable-tbb-debug\" resp. \"--enable-tbb-release\"." << std::endl; 
					#endif
				};
				
				end=time(NULL);

				std::cout << "Running time in seconds: " << difftime(end, start) << std::endl; 
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
			UInt RT_votes_cutoff_, real_RT_votes_cutoff_; ///<The number of subsequent scans a pattern must cover in order to be considered as signal 
		UInt RT_interleave_; ///<The number of scans we allow to be missed within RT_votes_cutoff_
			Int use_cmarr_, use_cuda_;
			bool use_tbb_, check_PPMs_;
			std::vector<UInt> gpu_ids_; ///< A list of all GPU devices that can be used
			
			#ifdef OPENMS_HAS_TBB_H
				tbb::atomic<int> progress_counter_;		
				Int device_num_, gpu_to_exclude_;	
			#else
				Int progress_counter_;
			#endif

		void updateMembers_() 
		{
				max_charge_ = this->param_.getValue ("max_charge");
				intensity_threshold_ = this->param_.getValue ("intensity_threshold");
				RT_votes_cutoff_ = this->param_.getValue ("rt_votes_cutoff");
				RT_interleave_ = this->param_.getValue ("rt_interleave");
			IsotopeWavelet::setMaxCharge(max_charge_);
				use_cmarr_ = this->param_.getValue ("use_cmarr");
				check_PPMs_ = ( (String)(this->param_.getValue("check_ppm"))=="true" );
				#ifdef OPENMS_HAS_CUDA 
					use_cuda_ = this->param_.getValue ("cuda:use_cuda");
					#ifdef OPENMS_HAS_TBB_H 
						use_tbb_ = ( (String)(this->param_.getValue("cuda:use_tbb"))=="true" );
						gpu_to_exclude_ = this->param_.getValue ("cuda:tbb:exclude_id");				
						for (Int i=0; i<device_num_; ++i)
						{
							if (i == gpu_to_exclude_)
							{
								continue;
							};
							gpu_ids_.push_back(i);
						};
					#else
						use_tbb_ = false;
					#endif	
				#else
					use_cuda_ = -1;
				#endif
		}
	};

} //namespace

#endif 
