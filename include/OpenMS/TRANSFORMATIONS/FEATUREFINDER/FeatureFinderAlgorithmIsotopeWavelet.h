// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <iostream>
#include <time.h>
#include <algorithm>
				
#ifdef OPENMS_HAS_TBB
	#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletParallelFor.h>
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
		#ifdef OPENMS_HAS_TBB
			friend class IsotopeWaveletParallelFor<PeakType, FeatureType>;
		#endif

	 public:
			
		typedef FeatureFinderAlgorithm<PeakType, FeatureType> Base;

		/** @brief Default Constructor */
		FeatureFinderAlgorithmIsotopeWavelet() 
		{ 
				this->defaults_.setValue ("max_charge", 3, "The maximal charge state to be considered.");
				this->defaults_.setMinInt ("max_charge", 1);

				this->defaults_.setValue ("intensity_threshold", 2., "The final threshold t' is build upon the formula: t' = av+t*sd," 
																	"where t is the intensity_threshold, av the average intensity within the wavelet transformed signal" 
																	"and sd the standard deviation of the transform. "
																	"If you set intensity_threshold=-1, t' will be zero.\n");
				
				this->defaults_.setValue ("check_ppm", "false", "Enables/disables a ppm test vs. the averagine model, i.e. "
																	"potential peptide masses are checked for plausibility.", StringList::create("advanced")); 
				this->defaults_.setValidStrings("check_ppm",StringList::create("true,false"));


				#if (defined(OPENMS_HAS_CUDA) || defined(OPENMS_HAS_TBB))
					this->defaults_.setValue ("parallel:use_gpus", "-1", "A comma-separated list of IDs corresponding to the GPU devices to use.\n"
																		"'-1' disables parallelization (CUDA/TBB) at all.\n");
				#endif

				this->defaults_.setValue ("sweep_line:rt_votes_cutoff", 5, "Defines the minimum number of "
																	"subsequent scans where a pattern must occur to be considered as a feature.", StringList::create("advanced"));
				this->defaults_.setMinInt("sweep_line:rt_votes_cutoff", 0);
				this->defaults_.setValue ("sweep_line:rt_interleave", 2, "Defines the maximum number of "
							"scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing. There is usually no reason to change the default value.", StringList::create("advanced"));
				this->defaults_.setMinInt ("sweep_line:rt_interleave", 0);
				this->defaultsToParam_();
		}


		/** @brief Destructor. */		
		virtual ~FeatureFinderAlgorithmIsotopeWavelet() 
		{
		}	


		/** @brief The working horse of this class. */
		void run ()
		{
				DoubleReal max_mz = this->map_->getMax()[1];
				DoubleReal min_mz = this->map_->getMin()[1];
					
				Size max_size=0;
				#ifdef OPENMS_HAS_CUDA 
					if (use_cuda_) //some preprocessing necessary for the GPU computation
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
					
				progress_counter_ = 0;
				this->ff_->startProgress (0, 2*this->map_->size()*max_charge_, "analyzing spectra");
			
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				time_t start=time(NULL), end;	
				#endif
				
				#if defined(OPENMS_HAS_TBB) && defined (OPENMS_HAS_CUDA)
					if (use_tbb_)
					{
						UInt num_gpus = this->gpu_ids_.size();
						tbb::task_scheduler_init init (num_gpus);
						std::vector<IsotopeWaveletTransform<PeakType>*> iwts (num_gpus); 
						for (UInt t=0; t<num_gpus; ++t)
						{
							iwts[t] = new IsotopeWaveletTransform<PeakType> (min_mz, max_mz, max_charge_, (UInt)max_size);
						};

						static tbb::affinity_partitioner ap;
						//The parallel execution over all available GPU devices
						tbb::parallel_for(tbb::blocked_range<size_t>(0, num_gpus, 1), IsotopeWaveletParallelFor<PeakType, FeatureType>(iwts, this), ap);
						
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "Merging."; std::cout.flush();
						#endif

						for (UInt t=1; t<num_gpus; ++t)
						{
							iwts[0]->mergeFeatures (iwts[t], RT_interleave_, RT_votes_cutoff_);	
						};

						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "Final mapping."; std::cout.flush();
						#endif
						*this->features_ = iwts[0]->mapSeeds2Features (*this->map_, real_RT_votes_cutoff_); 
						
						for (UInt t=0; t<num_gpus; ++t)
						{
							delete (iwts[t]);							
						};
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
					IsotopeWaveletTransform<PeakType> iwt (min_mz, max_mz, max_charge_, (UInt)max_size);
					#ifdef OPENMS_HAS_CUDA
						if (use_cuda_)
						{
							cudaSetDevice(gpu_ids_[0]);
						};
					#endif
	
					for (UInt i=0; i<this->map_->size(); ++i)
					{			
						const MSSpectrum<PeakType>& c_ref ((*this->map_)[i]);
					
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << ::std::fixed << ::std::setprecision(6) << "Spectrum " << i+1 << " (" << (*this->map_)[i].getRT() << ") of " << this->map_->size() << " ... " ; 
							std::cout.flush();
						#endif
						
						if (c_ref.size() <= 1) //unable to do transform anything
						{					
							#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
								std::cout << "scan empty or consisting of a single data point. Skipping." << std::endl; 
							#endif
							this->ff_->setProgress (progress_counter_+=2);
							continue;
						};

						if (!use_cuda_)
						{	
							iwt.initializeScan ((*this->map_)[i]);
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
										ofile << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << c_ref[k].getIntensity() << std::endl;
									};
									ofile.close();
								#endif

								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "transform O.K. ... "; std::cout.flush();
								#endif							
								this->ff_->setProgress (++progress_counter_);

								iwt.identifyCharge (c_trans, c_ref, i, c, intensity_threshold_, check_PPMs_);
						
								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "charge recognition O.K. ... "; std::cout.flush();
								#endif
								this->ff_->setProgress (++progress_counter_);
							};
						}	
						else
						{
							#ifdef OPENMS_HAS_CUDA
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
												ofile << c_trans.getMZ(k) << "\t" <<  c_trans.getTransIntensity(k) << "\t" << c_trans.getMZ(k) << "\t" << c_ref[k].getIntensity() << std::endl;
											};
											ofile.close();
										#endif					

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
										#endif
										this->ff_->setProgress (++progress_counter_);

										iwt.identifyChargeCuda (c_trans, i, c, intensity_threshold_, check_PPMs_);

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda charge recognition for charge " << c+1 << " O.K." << std::endl;
										#endif					
										this->ff_->setProgress (++progress_counter_);	
									};
									iwt.finalizeScanCuda();
								}
								else
								{
									std::cout << "Warning/Error generated at scan " << i << " (" << (*this->map_)[i].getRT() << ")." << std::endl;
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
					*this->features_ = iwt.mapSeeds2Features (*this->map_, real_RT_votes_cutoff_); 
				}
				else
				{
					#ifndef OPENMS_HAS_TBB
						std::cerr << "Error: You requested multi-threaded computation via threading building blocks, but OpenMS has not been configured for TBB usage." << std::endl;
						std::cerr << "Error: You need to rebuild OpenMS with -DENABLE_TBB=ON." << std::endl; 
					#endif
				};
				
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				end=time(NULL);
				#endif
				
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				std::cout << "Running time in seconds: " << difftime(end, start) << std::endl; 
				#endif
		}

		static const String getProductName()
		{ 
			return ("isotope_wavelet"); 
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
		UInt RT_votes_cutoff_, real_RT_votes_cutoff_, RT_interleave_; ///<The number of subsequent scans a pattern must cover in order to be considered as signal 
		String use_gpus_;
		bool use_tbb_, use_cuda_, check_PPMs_;
		std::vector<UInt> gpu_ids_; ///< A list of all GPU devices that can be used
			
		#if defined(OPENMS_HAS_TBB) && defined(OPENMS_HAS_CUDA)
			tbb::atomic<int> progress_counter_;		
			Int device_num_, gpu_to_exclude_;	
		#else
			Int progress_counter_;
		#endif

		void updateMembers_() 
		{
			max_charge_ = this->param_.getValue ("max_charge");
			intensity_threshold_ = this->param_.getValue ("intensity_threshold");
			RT_votes_cutoff_ = this->param_.getValue ("sweep_line:rt_votes_cutoff");
			RT_interleave_ = this->param_.getValue ("sweep_line:rt_interleave");
			IsotopeWavelet::setMaxCharge(max_charge_);
			check_PPMs_ = ( (String)(this->param_.getValue("check_ppm"))=="true" );
			#if defined(OPENMS_HAS_CUDA) || defined(OPENMS_HAS_TBB)
				use_gpus_ = this->param_.getValue ("parallel:use_gpus");
				std::vector<String> tokens; 
				use_gpus_.split(',', tokens);
				if (tokens.size()==0) tokens.push_back("-1");
					//Attention: updateMembers_ can be called several times!
				gpu_ids_.clear();
				if (tokens[0].trim().toInt() == -1) //no parallelization
				{
					use_cuda_ = false;
					use_tbb_ = false;
					return;
				};
				gpu_ids_.push_back(tokens[0].trim().toInt());
				use_cuda_ = true;
				use_tbb_ = false;
				for (UInt i=1; i<tokens.size(); ++i)
				{
					gpu_ids_.push_back(tokens[i].trim().toInt());
					use_tbb_ = true;
				};
			
			#else
				use_cuda_ = false;	
				use_tbb_ = false;
			#endif
		}
	};

} //namespace

#endif 
