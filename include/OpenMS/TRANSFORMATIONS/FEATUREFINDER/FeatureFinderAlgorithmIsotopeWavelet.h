// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzDataFile.h>
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

		The FeatureFinderAlgorithmIsotopeWavelet class has been designed for finding features in 1D or 2D MS data sets using the isotope wavelet. In the case of two dimensional data, the class provides additionally the sweep line algorithm. Please note that the algorithm implemented here is only marginally related to the algorithm presented in Schulz-Trieglaff et al. (2007, 2008), as no fitting procedure is applied anymore after the wavelet-based seeding step. The wavelet has been designed to extract even very lowly-abundant features (see Hussong et. al (2007, 2009)), usually featuring a very low signal-to-noise ratio. The wavelet in its current implementation is not able to resolve overlapping patterns (see also Hussong et. al (2009)) and slightly shifts masses to the right due to the construction of the wavelet.   
		
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

				this->defaults_.setValue ("intensity_threshold", -1., "The final threshold t' is build upon the formula: t' = av+t*sd, " 
																	"where t is the intensity_threshold, av the average intensity within the wavelet transformed signal " 
																	"and sd the standard deviation of the transform. "
																	"If you set intensity_threshold=-1, t' will be zero.\n"
																	"As the 'optimal' value for this parameter is highly data dependent, we would recommend to start "
																	"with -1, which will also extract features with very low signal-to-noise ratio. Subsequently, one "
																	"might increase the threshold to find an optimized trade-off between false positives and true positives. "
																	"Depending on the dynamic range of your spectra, suitable value ranges include: -1, [0:10], and if your data "
																	"features even very high intensity values, t can also adopt values up to around 30. "
																	"Please note that this parameter is not of an integer type, s.t. you can also use t:=0.1, e.g.");
				this->defaults_.setValue ("intensity_type", "ref", "Determines the intensity type returned for the identified features. 'ref' (default) returns the sum of the intensities of each isotopic peak within an isotope pattern. 'trans' refers to the intensity of the monoisotopic peak within the wavelet transform. 'corrected' refers also to the transformed intensity with an attempt to remove the effects of the convolution. While the latter ones might be preferable for qualitative analyses, 'ref' might be the best option to obtain quantitative results. Please note that intensity values might be spoiled (in particular for the option 'ref'), as soon as patterns overlap (see also the explanations given in the class documentation of FeatureFinderAlgorihtmIsotopeWavelet).", StringList::create("advanced"));
				this->defaults_.setValidStrings("intensity_type",StringList::create("ref,trans,corrected"));
	
				this->defaults_.setValue ("check_ppm", "false", "Enables/disables a ppm test vs. the averagine model, i.e. "
																	"potential peptide masses are checked for plausibility. In addition, a heuristic correcting potential mass shifts induced by the wavelet is applied.", StringList::create("advanced")); 
				this->defaults_.setValidStrings("check_ppm",StringList::create("true,false"));
				
				this->defaults_.setValue ("hr_data", "false", "Must be true in case of high-resolution data, i.e. "
							"for spectra featuring large m/z-gaps (present in FTICR and Orbitrap data, e.g.). Please check " 
							"a single MS scan out of your recording, if you are unsure.");
				this->defaults_.setValidStrings("hr_data",StringList::create("true,false"));

				#if (defined(OPENMS_HAS_CUDA) || defined(OPENMS_HAS_TBB))
					this->defaults_.setValue ("parallel:use_gpus", "-1", "A comma-separated list of IDs corresponding to the GPU devices to use.\n"
																		"'-1' disables parallelization (CUDA/TBB) at all.\n");
				#endif

				this->defaults_.setValue ("sweep_line:rt_votes_cutoff", 5, "Defines the minimum number of "
																	"subsequent scans where a pattern must occur to be considered as a feature.", StringList::create("advanced"));
				this->defaults_.setMinInt("sweep_line:rt_votes_cutoff", 0);
				this->defaults_.setValue ("sweep_line:rt_interleave", 1, "Defines the maximum number of "
							"scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing. There is usually no reason to change the default value.", StringList::create("advanced"));
				this->defaults_.setMinInt ("sweep_line:rt_interleave", 0);

				this->defaultsToParam_();
		}


		/** @brief Destructor. */		
		virtual ~FeatureFinderAlgorithmIsotopeWavelet() 
		{
		}	
	
		typename IsotopeWaveletTransform<PeakType>::TransSpectrum* prepareHRDataCuda (const UInt i, IsotopeWaveletTransform<PeakType>* iwt)
		{
			MSSpectrum<PeakType>* new_spec (createHRData(i));
			typename IsotopeWaveletTransform<PeakType>::TransSpectrum* c_trans = new typename IsotopeWaveletTransform<PeakType>::TransSpectrum (new_spec);
			iwt->initializeScanCuda (*new_spec);	

			return (c_trans);
		}

		MSSpectrum<PeakType>* createHRData (const UInt i)
		{
			MSSpectrum<PeakType> spec ((*this->map_)[i]);
			
			const MSSpectrum<PeakType>& specr ((*this->map_)[i]);
			
			for (UInt j=0; j<spec.size()-1; ++j)
			{
				spec[j].setMZ(-1*(specr[j+1].getMZ()-specr[j].getMZ()));
				spec[j].setIntensity((specr[j].getIntensity()+specr[j+1].getIntensity()));
			}
			spec[spec.size()-1].setMZ(-1); spec[spec.size()-1].setIntensity(-1);

			ConstRefVector<MSSpectrum<PeakType> > c_sorted_spec (spec.begin(), spec.end());	
			//Sort in ascending order according to the intensities present in the transform 	
			c_sorted_spec.sortByPosition();
		
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				std::ofstream ofilex ("spacings.trans");
				for (UInt j=0; j<spec.size()-1; ++j)
				{
					ofilex << ::std::setprecision(12) << std::fixed << spec[j].getMZ() << "\t" << spec[j].getIntensity() << std::endl; 
				};
				ofilex.close();
			#endif

			UInt pos=0;
			while (c_sorted_spec[pos].getIntensity() <= 0)
			{
				if (++pos >= c_sorted_spec.size())
				{
					std::cout << "Detected empty scan or a scan that cannot be interpolated with zeros in HR mode. " << std::endl;
					std::cout << "Please check scan # " << i << " of your data set." << std::endl;
					exit (-1);
				}
			};
			DoubleReal bound=-1*c_sorted_spec[pos].getMZ();
		
			if (bound > (1./max_charge_)/2.) 
			//that might be case for simulated spectra, 
			//which might show a very artificial spacing
			{
				bound = (1./max_charge_)/2./4.;
			};
	
			MSSpectrum<PeakType>* new_spec = new MSSpectrum<PeakType>;
			new_spec->reserve(200000); 
			new_spec->setRT(((*this->map_)[i]).getRT());
			PeakType p; p.setMZ(specr[0].getMZ()); p.setIntensity(specr[0].getIntensity());
			new_spec->push_back(p);

			UInt count;
			for (UInt j=0; j<spec.size()-1; ++j)
			{
				count=0;
				while (-spec[j].getMZ()-count*bound > bound)
				{
					++count;
					p.setMZ(specr[j].getMZ()+count*bound); p.setIntensity(0);
					new_spec->push_back(p);
				}
				p.setMZ(specr[j+1].getMZ()); p.setIntensity(specr[j+1].getIntensity());
				new_spec->push_back(p);
			};

			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				std::ofstream ofiley ("new_spec.trans");
				for (UInt j=0; j<new_spec->size(); ++j)
				{
					ofiley << ::std::setprecision(12) << std::fixed << (*new_spec)[j].getMZ() << "\t" << (*new_spec)[j].getIntensity() << std::endl; 
				};
				ofiley.close();
			#endif

			return (new_spec);
		};



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
					
				this->ff_->setLogType (ProgressLogger::CMD);
				progress_counter_ = 0;
				this->ff_->startProgress (0, 2*this->map_->size()*max_charge_, "analyzing spectra");
			
				#if defined(OPENMS_HAS_TBB) && defined (OPENMS_HAS_CUDA)
					if (use_tbb_)
					{
						UInt num_gpus = this->gpu_ids_.size();
						tbb::task_scheduler_init init (num_gpus);
						std::vector<IsotopeWaveletTransform<PeakType>*> iwts (num_gpus); 
						
						for (UInt t=0; t<num_gpus; ++t)
						{
							iwts[t] = new IsotopeWaveletTransform<PeakType> (min_mz, max_mz, max_charge_, max_size, true, hr_data_, intensity_type_);
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
					IsotopeWaveletTransform<PeakType>* iwt = new IsotopeWaveletTransform<PeakType> (min_mz, max_mz, max_charge_, max_size, use_cuda_, hr_data_, intensity_type_);
					#ifdef OPENMS_HAS_CUDA
						if (use_cuda_)
						{
							cudaSetDevice(gpu_ids_[0]);
							std::cout << "Using device with ID: " << gpu_ids_[0] << std::endl;
							cudaDeviceProp props;
							cudaGetDeviceProperties (&props, gpu_ids_[0]);
							std::cout << "This device is named: " << props.name << std::endl; 
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
							if (!hr_data_) //Lowres data
							{
								iwt->initializeScan ((*this->map_)[i]);
								for (UInt c=0; c<max_charge_; ++c)
								{
									MSSpectrum<PeakType> c_trans (c_ref);
									
									iwt->getTransform (c_trans, c_ref, c);
									
									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::stringstream stream;
										stream << "cpu_lowres_" << c_ref.getRT() << "_" << c+1 << ".trans\0"; 
										std::ofstream ofile (stream.str().c_str());
										for (UInt k=0; k < c_ref.size(); ++k)
										{
											ofile << ::std::setprecision(8) << std::fixed << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << c_ref[k].getIntensity() << std::endl;
										};
										ofile.close();
									#endif

									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "transform O.K. ... "; std::cout.flush();
									#endif							
									this->ff_->setProgress (++progress_counter_);

									iwt->identifyCharge (c_trans, c_ref, i, c, intensity_threshold_, check_PPMs_);
							
									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "charge recognition O.K. ... "; std::cout.flush();
									#endif
									this->ff_->setProgress (++progress_counter_);
								};
							}
							else //Highres data
							{
								MSSpectrum<PeakType>* new_spec (NULL);
								for (UInt c=0; c<max_charge_; ++c)
								{
									new_spec = createHRData (i);
									iwt->initializeScan (*new_spec, c);
									MSSpectrum<PeakType> c_trans (*new_spec);
									
									iwt->getTransformHighRes (c_trans, *new_spec, c);
									
									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::stringstream stream;
										stream << "cpu_highres_" << new_spec->getRT() << "_" << c+1 << ".trans\0"; 
										std::ofstream ofile (stream.str().c_str());
										for (UInt k=0; k < new_spec->size(); ++k)
										{
											ofile << ::std::setprecision(8) << std::fixed << c_trans[k].getMZ() << "\t" << c_trans[k].getIntensity() << "\t" << (*new_spec)[k].getIntensity() << std::endl;
										};
										ofile.close();
									#endif

									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "transform O.K. ... "; std::cout.flush();
									#endif							
									this->ff_->setProgress (++progress_counter_);

									iwt->identifyCharge (c_trans, *new_spec, i, c, intensity_threshold_, check_PPMs_);
							
									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "charge recognition O.K. ... "; std::cout.flush();
									#endif
									this->ff_->setProgress (++progress_counter_);
									
									delete (new_spec); new_spec=NULL;
								};
							};
						}	
						else //use CUDA but not TBB
						{
							#ifdef OPENMS_HAS_CUDA
								bool success=true;
								typename IsotopeWaveletTransform<PeakType>::TransSpectrum* c_trans (NULL); MSSpectrum<PeakType>* new_spec (NULL);
								if (!hr_data_) //LowRes data
								{
									c_trans = new typename IsotopeWaveletTransform<PeakType>::TransSpectrum (&(*this->map_)[i]);
									success = iwt->initializeScanCuda ((*this->map_)[i]) == Constants::CUDA_INIT_SUCCESS;

									if (success)
									{		
										for (UInt c=0; c<max_charge_; ++c)
										{	
											iwt->getTransformCuda (*c_trans, c);

											#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
												std::stringstream stream;
												stream << "gpu_lowres_" << ((*this->map_)[i]).getRT() << "_" << c+1 << ".trans\0"; 
												std::ofstream ofile (stream.str().c_str());
												for (UInt k=0; k < c_trans->size(); ++k)
												{
													ofile << ::std::setprecision(8) << std::fixed << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k) << "\t" << c_trans->getRefIntensity(k) << std::endl;
												};
												ofile.close();
											#endif					

											#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
												std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
											#endif
											this->ff_->setProgress (++progress_counter_);

											iwt->identifyChargeCuda (*c_trans, i, c, intensity_threshold_, check_PPMs_);

											#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
												std::cout << "cuda charge recognition for charge " << c+1 << " O.K." << std::endl;
											#endif					
											this->ff_->setProgress (++progress_counter_);	
										};
										iwt->finalizeScanCuda();
									}
									else
									{
										std::cout << "Warning/Error generated at scan " << i << " (" << ((*this->map_)[i]).getRT() << ")." << std::endl;
									};
								}
								else //HighRes data
								{	
									c_trans = prepareHRDataCuda (i, iwt);
									for (UInt c=0; c<max_charge_; ++c)
									{	
										iwt->getTransformCuda (*c_trans, c);

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::stringstream stream;
											stream << "gpu_highres_" << ((*this->map_)[i]).getRT() << "_" << c+1 << ".trans\0"; 
											std::ofstream ofile (stream.str().c_str());
											for (UInt k=0; k < c_trans->size(); ++k)
											{
												ofile << ::std::setprecision(8) << std::fixed << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k)  << "\t" << c_trans->getRefIntensity(k) << std::endl;
											};
											ofile.close();
										#endif					

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
										#endif
										this->ff_->setProgress (++progress_counter_);

										iwt->identifyChargeCuda (*c_trans, i, c, intensity_threshold_, check_PPMs_);

										#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
											std::cout << "cuda charge recognition for charge " << c+1 << " O.K." << std::endl;
										#endif					
										this->ff_->setProgress (++progress_counter_);
									};									
									c_trans->destroy();
									iwt->finalizeScanCuda();
								};
							
								delete (new_spec); new_spec=NULL;
								delete (c_trans); c_trans=NULL;

							#else
								std::cerr << "Error: You requested computation on GPU, but OpenMS has not been configured for CUDA usage." << std::endl;
								std::cerr << "Error: You need to rebuild OpenMS using the configure flag \"--enable-cuda\"." << std::endl; 
							#endif
						};
		
						iwt->updateBoxStates(*this->map_, i, RT_interleave_, real_RT_votes_cutoff_);
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "updated box states." << std::endl;
						#endif

						std::cout.flush();
					};

					this->ff_->endProgress();

					//Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
					iwt->updateBoxStates(*this->map_, INT_MAX, RT_interleave_, real_RT_votes_cutoff_);
									
					#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
						std::cout << "Final mapping."; std::cout.flush();
					#endif
					*this->features_ = iwt->mapSeeds2Features (*this->map_, real_RT_votes_cutoff_); 

					delete (iwt);
				}
				else
				{
					#ifndef OPENMS_HAS_TBB
						std::cerr << "Error: You requested multi-threaded computation via threading building blocks, but OpenMS has not been configured for TBB usage." << std::endl;
						std::cerr << "Error: You need to rebuild OpenMS with -DENABLE_TBB=ON." << std::endl; 
					#endif
				};
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
		String use_gpus_, intensity_type_;
		bool use_tbb_, use_cuda_, check_PPMs_, hr_data_;
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
			hr_data_ = ( (String)(this->param_.getValue("hr_data"))=="true" );
			intensity_type_ = ( (String)(this->param_.getValue("intensity_type")) );
			#if defined(OPENMS_HAS_CUDA) || defined(OPENMS_HAS_TBB)
				use_gpus_ = this->param_.getValue ("parallel:use_gpus");
				std::vector<String> tokens; 
				if (!use_gpus_.split(',', tokens))	
				{
					tokens.push_back(use_gpus_);
				};			
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
					if (tokens[i].trim().toInt() == (Int) gpu_ids_[i-1])
					{
						continue;
					};
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
