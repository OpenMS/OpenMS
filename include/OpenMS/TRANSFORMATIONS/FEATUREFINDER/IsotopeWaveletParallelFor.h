// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETPARALLELFOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETPARALLELFOR_H

#include <tbb/blocked_range.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>

#ifndef NULL
#define NULL 0
#endif


namespace OpenMS
{
	template <typename PeakType, typename FeatureType>
	class FeatureFinderAlgorithmIsotopeWavelet;

	template <typename PeakType, typename FeatureType>
	class IsotopeWaveletParallelFor
	{
		public:

			IsotopeWaveletParallelFor (std::vector<IsotopeWaveletTransform<PeakType>*>& iwts, FeatureFinderAlgorithmIsotopeWavelet<PeakType, FeatureType>* ff)
				: iwts_(iwts), ff_(ff)
			{
			};

			void operator() (const tbb::blocked_range<size_t>& r) const
			{
				for (size_t t=r.begin(); t!=r.end(); ++t) //this will be essentially one iteration
				{
					cudaSetDevice(ff_->gpu_ids_[t]);
					IsotopeWaveletTransform<PeakType>* c_iwt = iwts_[t];

					UInt num_gpus = ff_->gpu_ids_.size();
					UInt block_size = (int)(ff_->map_->size() / num_gpus); UInt additional= ff_->map_->size()-num_gpus*block_size;
					UInt from = t*block_size;
					UInt up_to = (t>=num_gpus-1) ? from+block_size+additional : from+block_size;
					
					for (UInt i=from; i<up_to; ++i)
					{
						const MSSpectrum<PeakType>& c_ref ((*ff_->map_)[i]);						
						if (c_ref.size() <= 1) //unable to transform anything
						{					
							//need to do this atomic
							ff_->ff_->setProgress (ff_->progress_counter_+=2);
							continue;
						};

						bool success=true;
						typename IsotopeWaveletTransform<PeakType>::TransSpectrum* c_trans (NULL); 
						if (!ff_->hr_data_) //LowRes data
						{
							c_trans = new typename IsotopeWaveletTransform<PeakType>::TransSpectrum (&(*ff_->map_)[i]);
							success = c_iwt->initializeScanCuda ((*ff_->map_)[i]) == Constants::CUDA_INIT_SUCCESS;

							if (success)
							{		
								for (UInt c=0; c<ff_->max_charge_; ++c)
								{	
									c_iwt->getTransformCuda (*c_trans, c);

									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::stringstream stream;
										stream << "gpu_lowres_" << ((*ff_->map_)[i]).getRT() << "_" << c+1 << ".trans\0"; 
										std::ofstream ofile (stream.str().c_str());
										for (UInt k=0; k < c_trans->size(); ++k)
										{
											ofile << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k) << "\t" << c_trans->getMZ(k) << "\t" << c_trans->getRefIntensity(k) << std::endl;
										};
										ofile.close();
									#endif					

									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
									#endif
									ff_->ff_->setProgress (++ff_->progress_counter_);

									c_iwt->identifyChargeCuda (*c_trans, i, c, ff_->intensity_threshold_, ff_->check_PPMs_);

									#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
										std::cout << "cuda charge recognition for charge " << c+1 << " O.K." << std::endl;
									#endif					
									ff_->ff_->setProgress (++ff_->progress_counter_);	
								};
								c_iwt->finalizeScanCuda();
							}
							else
							{
								std::cout << "Warning/Error generated at scan " << i << " (" << ((*ff_->map_)[i]).getRT() << ")." << std::endl;
							};
						}
						else //HighRes data
						{
							c_trans = ff_->prepareHRDataCuda (i, 0, c_iwt);
							for (UInt c=0; c<ff_->max_charge_; ++c)
							{	
								//c_trans = ff_->prepareHRDataCuda (i, c, c_iwt);
							
								c_iwt->getTransformCuda (*c_trans, c);

								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::stringstream stream;
									stream << "gpu_highres_" << ((*ff_->map_)[i]).getRT() << "_" << c+1 << ".trans\0"; 
									std::ofstream ofile (stream.str().c_str());
									for (UInt k=0; k < c_trans->size(); ++k)
									{
										ofile << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k) << "\t" << c_trans->getMZ(k) << "\t" << c_trans->getRefIntensity(k) << std::endl;
									};
									ofile.close();
								#endif					

								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
								#endif
								ff_->ff_->setProgress (++ff_->progress_counter_);

								c_iwt->identifyChargeCuda (*c_trans, i, c, ff_->intensity_threshold_, ff_->check_PPMs_);

								#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
									std::cout << "cuda charge recognition for charge " << c+1 << " O.K." << std::endl;
								#endif					
								ff_->ff_->setProgress (++ff_->progress_counter_);
							
								//c_trans->destroy();
								//c_iwt->finalizeScanCuda();
							};								

							c_trans->destroy();
							c_iwt->finalizeScanCuda();
						}
					
						delete (c_trans); c_trans=NULL;
						
						c_iwt->updateBoxStates(*ff_->map_, i, ff_->RT_interleave_, ff_->real_RT_votes_cutoff_, from, up_to-1);
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
							std::cout << "updated box states." << std::endl;
						#endif

						std::cout.flush();
					};
					
					c_iwt->updateBoxStates(*ff_->map_, INT_MAX, ff_->RT_interleave_, ff_->real_RT_votes_cutoff_);
				};
			};

		private:
			
			std::vector<IsotopeWaveletTransform<PeakType>*>& iwts_;
			FeatureFinderAlgorithmIsotopeWavelet<PeakType, FeatureType>* ff_;

	}; //class


}//namespace
#endif
