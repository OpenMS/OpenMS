// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

//#define DINGSBUMS


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
				for (size_t t=r.begin(); t!=r.end(); ++t) //this will be essentially one iteratoin
				{
					cudaSetDevice(t);
					IsotopeWaveletTransform<PeakType>* c_iwt = iwts_[t];

					/*div_t help = div((int)ff_->map_->size(), (int)iwts_.size());
					UInt block_size = help.quot; UInt additional= rint(help.rem*ff_->map_->size());
					UInt from = t*block_size;
					UInt up_to = (t>=iwts_.size()-1) ? from+block_size+additional : from+block_size;*/
					DoubleReal gpu0_fac = 0.145, gpu1_fac=0.855;
					UInt from, up_to;
					if (t==0)
					{
						from = 0;
						up_to = trunc(gpu0_fac*ff_->map_->size());
					}
					else
					{
						from = trunc(gpu0_fac*ff_->map_->size())+1;
						up_to = ff_->map_->size();
					};
					
					for (UInt i=from; i<up_to; ++i)
					{
						const MSSpectrum<PeakType>& c_ref ((*ff_->map_)[i]);						
						if (c_ref.size() <= 1) //unable to do transform anything
						{					
							//need to do this atomic
							ff_->ff_->setProgress (ff_->progress_counter_+=2);
							continue;
						};

						typename IsotopeWaveletTransform<PeakType>::TransSpectrum c_trans (&c_ref);
						if (c_iwt->initializeCudaScan (c_ref) == Constants::CUDA_INIT_SUCCESS)
						{
							for (UInt c=0; c<ff_->max_charge_; ++c)
							{
								c_iwt->getCudaTransforms (c_trans, c);

								#ifdef DINGSBUMS
									std::stringstream stream;
									stream << "gpu_" << (*ff_->map_)[i].getRT() << "_" << c+1 << ".trans\0"; 
									std::ofstream ofile (stream.str().c_str());
									for (UInt k=0; k < c_trans.size(); ++k)
									{
										ofile << c_trans.getMZ(k) << "\t" <<  c_trans.getTransIntensity(k) << "\t" << c_ref[k].getIntensity() << std::endl;
									};
									ofile.close();
								#endif					

								#ifdef DINGSBUMS
									std::cout << "cuda transform for charge " << c+1 << "  O.K. ... "; std::cout.flush();
								#endif
								ff_->ff_->setProgress (++ff_->progress_counter_);

								c_iwt->identifyCudaCharges (c_trans, c_ref, i, c, ff_->intensity_threshold_, (ff_->use_cmarr_>0) ? true : false);

								#ifdef DINGSBUMS
									std::cout << "cuda charge recognition for charge " << c+1 << " O.K. ... "; std::cout.flush();
								#endif					
								ff_->ff_->setProgress (++ff_->progress_counter_);	
							};
							c_iwt->finalizeCudaScan();
						};
						
						c_iwt->updateBoxStates(*ff_->map_, i, ff_->RT_interleave_, ff_->real_RT_votes_cutoff_);
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
