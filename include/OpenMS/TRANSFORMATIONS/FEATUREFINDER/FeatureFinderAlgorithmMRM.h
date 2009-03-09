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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMMRM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMMRM_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

namespace OpenMS
{
	/** 
		@brief FeatureFinderAlgorithm for MRM experiments.

    @htmlinclude OpenMS_FeatureFinderAlgorithmMRM.parameters

		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithmMRM 
		: public FeatureFinderAlgorithm<PeakType, FeatureType>,
			public FeatureFinderDefs
	{
		public:
			///@name Type definitions
			//@{
			typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType;
			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename SpectrumType::MetaDataArrays MetaDataArrays;
			//@}
			
			using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::ff_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::defaults_;
				
		public:			
			/// default constructor 
			FeatureFinderAlgorithmMRM() 
				: FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
				//Isotopic pattern search parameters
				//defaults_.setValue("isotopic_pattern:charge_low",1,"Lowest charge to search for.");
				//defaults_.setMinInt("isotopic_pattern:charge_low",1);
				
				this->defaultsToParam_();
			}
			
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				//-------------------------------------------------------------------------
				//General initialization
				//-------------------------------------------------------------------------

				ff_->startProgress(0, FeatureFinderAlgorithm<PeakType, FeatureType>::map_->size(), "Finding MRM features.");	
				for (Size i = 0; i != FeatureFinderAlgorithm<PeakType, FeatureType>::map_->size(); ++i)
				{
					ff_->setProgress(i);
				}
				 ff_->endProgress();
				// Split the whole map into traces (== MRM transitions)


				// 
			}			
	
			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmMRM();
			}

			static const String getProductName()
			{
				return "MRM";
			}
	
		protected:

			//Docu in base class
			virtual void updateMembers_()
			{
			}
	};

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMMRM_H
