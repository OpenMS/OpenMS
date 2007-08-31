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
// $Maintainer: Marc Sturm, Marcel Grunert$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H

#include<OpenMS/KERNEL/MSExperiment.h>
#include<OpenMS/KERNEL/FeatureMap.h>
#include<OpenMS/CONCEPT/FactoryProduct.h>
#include<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

namespace OpenMS
{
  /** 
  	@brief Abstract base class for FeatureFinder algorithms 
   
    @ingroup FeatureFinder
  */
	template<class PeakType, class FeatureType>
  class FeatureFinderAlgorithm
		: public FactoryProduct
  {
	  public:	  	
	   	/// Type of the calling FeatureFinder
			typedef FeatureFinder<PeakType,FeatureType> FeatureFinderType;
			///Index to peak consisting of two UInts (scan index / peak index)	
			typedef typename FeatureFinderType::IDX IDX;
			///Index to peak consisting of two UInts (scan index / peak index) and charge information	
			typedef typename FeatureFinderType::ChargedIndexSet ChargedIndexSet;
			///A set of peak indices
			typedef typename FeatureFinderType::IndexSet IndexSet;
			/// Input map type
			typedef typename FeatureFinderType::MapType MapType;
			/// Coordinate/Position type of peaks
			typedef typename FeatureFinderType::CoordinateType CoordinateType;
			/// Intensity type of peaks
			typedef typename FeatureFinderType::IntensityType IntensityType;
			/// Output feature type
			typedef typename FeatureFinderType::FeatureMapType FeatureMapType;
		 	/// Flags that indicate if a peak is alread used in a feature
			enum Flag { UNUSED, USED };
			///Exception that is thrown if a method a invalid IDX is given
			typedef typename FeatureFinderType::NoSuccessor NoSuccessor;
			
			/// default constructor 
	    FeatureFinderAlgorithm()
				: FactoryProduct("FeatureFinderAlgorithm"),
					ff_(0)
			{
			};
	
	    /// destructor 
	    virtual ~FeatureFinderAlgorithm()
			{
			};
	
	    /// register all derived classes here 
	    static void registerChildren();

			/// Main method for actual FeatureFinder
			virtual void run()=0;
			
			/// Sets a reference to the calling FeatureFinder
			void setFeatureFinder(FeatureFinderType& ff)
			{
			  ff_ = &ff;
			}
		protected:
			///Pointer to the calling FeatureFinder that is used to access the data
			FeatureFinderType* ff_;
			
			/// @name Wrapper methods that allow the derived classes to access the data of the FeatureFinder friend class.
			//@{
			inline const FeatureMapType& getFeatureMap_() const 
			{ 
				return ff_->getFeatureMap_();
			}
			inline FeatureMapType& getFeatureMap_() 
			{ 
				return ff_->getFeatureMap_();
			}
			inline const MapType& getData_() const 
			{ 
				return ff_->getData_();
			}
	    inline const Flag& getPeakFlag_(const IDX& index) const
	    {
				return ff_->getPeakFlag_(index);
	    }
	    inline Flag& getPeakFlag_(const IDX& index) 
	    { 
				return ff_->getPeakFlag_(index);
	    }
	    inline IntensityType getPeakIntensity_(const IDX& index) const
	    { 
				return ff_->getPeakIntensity_(index);
	    }
	    inline CoordinateType getPeakMz_(const IDX& index) const
	    { 
				return ff_->getPeakMz_(index);
	    }
	    inline CoordinateType getPeakRt_(const IDX& index) const
	    { 
				return ff_->getPeakRt_(index);
			}			
	    inline void getNextMz_(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
				return ff_->getNextMz_(index);
	    }
	    inline void getPrevMz_(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
				return ff_->getPrevMz_(index);
	    }
  		void getNextRt_(IDX& index) throw (NoSuccessor, Exception::Precondition)
  		{
				return ff_->getNextRt_(index);
			}
			void getPrevRt_(IDX& index) throw (NoSuccessor, Exception::Precondition)
  		{
				return ff_->getPrevRt_(index);
			}
			void addConvexHull_(const IndexSet& set, Feature& f) const
			{
				return ff_->addConvexHull_(index, f);
			}
		//@}

		private:
			// not implemented -> private
			FeatureFinderAlgorithm& operator=(const FeatureFinderAlgorithm&);
			// not implemented -> private
			FeatureFinderAlgorithm(const FeatureFinderAlgorithm&);
	};
}

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimple.h>
namespace OpenMS
{
			template<class PeakType, class FeatureType>
	    void FeatureFinderAlgorithm<PeakType,FeatureType>::registerChildren()
			{
				Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct(FeatureFinderAlgorithmSimple<PeakType,FeatureType>::getProductName(), &FeatureFinderAlgorithmSimple<PeakType,FeatureType>::create);
			}
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
