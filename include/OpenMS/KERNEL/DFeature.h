// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $ 
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DFEATURE_H
#define OPENMS_KERNEL_DFEATURE_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>

namespace OpenMS
{

	/**	
		@brief A D-dimensional feature.
		
		A feature represents a subset of peaks in a map.  In general, it
		summarizes all peaks related to a specific peptide or chemical entity and
		thus reduces partitions of the LCMS dataset to a more meaningful entity.
		Peaks can be converted to features through the FeatureFinder.

		Features are usually contained in DFeatureMap.  Features themselves can
		either contain features again (composite design pattern) or their
		constituent peaks.

		Hierarchical relationships in features (i.e. features containing features
		containing features...) can be used to express group relationships.  For
		example, one might group the features corresponding to an ICAT pair into
		an aggregate ICAT feature.

		Features are derived from DPeak, as they inherit most of their properties.
		In particular, a feature has a position and an intensity (area).  The
		position of a feature is defined as maximum position of the model for the
		retention time dimension and the mass of the monoisotopic peak for the m/z
		dimension.
		
		@ingroup Kernel, Serialization
	*/
	template <Size D, typename Traits = KernelTraits>
	class DFeature 
		: public DPeak<D, Traits>
	{
		public:
		
		/** 
			@name Type definitions
		*/
		//@{	
		enum { DIMENSION = D };
		typedef Traits TraitsType;
		typedef DRange<D, TraitsType> BoundingBoxType;		
		typedef std::vector< DPosition<D, TraitsType>	> ConvexHullType;	
		typedef std::vector<ConvexHullType> ConvexHullVector;			
		typedef	typename Traits::RealType QualityType;
		/// 
		typedef	typename Traits::ChargeType ChargeType;
		//@}

		/** @name Constructors and Destructor
		*/
		//@{
		/// Default constructor 
		DFeature() 
			: DPeak<D, TraitsType>(),
				overall_quality_(),
				convex_hulls_(),
				charge_(0)
		{
			std::fill(qualities_,qualities_+D,0);
		}
		
		/// Copy constructor
		inline DFeature(const DFeature& feature) 
			: DPeak<D, TraitsType>(feature),
				overall_quality_(feature.overall_quality_),
				model_desc_(feature.model_desc_),
				convex_hulls_(feature.convex_hulls_),
				charge_(feature.charge_)
		{
			std::copy(feature.qualities_,feature.qualities_+D,qualities_);
		}

		/// Destructor
		~DFeature() {}
		//@}
		
		/**	@name Accessors	*/
		//@{
		/// Non-mutable access to the bounding box
		BoundingBoxType getBoundingBox() const;

		/// Non-mutable access to the overall quality
		inline const QualityType& getOverallQuality() const { return overall_quality_; }
		/// Mutable access to the overall quality
		inline QualityType& getOverallQuality() { return overall_quality_; }
		/// Set the overall quality
		inline void setOverallQuality(const QualityType& q) { overall_quality_ = q; }

		/// Non-mutable access to the quality in dimension c
		inline const QualityType& getQuality(const Position index) const
		{
			OPENMS_PRECONDITION(index < D, "DFeature<D>:getQuality(Position): index overflow!")
			return qualities_[index]; 
		}

		/// Mutable access to the quality in dimension c
		inline QualityType& getQuality(const Position index)
		{
			OPENMS_PRECONDITION(index < D, "DFeature<D>:getQuality(Position): index overflow!")
			return qualities_[index]; 
		}

		/// Set the quality in dimension c
		inline void setQuality(const Position index, const QualityType& q)
		{
			OPENMS_PRECONDITION(index < D, "DFeature<D>:setQuality(Position): index overflow!")
			qualities_[index] = q; 
		}
		
		/// Non-mutable access to charge state
		inline const ChargeType& getCharge() const { return charge_; }
		/// Mutable access to charge state
		inline ChargeType& getCharge() { return charge_; }
		/// Set charge state
		inline void setCharge(const ChargeType ch) { charge_ = ch; }

		/// Non-mutable access to the model description
		inline const ModelDescription<D,Traits>& getModelDescription() const { return model_desc_; }
		/// Mutable access to the model description
		inline ModelDescription<D,Traits>& getModelDescription() { return model_desc_; }
		/// Set the model description
		inline void setModelDescription(const ModelDescription<D,Traits>& q) { model_desc_ = q; }

		/// Non-mutable access to the convex hulls
		inline const ConvexHullVector& getConvexHulls() const { return convex_hulls_; }
		/// Mutable access to the convex hulls
		inline ConvexHullVector& getConvexHulls() { return convex_hulls_; }
		/// Set the convex hulls
		inline void setConvexHulls(const ConvexHullVector& hulls) { convex_hulls_ = hulls; }
		//@}

		/// Assignment operator
		DFeature& operator = (const DFeature& rhs);

		/// Equality operator
		bool operator == (const DFeature& rhs) const;
		
		/// Compare by getOverallQuality()
		struct OverallQualityLess
			: std::binary_function < DFeature, DFeature, bool >
		{
			bool operator () ( DFeature const & left, DFeature const & right ) const
			{
				return ( left.getOverallQuality() < right.getOverallQuality() );
			}
			bool operator () ( DFeature const & left, QualityType const & right ) const
			{
				return ( left.getOverallQuality() < right );
			}
			bool operator () ( QualityType const & left, DFeature const & right ) const
			{
				return ( left < right.getOverallQuality() );
			}
			bool operator () ( QualityType const & left, QualityType const & right ) const
			{
				return ( left < right );
			}
		};

		protected:
		/// Overall quality measure of the feature
		QualityType overall_quality_;
		/// Quality measures for each dimension
		QualityType qualities_[D];
		/// Description of the theoretical model the feature was constructed with
		ModelDescription<D,Traits> model_desc_;
		/// Array of convex hulls of the feature areas
		ConvexHullVector convex_hulls_;		
		/// Charge of the peptide represented by this feature
		ChargeType charge_;

		/**@name Serialization
		 */		
		//@{
	 public:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{
			ar & boost::serialization::make_nvp("dpeak",boost::serialization::base_object<DPeak<D,Traits> >(*this));
      ar & boost::serialization::make_nvp("overall_quality",overall_quality_);
      ar & boost::serialization::make_nvp("qualities",qualities_);
			// TODO: serialization of model_desc_
			// TODO: serialization of convex_hulls_
      ar & boost::serialization::make_nvp("charge",charge_);
		}
		//@}

		/// Serialization
		friend class boost::serialization::access;

	};


	template <Size D, typename Traits>
	DFeature<D, Traits>& DFeature<D, Traits>::operator = (const DFeature<D, Traits>& rhs)
	{
		if (this==&rhs) return *this;
		
		DPeak<D, TraitsType>::operator = (rhs);
		overall_quality_  = rhs.overall_quality_;
		std::copy(rhs.qualities_,rhs.qualities_+D,qualities_);
		model_desc_       = rhs.model_desc_;
		convex_hulls_     = rhs.convex_hulls_;
		charge_           = rhs.charge_;
		
		return *this;
	}

	template <Size D, typename Traits>
	bool DFeature<D, Traits>::operator == (const DFeature<D, Traits>& rhs) const
	{
		return (DPeak<D, TraitsType>::operator == (rhs) 
						&& (overall_quality_   == rhs.overall_quality_)
						&& (charge_ == rhs.charge_)
						&& std::equal(qualities_, qualities_+D, rhs.qualities_)
						&& (model_desc_ == rhs.model_desc_)
						&& (convex_hulls_ == rhs.convex_hulls_));
	}

	template <Size D, typename Traits>
	typename DFeature<D, Traits>::BoundingBoxType DFeature<D, Traits>::getBoundingBox() const
	{
		if (convex_hulls_.size()==0 || convex_hulls_[0].size()==0)
		{
			return DRange<D>();
		}
		else
		{
			DPosition<D> min(std::numeric_limits< typename DPeak<2,Traits>::CoordinateType>::max());
			DPosition<D> max(std::numeric_limits< typename DPeak<2,Traits>::CoordinateType>::min());
			
			// find min/max in each dimension of all convex hull points
			for (typename ConvexHullType::const_iterator	it=convex_hulls_[0].begin(); it!=convex_hulls_[0].end(); ++it)
			{
				for (UnsignedInt dim=0; dim<D; ++dim)
				{
					if ((*it)[dim]>max[dim]) max[dim] = (*it)[dim];
					if ((*it)[dim]<min[dim]) min[dim] = (*it)[dim];
				}
			}
			return DRange<D>(min,max);
		}
	}
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATURE_H
