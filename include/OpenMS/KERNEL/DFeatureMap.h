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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DFEATUREMAP_H
#define OPENMS_KERNEL_DFEATUREMAP_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/DATASTRUCTURES/RangeManager.h>

#include <algorithm>
#include <vector>

namespace OpenMS
{

	/**	
		@brief A container for (composite) features.
		
		A map is a container holding D-dimensional features,
		which in turn represent chemical entities (peptides, proteins, etc.) found
		in a D-dimensional experiment.
		Maps are implemented as vectors of features and have basically the same interface
		as an STL vector has (model of Random Access Container and Back Insertion Sequence).
		Maps are typically created from peak data of 2D runs through the FeatureFinder.
				
		@todo post-release: derive Dimension D from FeatureT and thus eliminate 1st template argument (Marcel)

		@ingroup Kernel, Serialization
	*/
	template <Size D, typename FeatureT = DFeature<D> >
	class DFeatureMap
		: public std::vector<FeatureT>,
			public RangeManager<D, typename FeatureT::TraitsType>,
			public ExperimentalSettings
	{
	 public:
			/**	
				 @name Type definitions
			*/
			//@{
			typedef typename FeatureT::TraitsType TraitsType;
			typedef FeatureT FeatureType;
			typedef RangeManager<D, TraitsType> RangeManagerType;
			typedef std::vector<FeatureType> Base;
			typedef typename Base::iterator Iterator;
			typedef typename Base::const_iterator ConstIterator;
			typedef typename Base::reverse_iterator ReverseIterator;
			typedef typename Base::const_reverse_iterator ConstReverseIterator;
			typedef FeatureType& Reference;
			typedef const FeatureType& ConstReference;
	
			//@}
		
    enum { DIMENSION = FeatureType::DIMENSION };
		
			/**	
				 @name Constructors and Destructor
			*/
			//@{
			
			/// Default constructor
			DFeatureMap()
				: Base(),
					RangeManagerType(),
					ExperimentalSettings()
			{
				
			}
			
			/// Copy constructor
			DFeatureMap(const DFeatureMap& map) 
				: Base(map),
					RangeManagerType(map),
					ExperimentalSettings(map)
			{
			
			}
			
			/// Destructor
			virtual ~DFeatureMap()
			{
				
			}
			
			//@}
				
			/// Assignment operator
			DFeatureMap& operator = (const DFeatureMap& rhs)
			{
				if (&rhs==this) return *this;
					
				Base::operator=(rhs);
				RangeManagerType::operator=(rhs);
				ExperimentalSettings::operator=(rhs);
				
				return *this;
			}
	
			/// Equality operator
			bool operator == (const DFeatureMap& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					RangeManagerType::operator==(rhs) &&
					ExperimentalSettings::operator==(rhs) 
					;				
			}
				
			/// Equality operator
			bool operator != (const DFeatureMap& rhs) const
			{
				return !(operator==(rhs));
			}
				
			/** @brief Sort features by intensity. */
			void sortByIntensity() 
			{ 
				typename DFeatureMap<D>::iterator beg = this->begin();
				typename DFeatureMap<D>::iterator ed  = this->end();
				std::sort(beg, ed, typename FeatureType::IntensityLess() ); 
			}
				
			/** @brief Sort features by position.
				
				Lexicographical sorting from dimention 0 to dimension D is performed.			
			*/
			void sortByPosition() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::PositionLess() );
			}
			
			/** @brief Sort features by position @p i.
				
			   Features are only sorted by position @p i.		
			*/
			void sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented);
				
			void sortByOverallQuality()
			{
				typename DFeatureMap<D>::iterator beg = this->begin();
				typename DFeatureMap<D>::iterator ed  = this->end();
				std::sort(beg, ed, typename FeatureType::OverallQualityLess() ); 
			}
	
			/// Serialization interface
			template<class Archive>
			void serialize(Archive & ar, const unsigned int /* version */ )
			{
				ar & boost::serialization::make_nvp("vector",boost::serialization::base_object<std::vector<FeatureT> >(*this));
				// TODO: serialization of base object ExperimentalSettings
			}
			/// Serialization
			friend class boost::serialization::access;
			
			// Docu in base class
			void updateRanges()
			{
				this->clearRanges();
				updateRanges_(this->begin(),this->end());
			}
	};
	
	/// Print content of a feature map to a stream.
	template <Size D, typename FeatureType >
	std::ostream& operator << (std::ostream& os, const DFeatureMap<D, FeatureType>& map)
	{
		os << "# -- DFEATUREMAP BEGIN --"<< std::endl;
		os << "# POSITION \tINTENSITY\tOVERALLQUALITY\tCHARGE" << std::endl; 
		for (typename DFeatureMap<D>::const_iterator iter = map.begin(); iter!=map.end(); iter++)
		{
			os << iter->getPosition() << '\t'
				 << iter->getIntensity() << '\t'
				 << iter->getOverallQuality() << '\t'
				 << iter->getCharge()
				 << std::endl;
		}
		os << "# -- DFEATUREMAP END --"<< std::endl;
		return os;
	}
	
	template <Size D, typename FeatureType > 
	void DFeatureMap<D,FeatureType>::sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented)
	{ 
		OPENMS_PRECONDITION(i < Index(D), "illegal dimension")
		if (i==0)
		{
			std::sort(Base::begin(), Base::end(), typename FeatureType::template NthPositionLess<0>() );
		}
		else if (i==1)
		{
			std::sort(Base::begin(), Base::end(), typename FeatureType::template NthPositionLess<1>() );
		}
		else if (i==2)
		{
			std::sort(Base::begin(), Base::end(), typename FeatureType::template NthPositionLess<2>() );
		}
		else
		{
			throw Exception::NotImplemented(__FILE__,__LINE__,__FUNCTION__);
		}
	}
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
