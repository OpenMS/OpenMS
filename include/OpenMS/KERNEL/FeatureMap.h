// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_FEATUREMAP_H
#define OPENMS_KERNEL_FEATUREMAP_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <algorithm>
#include <vector>

namespace OpenMS
{

	/**	
		@brief A container for (composite) features.
		
		A map is a container holding 2-dimensional features,
		which in turn represent chemical entities (peptides, proteins, etc.) found
		in a 2-dimensional experiment.
		Maps are implemented as vectors of features and have basically the same interface
		as an STL vector has (model of Random Access Container and Back Insertion Sequence).
		Maps are typically created from peak data of 2D runs through the FeatureFinder.
		
		@ingroup Kernel
	*/
	template <typename FeatureT = Feature >
	class FeatureMap
		: public std::vector<FeatureT>,
			public RangeManager<2>,
			public ExperimentalSettings
	{
	 public:
			/**	
				 @name Type definitions
			*/
			//@{
			typedef FeatureT FeatureType;
			typedef RangeManager<2> RangeManagerType;
			typedef std::vector<FeatureType> Base;
			typedef typename Base::iterator Iterator;
			typedef typename Base::const_iterator ConstIterator;
			typedef typename Base::reverse_iterator ReverseIterator;
			typedef typename Base::const_reverse_iterator ConstReverseIterator;
			typedef FeatureType& Reference;
			typedef const FeatureType& ConstReference;
	
			//@}
			/**	
				 @name Constructors and Destructor
			*/
			//@{
			
			/// Default constructor
			FeatureMap()
				: Base(),
					RangeManagerType(),
					ExperimentalSettings()
			{
				
			}
			
			/// Copy constructor
			FeatureMap(const FeatureMap& map) 
				: Base(map),
					RangeManagerType(map),
					ExperimentalSettings(map)
			{
			
			}
			
			/// Destructor
			virtual ~FeatureMap()
			{
				
			}
			
			//@}
				
			/// Assignment operator
			FeatureMap& operator = (const FeatureMap& rhs)
			{
				if (&rhs==this) return *this;
					
				Base::operator=(rhs);
				RangeManagerType::operator=(rhs);
				ExperimentalSettings::operator=(rhs);
				
				return *this;
			}
	
			/// Equality operator
			bool operator == (const FeatureMap& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					RangeManagerType::operator==(rhs) &&
					ExperimentalSettings::operator==(rhs) 
					;				
			}
				
			/// Equality operator
			bool operator != (const FeatureMap& rhs) const
			{
				return !(operator==(rhs));
			}
				
			/// Sorts the peaks according to ascending intensity
			void sortByIntensity(bool reverse=false) 
			{ 
				if (reverse)
				{
					std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::IntensityLess()) );
				}
				else
				{
					std::sort(this->begin(), this->end(), typename FeatureType::IntensityLess() ); 
				}
			}
				
			///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
			void sortByPosition() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::PositionLess() );
			}
			
			///Sort features by RT position.
			void sortByRT() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::RTLess() );
			}

			///Sort features by m/z position.
			void sortByMZ() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::MZLess() );
			}
			
			///Sort features by overall quality.
			void sortByOverallQuality()
			{
				std::sort(this->begin(), this->end(), typename FeatureType::OverallQualityLess() ); 
			}
			
			// Docu in base class
			void updateRanges()
			{
				this->clearRanges();
				updateRanges_(this->begin(),this->end());
				
				//enlarge the range by the convex hull points
				for (UInt i=0; i<this->size(); ++i)
				{
					DBoundingBox<2> box = this->operator[](i).getConvexHull().getBoundingBox();
					if (!box.isEmpty())
					{
						//update RT
						if (box.min()[Peak2D::RT] < this->pos_range_.min()[Peak2D::RT])
						{
							this->pos_range_.setMinX(box.min()[Peak2D::RT]);
						}
						if (box.max()[Peak2D::RT] > this->pos_range_.max()[Peak2D::RT])
						{
							this->pos_range_.setMaxX(box.max()[Peak2D::RT]);
						}
						//update m/z
						if (box.min()[Peak2D::MZ] < this->pos_range_.min()[Peak2D::MZ])
						{
							this->pos_range_.setMinY(box.min()[Peak2D::MZ]);
						}
						if (box.max()[Peak2D::MZ] > this->pos_range_.max()[Peak2D::MZ])
						{
							this->pos_range_.setMaxY(box.max()[Peak2D::MZ]);
						}
					}
				}
			}

			/// Swaps the content of this map with the content of @p from
			void swap(FeatureMap& from)
			{
				FeatureMap tmp;
				
				//swap range information
				tmp.RangeManagerType::operator=(*this);
				this->RangeManagerType::operator=(from);
				from.RangeManagerType::operator=(tmp);
				
				//swap experimental settings
				tmp.ExperimentalSettings::operator=(*this);
				this->ExperimentalSettings::operator=(from);
				from.ExperimentalSettings::operator=(tmp);
				
				//swap features
				Base::swap(from);
			}

	};
	
	/// Print content of a feature map to a stream.
	template <typename FeatureType >
	std::ostream& operator << (std::ostream& os, const FeatureMap<FeatureType>& map)
	{
		os << "# -- DFEATUREMAP BEGIN --"<< std::endl;
		os << "# POSITION \tINTENSITY\tOVERALLQUALITY\tCHARGE" << std::endl; 
		for (typename FeatureMap<FeatureType>::const_iterator iter = map.begin(); iter!=map.end(); iter++)
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
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
