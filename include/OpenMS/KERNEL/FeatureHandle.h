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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATUREHANDLE_H
#define OPENMS_KERNEL_FEATUREHANDLE_H

#include <iostream>
#include <vector>

#include <OpenMS/KERNEL/Feature.h>

namespace OpenMS
{
	class ConsensusFeature;
  /**
    @brief Representation of a Peak2D, RichPeak2D or Feature .
    
    The position and the intensity of the referenced feature are stored in the base class Peak2D.
    The original datapoint is referenced by the map and unique id.
  	
  	@ingroup Kernel
  */
  class OPENMS_DLLAPI FeatureHandle
  	: public Peak2D,
  	  public UniqueIdInterface
  {
  	
	 public:

		class FeatureHandleMutable_;

    ///@name Constructors and destructor
    //@{
		/// Default constructor
		FeatureHandle()
			: Peak2D(),
			  UniqueIdInterface(),
				map_index_(0),
				charge_(0)
		{
		}

		/// Constructor with map index, element index and position
		FeatureHandle(UInt64 map_index,const Peak2D& point, UInt64 element_index )
			: Peak2D(point),
			  map_index_(map_index),
				charge_(0)
		{
      setUniqueId(element_index);
		}

		/// Constructor from map index and basic feature
		FeatureHandle(UInt64 map_index, const BaseFeature& point)
			: Peak2D(point),
			  UniqueIdInterface(point),
				map_index_(map_index),
				charge_(point.getCharge())
		{
		}

		/// Copy constructor
		FeatureHandle(const FeatureHandle& rhs)
			: Peak2D(rhs),
        UniqueIdInterface(rhs),
				map_index_(rhs.map_index_),
				charge_(rhs.charge_)
		{
		}

		/// Assignment operator
		FeatureHandle& operator = (const FeatureHandle& rhs)
		{
			Peak2D::operator=(rhs);
			UniqueIdInterface::operator=(rhs);
			map_index_ = rhs.map_index_;
			charge_ = rhs.charge_;
			
			return *this;
		}
		/// Destructor
		virtual ~FeatureHandle()
		{
		}
		
		/**@brief Override (most of all) constness.

		We provide this such that you can modify instances FeatureHandle which are
		stored within a ConsensusFeature.  Note that std::set does not provide
		non-const iterators, because these could be used to change the relative
		ordering of the elements, and iterators are (by design/concept) unaware of
		their containers.  Since ConsensusFeature uses the ordering by IndexLess
		(which see), you <i>must not</i> modify the map index of element index if
		there is more than one FeatureHandle stored in a ConsensusFeature.
		Consequently, we "disable" setMapIndex() or setElementIndex() even within
		FeatureHandleMutable_.  On the other hand, it is perfectly safe to apply
		FeatureHandle::setRT(), FeatureHandle::setMZ(),
		FeatureHandle::setIntensity(), FeatureHandle::setCharge(), etc..
		*/
		FeatureHandleMutable_& asMutable() const;
    //@}
    
    ///@name Accessors
    //@{
		/// Returns the map index
		UInt64 getMapIndex() const
		{
			return map_index_;
		}
		/// Set the map index
		void setMapIndex(UInt64 i)
		{
			map_index_ = i;
		}

		/// Sets the charge
		void setCharge(Int charge)
		{
			charge_ = charge;
		}
		/// Returns the charge
		Int getCharge() const
		{
			return charge_;
		}
		//@}
				
		/// Equality operator
		bool operator == (const FeatureHandle& i) const
		{
			return  (Peak2D::operator==(i)) &&
			    (UniqueIdInterface::operator==(i)) &&
			    (map_index_ == i.map_index_) &&
			    (charge_ == i.charge_);
		}

		/// Equality operator
		bool operator != (const FeatureHandle& i) const
		{
			return !(operator==(i));
		}
			
		///Comparator by map and unique id 
		struct IndexLess
			: std::binary_function < FeatureHandle, FeatureHandle, bool >
		{
			bool operator () ( FeatureHandle const & left, FeatureHandle const & right ) const
			{
				//if map indices are equal, use unique ids
				if ( left.map_index_ == right.map_index_)
				{
					return left.getUniqueId() < right.getUniqueId();
				}
				//else use map indices
				return ( left.map_index_ < right.map_index_ );
			}
		};

	 protected:
    	
		/// Index of the element's container
    UInt64 map_index_;
		/// Charge of the feature
		Int charge_;
  };

	/**@brief Helper class returned by FeatureHandle::asMutable(), which see.

	Note that the mutators for unique id and map index are declared private.
	This is done because these are used by IndexLess comparator.  This way it is
	a bit harder to use FeatureHandle::asMutable() for illegal purposes ;-)
	 */
	class OPENMS_DLLAPI FeatureHandle::FeatureHandleMutable_ : public FeatureHandle
	{
	 private:
		FeatureHandle::setUniqueId;
		FeatureHandle::setMapIndex;
		FeatureHandleMutable_();
		FeatureHandleMutable_(const FeatureHandleMutable_&);
	};

	inline FeatureHandle::FeatureHandleMutable_& FeatureHandle::asMutable() const
	{
		// the const cast is to remove constness, but note that FeatureHandleMutable_ lacks some mutators
		return static_cast<FeatureHandleMutable_&>(const_cast<FeatureHandle&>(*this));
	}
	
  ///Print the contents of a FeatureHandle to a stream.
  OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const FeatureHandle& cons);
  	
} // namespace OpenMS

#endif // OPENMS_KERNEL_FEATUREHANDLE_H
