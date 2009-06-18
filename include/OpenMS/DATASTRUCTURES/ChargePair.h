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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CHARGEPAIR_H
#define OPENMS_DATASTRUCTURES_CHARGEPAIR_H

#include <iostream>
#include <vector>

#include <OpenMS/KERNEL/FeatureMap.h>


namespace OpenMS
{

  /**
    @brief Representation of a (putative) link between two Features, which stem from the same compound
	  but have different charge (including different adduct ions (H+, Na+, ..)
    
    A ChargePair represents an edge between two Features and specifies their respective charge and adducts,
    so that when decharged they can be explained as stemming from the same compound.
    
  	
  	@ingroup Datastructures
  */
  class OPENMS_DLLAPI ChargePair
  {
  	
	 public:
    ///@name Constructors and destructor
    //@{
		/// Default constructor
		ChargePair()
			: feature0_index_(0),
				feature1_index_(0),
				feature0_charge_(0),
				feature1_charge_(0),
				compomer_id_(),
				mass_diff_(0),
				is_active_(false)
		{
		}
		/// Constructor from map index, element index and Feature
		ChargePair(const Size& index0,
							 const Size& index1,
							 const Int& charge0,
							 const Int& charge1,
							 const Size& compomer_id, 
							 const DoubleReal& mass_diff, 
						   const bool active)
			:	feature0_index_(index0),
				feature1_index_(index1),
				feature0_charge_(charge0),
				feature1_charge_(charge1),
				compomer_id_(compomer_id),
				mass_diff_(mass_diff),
				is_active_(active)
		{
		}
		/// Copy constructor
		ChargePair(const ChargePair& rhs)
			: feature0_index_(rhs.feature0_index_),
				feature1_index_(rhs.feature1_index_),
				feature0_charge_(rhs.feature0_charge_),
				feature1_charge_(rhs.feature1_charge_),
				compomer_id_(rhs.compomer_id_),
				mass_diff_(rhs.mass_diff_),
				is_active_(rhs.is_active_)
		{
		}
		/// Assignment operator
		ChargePair& operator = (const ChargePair& rhs)
		{
				if (&rhs == this) return *this;

				feature0_index_ = rhs.feature0_index_;
				feature1_index_ = rhs.feature1_index_;
				feature0_charge_ = rhs.feature0_charge_;
				feature1_charge_ = rhs.feature1_charge_;
				compomer_id_ = rhs.compomer_id_;
				mass_diff_ = rhs.mass_diff_;
				is_active_ = rhs.is_active_;

				return *this;
		}
		
		/// Destructor
		virtual ~ChargePair()
		{
		}
    //@}
    
    //@name Accessors
    //@{
		/// Returns the charge (for element 0 or 1)
		Int getCharge(UInt pairID) const
		{
			if (pairID == 0) return feature0_charge_;
			else return feature1_charge_;
		}
		/// Set the charge (for element 0 or 1)
		void setCharge(UInt pairID, Int e)
		{
			if (pairID == 0) feature0_charge_ = e;
			else feature1_charge_ = e;
		}
				
		/// Returns the element index (for element 0 or 1)
		Size getElementIndex(UInt pairID) const
		{
			if (pairID == 0) return feature0_index_;
			else return feature1_index_;
		}
		/// Set the element index (for element 0 or 1)
		void setElementIndex(UInt pairID, Size e)
		{
			if (pairID == 0) feature0_index_ = e;
			else feature1_index_ = e;
		}

		
		/// Returns the Id of the compomer that explains the mass difference
		Size getCompomerId() const
		{
				return compomer_id_;
				// TODO throw exception				
		}
		/// Set the compomer id
		void setCompomerId(Size compomer_id)
		{
			compomer_id_ = compomer_id;
			// TODO throw exception		
		}
		
		/// Returns the mass difference
		DoubleReal getMassDiff() const
		{
			return mass_diff_;
		}
		/// Sets the mass difference
		void setMassDiff(DoubleReal mass_diff)
		{
			mass_diff_ = mass_diff;
		}

		/// is this pair realized?
		bool isActive() const
		{
				return is_active_;
		}
		
		void setActive(const bool active) 
		{
				is_active_ = active;
		}
		
		//@}
				
		/// Equality operator
		virtual bool operator == (const ChargePair& i) const
		{
			return ((feature0_index_ == i.feature0_index_) && 
							(feature1_index_ == i.feature1_index_) && 
							(feature0_charge_ == i.feature0_charge_) && 
							(feature1_charge_ == i.feature1_charge_) && 
							(compomer_id_ == i.compomer_id_) && 
							(mass_diff_ == i.mass_diff_) &&
							(is_active_ == i.is_active_) );
		}

		/// Equality operator
		virtual bool operator != (const ChargePair& i) const
		{
			return !(this->operator ==(i));
		}
			

	 protected:
    	
		/// Int of the first element within the FeatureMap
		Size feature0_index_;
		/// Int of the second element within the FeatureMap
		Size feature1_index_;
		/// Assumed charge of the first feature
		Int feature0_charge_;
		/// Assumed charge of the second feature
		Int feature1_charge_;
		/// id of compomer that explains the mass difference
		Size compomer_id_;
		/// mass difference (after explanation by compomer)
		DoubleReal mass_diff_;
		/// was this pair realized by ILP?
		bool is_active_;
  };

  ///Print the contents of a ChargePair to a stream.
  OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ChargePair& cons);
  	
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CHARGEPAIR_H
 
