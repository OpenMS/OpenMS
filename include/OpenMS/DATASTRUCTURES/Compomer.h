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

#ifndef OPENMS_DATASTRUCTURES_COMPOMER_H
#define OPENMS_DATASTRUCTURES_COMPOMER_H

#include <vector>
#include <map>
#include <set>
#include <cstdlib>

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{

/**
	@brief TODO
	
	@todo Add docu here (Chris)
*/
class OPENMS_DLLAPI Compomer
{
public:
	/// side of compomer (LEFT ^ substract; RIGHT ^ add)
	enum SIDE {LEFT, RIGHT, BOTH};
	
	typedef Map<String,Adduct> CompomerSide; /// adducts and their abundance etc
	typedef std::vector< CompomerSide > CompomerComponents; /// container for the two sides [0]=left, [1]=right
	
	/// Default Constructor
	Compomer()
		: cmp_(2),
		net_charge_(0),
		mass_(0),
		pos_charges_(0),
		neg_charges_(0),
		log_p_(0),
		id_(0)
	{
	}

	/// Constructor with net-charge and mass
	Compomer(Int net_charge, DoubleReal mass, DoubleReal log_p)
		: cmp_(2),
		net_charge_(net_charge),
		mass_(mass),
		pos_charges_(0),
		neg_charges_(0),
		log_p_(log_p),
		id_(0)
	{
	}
	
	/// Copy C'tor
	Compomer(const Compomer& p)
		: cmp_(p.cmp_),
			net_charge_(p.net_charge_),
			mass_(p.mass_),
			pos_charges_(p.pos_charges_),
			neg_charges_(p.neg_charges_),
			log_p_(p.log_p_),
			id_(p.id_)
	{
	}
	
	/// Assignment Operator
	Compomer& operator=(const Compomer& source)
	{
		if (&source == this) return *this;
		cmp_ = source.cmp_;
		net_charge_ = source.net_charge_;
		mass_ = source.mass_;
		pos_charges_ = source.pos_charges_;
		neg_charges_ = source.neg_charges_;
		log_p_ = source.log_p_;
		id_ = source.id_;
		
		return *this;
	}	

	/// Add a.amount of Adduct @param a to Compomer's @param side and update its properties
	void add(const Adduct& a, UInt side)
	{
		if (side >= BOTH) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Compomer::add() does not support this value for 'side'!", String(side));

		if (a.getAmount() <=0) 
		{	std::cerr << "Compomer::add() was given adduct with negative amount! Are you sure this is what you want?!\n"; }
		if (a.getCharge() <=0) 
		{	std::cerr << "Compomer::add() was given adduct with negative charge! Are you sure this is what you want?!\n"; }

    if (cmp_[side].count(a.getFormula())==0)
    {
      cmp_[side][a.getFormula()] = a;
    }
    else
    {
      cmp_[side][a.getFormula()] += a;//update adduct´s amount
    }
		int mult[] = {-1,1};
		net_charge_ += a.getAmount()*a.getCharge()*mult[side];
		mass_ += a.getAmount()*a.getSingleMass()*mult[side];
		pos_charges_ +=  std::max(a.getAmount()*a.getCharge()*mult[side],0);
		neg_charges_ -=  std::min(a.getAmount()*a.getCharge()*mult[side],0);
		log_p_ += std::abs((Real)a.getAmount())*a.getLogProb();
	}

	
	/**
	 *  indicates if these two compomers can coexist for one feature
	 * @param cmp The other Compomer we compare to
	 * @param side_this Indicates which "side"(negative or positive adducts) we are looking at. Negative adducts belong to the left side of the ChargePair.
	 * @param side_other See above.
	 */
	bool isConflicting(const Compomer& cmp, UInt side_this, UInt side_other) const
	{
		if (side_this  >= BOTH) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Compomer::isConflicting() does not support this value for 'side_this'!", String(side_this));
		if (side_other >= BOTH) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Compomer::isConflicting() does not support this value for 'side_other'!", String(side_other));
		
		bool conflict_found = false;
		
		// size is equal - we need to check more thorough...
		if (cmp_[side_this].size()==cmp.getComponent()[side_other].size())
		{
			for (CompomerSide::const_iterator it = cmp_[side_this].begin(); it!=cmp_[side_this].end(); ++it)
			{
				if ((!cmp.getComponent()[side_other].has(it->first)) || cmp.getComponent()[side_other][it->first].getAmount() != it->second.getAmount())
				{
					conflict_found = true;
					break;
				}
			}
		}
		else conflict_found = true;
    //
		//if (conflict_found) std::cout << "found conflict!! between \n" << (*this) << "and\n" << cmp << " at sides i:" << (left_this?"left":"right") << " and j:" << (left_other?"left":"right") << "\n"
		//										<< "with implicits  i:" << implicit_this.getAmount() << " && j: " << implicit_other.getAmount() << "\n";
		return conflict_found; 
	}

	/// set an Id which allows unique identification of a compomer
	void setID(const Size& id)
	{
		id_ = id;
	}
	/// return Id which allows unique identification of this compomer
	const Size& getID() const
	{
		return id_;
	}


	/// return Id which allows unique identification of this compomer
	const CompomerComponents& getComponent() const
	{
		return cmp_;
	}

	/// net charge of compomer (i.e. difference between left and right side of compomer)
	const Int& getNetCharge() const
	{
		return net_charge_;
	}

	/// mass of all contained adducts
	const DoubleReal& getMass() const
	{
		return mass_;
	}

	/// summed positive charges of contained adducts
	const Int& getPositiveCharges() const
	{
		return pos_charges_;
	}
	
	/// summed negative charges of contained adducts
	const Int& getNegativeCharges() const
	{
		return neg_charges_;
	}
	
	/// return log probability
	const DoubleReal& getLogP() const
	{
		return log_p_;
	}	

	/// get adducts with their abundance as compact string for both sides
	String getAdductsAsString() const
	{
		return "(" + getAdductsAsString(LEFT) + ")-(" + getAdductsAsString(RIGHT) + ")";
	}

	/// get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
	/// @param side Use LEFT for left, RIGHT for right
	String getAdductsAsString(UInt side) const
	{
		if (side >= BOTH) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Compomer::getAdductsAsString() does not support this value for 'side'!", String(side));
		
		String r;
		CompomerSide::const_iterator it=cmp_[side].begin();
		for (; it!=cmp_[side].end(); ++it)
		{
			Int f = it->second.getAmount();
			
			if (it->first.has('+')) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "An Adduct contains ímplicit charge. This is not allowed!", it->first);
			
			EmpiricalFormula ef(it->first);
			ef = ef * f;
		  r += ef.getString();
		}

		return r;
	}
	
	/*
	/// check if Compomer only contains a single adduct
	bool isSimpleAdduct(Adduct& a)
	{
		if ((cmp_[LEFT].size() > 1) || (cmp_[RIGHT].size() > 1) return false;
		
		if ( ((cmp_[LEFT].size()==1) && (cmp_[LEFT].count(a.getFormula())==0)) ||
				 ((cmp_[RIGHT].size()==1) && (cmp_[RIGHT].count(a.getFormula())==0)) ) return false;
				
		return true;
	}
	*/
	
	/**
		@brief Remove all adducts of type @p a
	
		Remove ALL instances of the given adduct,
		BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
	**/
	Compomer removeAdduct(const Adduct& a) const
	{
		Compomer tmp = removeAdduct(a, LEFT);
		tmp = tmp.removeAdduct(a, RIGHT);
		return tmp;
	}
	
	/**
		@brief Remove all adducts of type @p a from @p side (LEFT or RIGHT)
		
		remove ALL instances of the given adduct from the given side (LEFT or RIGHT),
		BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
	*/
	Compomer removeAdduct(const Adduct& a, const UInt side) const
	{
		if (side >= BOTH) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Compomer::removeAdduct() does not support this value for 'side'!", String(side));
	
		Compomer tmp(*this);
		if (tmp.cmp_[side].count(a.getFormula())>0)
		{
			{ // how many instances does this side contain?
				Int amount = tmp.cmp_[side][a.getFormula()].getAmount();
				int mult[] = {-1,1};
				//const Adduct &to_remove = tmp.cmp_[side][a.getFormula()];
				tmp.net_charge_ -= amount *a.getCharge()*mult[side];
				tmp.mass_ -= amount *a.getSingleMass()*mult[side];
				tmp.pos_charges_ -=  std::max(amount *a.getCharge()*mult[side],0);
				tmp.neg_charges_ -= -std::min(amount *a.getCharge()*mult[side],0);
				tmp.log_p_ -= std::abs((Real)amount) *a.getLogProb();
			}
			// remove entry from map
			tmp.cmp_[side].erase(a.getFormula());
		}
		
		return tmp;
	}
	
	/// Adds @p add_side to this compomer.
	void add(const CompomerSide& add_side, UInt side)
	{
		for (CompomerSide::const_iterator it=add_side.begin(); it!=add_side.end(); ++it)
		{
			this->add(it->second, side);
		}
	}
	
/*	Size augmentSide_(const CompomerSide& to_augment, const Adduct& to_replace, UInt side)
	{
		Size rvalue = 0;
		Size free_charges;
		Size free_charges_used = 0;
		if (cmp_[side].count(to_replace.getFormula())==0)  free_charges = 0;
    else free_charges = cmp_[side][a.getFormula()].getAmount() * cmp_[side][a.getFormula()].getCharge();
    
    for (CompomerSide::iterator it=to_augment.begin(); it!=to_augment.end(); ++it)
		{
			Size required_amount = it->second.getAmount();
			if (cmp_[side].count(it->first)>0) required_amount = std::max(0, required_amount - cmp_[side][it->first]);
			
			if (required_amount>0)
			{
				*this = removeAdduct(cmp_[side][it->first], side);
				add( it->second, side);
				rvalue = 1;
			}
			
			free_charges_used += required_amount * it->second.getCharge();
		}
		if (free_charges_used > free_charges) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "More replacements requested than available from replaceable adduct while augmenting compomer " + String(free_charges_used) + " > " + String(free_charges) << "!", free_charges_used);

		// remove replacement adducts
		if (free_charges_used>0)
		{
			if ((free_charges - free_charges_used) % to_replace.getCharge() != 0) return -1; // variable adduct cannot exactly compensate the augmented adducts
	
			Size amount = (free_charges - free_charges_used) / to_replace.getCharge();

			*this = removeAdduct(to_replace, side);
			Adduct refill(to_replace);
			refill.setAmount(amount);
			add( refill, side);
			rvalue = 1;
		}
		
		return rvalue;
	}*/
	
	
	
	/// Sort compomer by (in order of importance): net-charge, mass, probability
	friend OPENMS_DLLAPI bool operator< (const Compomer &c1, const Compomer &c2); 

	/// Print the contents of a Compomer to a stream.
	friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Compomer & cmp);

	/// Comparator
	friend OPENMS_DLLAPI bool operator==(const Compomer& a, const  Compomer& b);
	
private:	
	
	CompomerComponents cmp_;
	Int net_charge_;
	DoubleReal mass_;
	Int pos_charges_;
	Int neg_charges_;
	DoubleReal log_p_;   // log probability of compomer
	Size id_;

}; // \Compomer

} // namespace OpenMS

#endif //OPENMS_DATASTRUCTURES_COMPOMER_H
