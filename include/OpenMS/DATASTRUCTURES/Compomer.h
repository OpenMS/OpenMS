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

#ifndef OPENMS_ANALYSIS_DATASTRUCTURES_COMPOMER_H
#define OPENMS_ANALYSIS_DATASTRUCTURES_COMPOMER_H

#include <vector>
#include <map>
#include <set>
#include <cstdlib>

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace std;

namespace OpenMS {

class OPENMS_DLLAPI Compomer
{
public:
	enum {LEFT=0, RIGHT};

	typedef std::map<String,Adduct> CompomerComponents;
	typedef std::set< std::pair<String, Int> > CompomerSide;
	
	/// default Constructor
	Compomer()
		: cmp_(),
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
		: cmp_(),
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

	/// Add a.amount of Adduct @param a to Compomer and update its properties
	void add(const Adduct& a)
	{
    if (cmp_.count(a.getFormula())==0)
    {
      cmp_[a.getFormula()] = a;
    }
    else
    {
      cmp_[a.getFormula()] += a;//update adduct´s amount
    }
		//
		net_charge_ += a.getAmount()*a.getCharge();
		mass_ += a.getAmount()*a.getSingleMass();
		pos_charges_ +=  std::max(a.getAmount()*a.getCharge(),0);
		neg_charges_ -=  std::min(a.getAmount()*a.getCharge(),0);
		log_p_ += std::abs((Real)a.getAmount())*a.getLogProb();
	}

	
	/**
	 *  indicates if these two compomers can coexist for one feature
	 * @param cmp The other Compomer we compare to
	 * @param left_this Indicates which "side"(negative or positive adducts) we are looking at. Negative adducts belong to the left side of the ChargePair.
	 * @param left_other See above.
	 */
	bool isConflicting(const Compomer& cmp, bool left_this, bool left_other, const Adduct& implicit_this, const Adduct& implicit_other) const
	{
		// build two sets (one for each compomer) and look at the difference of both. It should be empty!
		CompomerSide side_this = getCompomerSide_(left_this, implicit_this);
		CompomerSide side_other = cmp.getCompomerSide_(left_other, implicit_other);
		
		// size is equal - we need to check more thorough...
		if (side_this.size()==side_other.size())
		{
			side_this.insert(side_other.begin(),side_other.end());
			if (side_this.size()==side_other.size()) return false;
		}
    //
		std::cout << "found conflict!! between \n" << (*this) << "and\n" << cmp << " at sides i:" << (left_this?"left":"right") << " and j:" << (left_other?"left":"right") << "\n"
							<< "with implicits  i:" << implicit_this.getAmount() << " && j: " << implicit_other.getAmount() << "\n";
		return true; 
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

	/// get adducts with their abundance as compact string (amounts are absolute unless side=0)
	/// @param side Use -1 for left, +1 for right and 0 for both
	String getAdductsAsString(Int side=0) const
	{
		String r;
		std::map<String,Adduct>::const_iterator it=cmp_.begin();
		for (; it!=cmp_.end(); ++it)
		{
			Int f = it->second.getAmount();
			//if (it!=cmp_.begin()) r+= " ";
			EmpiricalFormula ef(it->first);
			ef = ef * std::abs(f);

			if ( (f < 0 && side==-1) ||
					 (f > 0 && side==+1))
			{
			  r += ef.getString();
			}
			else if (side==0)
			{
				r += String(it->second.getAmount()) + "(" + it->first + ")";
			}
		}

		return r;
	}
	
	/// Sort compomer by (in order of importance): net-charge, mass, probability
	friend OPENMS_DLLAPI bool operator< (const Compomer &c1, const Compomer &c2); 

	/// Print the contents of a Compomer to a stream.
	friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Compomer & cmp);

	/// Comparator
	friend OPENMS_DLLAPI bool operator==(const Compomer& a, const  Compomer& b);
	
private:	
	
	/// get the left(negative amount) or right(positive amount) side subset of the compomers adduct set
	/// @note The result will give absolute amounts
	CompomerSide getCompomerSide_(const bool left, const Adduct& implicit) const
	{
		CompomerSide side; //set()
		CompomerComponents::const_iterator it=cmp_.begin(); // vec<String, Adduct>
    Int extra_amount, normal_amount;

		for (; it!=cmp_.end(); ++it)
		{
			if (it->first == implicit.getFormula())
			{
				extra_amount = implicit.getAmount();
        OPENMS_POSTCONDITION(extra_amount >= 0, "Compomer.h::getCompomerSide_() has invalid implicit adduct!");
			}
			else extra_amount=0;
				
			normal_amount=0;
			if (left && (it->second.getAmount()<0))
			{
        normal_amount = -it->second.getAmount();
			}
			else if((!left) && (it->second.getAmount()>0))
			{
        normal_amount = it->second.getAmount();
			}
      else if(it->second.getAmount()==0)
      {
        OPENMS_POSTCONDITION(extra_amount > 0, "Compomer.h::getCompomerSide_() has invalid adduct amount of 0!");
      }

      if (normal_amount + extra_amount > 0)
      {
				side.insert(std::make_pair(it->first, normal_amount + extra_amount));
			}

		}

		return side;
	}
	
	CompomerComponents cmp_;
	Int net_charge_;
	DoubleReal mass_;
	Int pos_charges_;
	Int neg_charges_;
	DoubleReal log_p_;   // log probability of compomer
	Size id_;

}; // \Compomer

} // namespace OpenMS

#endif //OPENMS_ANALYSIS_DATASTRUCTURES_COMPOMER_H
