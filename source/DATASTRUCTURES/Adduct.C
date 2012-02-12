// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/Adduct.h>

using namespace std;

namespace OpenMS
{

  Adduct::Adduct()
  : charge_(0),
		amount_(0),
		singleMass_(0),
		log_prob_(0),
		formula_(),
		rt_shift_(0),
		label_()
  {
  }

  Adduct::Adduct(Int charge)
  : charge_(charge),
		amount_(0),
		singleMass_(0),
		log_prob_(0),
		formula_(),
		rt_shift_(0),
		label_()
  {
  }

  Adduct::Adduct(Int charge, Int amount, DoubleReal singleMass, String formula, DoubleReal log_prob, DoubleReal rt_shift, const String label)
  : charge_(charge),
		amount_(amount),
		singleMass_(singleMass),
		log_prob_(log_prob),
		rt_shift_(rt_shift),
		label_(label)
  {
		if (amount < 0) std::cerr << "Attention: Adduct received negative amount! (" << amount << ")\n";
		formula_ = checkFormula_(formula);
  }

  Adduct Adduct::operator *(const Int m) const
  {
      Adduct a = *this;
      a.amount_ *= m;
      return a;
  }
  Adduct Adduct::operator +(const Adduct& rhs)
  {
      if (this->formula_ != rhs.formula_)
      {
        throw "Adduct::Operator +()  tried to add incompatible adduct!";
      }
      Adduct a = *this;
      a.amount_ += rhs.amount_;
      return a;
  }

  void Adduct::operator +=(const Adduct& rhs)
  {
    if (this->formula_ != rhs.formula_)
    {
      throw "Adduct::Operator +=()  tried to add incompatible adduct!";
    }
    this->amount_ += rhs.amount_;
  }

  //@{ Accessors
  const Int& Adduct::getCharge() const
  {
    return charge_;
  }
  
  void Adduct::setCharge(const Int& charge)
  {
    charge_ = charge;
  }

  const Int& Adduct::getAmount() const
  {
    return amount_;
  }
  
  void Adduct::setAmount(const Int& amount)
  {
		if (amount < 0) std::cerr << "Warning: Adduct received negative amount! (" << amount << ")\n";
    amount_ = amount;
  }    

  const DoubleReal& Adduct::getSingleMass() const
  {
    return singleMass_;
  }
  
  void Adduct::setSingleMass(const DoubleReal& singleMass)
  {
    singleMass_ = singleMass;
  }       

  const DoubleReal& Adduct::getLogProb() const
  {
    return log_prob_;
  }
  
  void Adduct::setLogProb(const DoubleReal& log_prob)
  {
    log_prob_ = log_prob;
  }  

  const String& Adduct::getFormula() const
  {
    return formula_;
  }
  void Adduct::setFormula(const String& formula)
  {
		formula_ = checkFormula_(formula);
  }      
  
	const DoubleReal& Adduct::getRTShift() const
	{
		return rt_shift_;
	}
  const String& Adduct::getLabel() const
  {
    return label_;
  }

	String Adduct::checkFormula_(const String & formula)
	{
		EmpiricalFormula ef(formula);
		if (ef.getCharge()!=0)
		{
			std::cerr << "Warning: Adduct contains explicit charge (alternating mass)! (" << formula << ")\n";
		}
		if (ef.isEmpty())
		{
			std::cerr << "Warning: Adduct was given empty formula! (" << formula << ")\n";
		}
		if ( (ef.getNumberOfAtoms()>1) && (std::distance(ef.begin(),ef.end())==1) ) 
		{
			std::cerr << "Warning: Adduct was given only a single element but with an abundance>1. This might lead to errors! (" << formula << ")\n";
		}
		
		return ef.getString();			
	}

	///Print the contents of an Adduct to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Adduct& a)
	{
					os  << "---------- Adduct -----------------\n";
					os << "Charge: " << a.charge_ << std::endl;
					os << "Amount: " << a.amount_ << std::endl;
					os << "MassSingle: " << a.singleMass_ << std::endl;
					os << "Formula: " << a.formula_ << std::endl;
					os << "log P: " << a.log_prob_ << std::endl;
					return os;
	}
	
	OPENMS_DLLAPI bool operator==(const Adduct& a, const  Adduct& b)
	{
		return (a.charge_ == b.charge_ 
						&& a.amount_ == b.amount_
						&& a.singleMass_ == b.singleMass_
						&& a.log_prob_ == b.log_prob_
						&& a.formula_ == b.formula_);
				
	}

}

