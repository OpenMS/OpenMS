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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_ADDUCT_H
#define OPENMS_DATASTRUCTURES_ADDUCT_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS {
    
class OPENMS_DLLAPI Adduct
{
public:

		typedef std::vector< Adduct > AdductsType;

		/// Default C'tor
    Adduct();

		/// C'tor with initial charge
    Adduct(Int charge);

		/// C'tor for all members
    Adduct(Int charge, Int amount, DoubleReal singleMass, String formula, DoubleReal log_prob, DoubleReal rt_shift, const String label="");

		/// Increase amount of this adduct by factor @param m
    Adduct operator *(const Int m) const;
    /// Add two adducts amount if they are equal (defined by equal formula)
    Adduct operator +(const Adduct& rhs);
		/// Add other adducts amount to *this (equal formula required!)
    void operator +=(const Adduct& rhs);

		
		/// Print the contents of an Adduct to a stream.
    friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Adduct& a);

		/// Comparator
		friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b);
	
    //@{ Accessors
    const Int& getCharge() const;
    
    void setCharge(const Int& charge);

    const Int& getAmount() const;
    void setAmount(const Int& amount);
    
    const DoubleReal& getSingleMass() const;
    void setSingleMass(const DoubleReal& singleMass);

    const DoubleReal& getLogProb() const;
    void setLogProb(const DoubleReal& log_prob);
    
    const String& getFormula() const;
    void setFormula(const String& formula);
    
		const DoubleReal& getRTShift() const;
    const String& getLabel() const;
		//}
    
private:
    Int charge_; //< usually +1
    Int amount_; //< number of entities
    DoubleReal singleMass_; //< mass of a single entity
    DoubleReal log_prob_;   //< log probability of observing a single entity of this adduct
    String formula_;   //< chemical formula (parsable by EmpiricalFormula)
		DoubleReal rt_shift_; //< RT shift induced by a single entity of this adduct (this is for adducts attached prior to ESI, e.g. labeling)
		String label_; //< Label for this adduct (can be used to indicate heavy labels)

		String checkFormula_(const String & formula);

};

} // namespace OpenMS


#endif

