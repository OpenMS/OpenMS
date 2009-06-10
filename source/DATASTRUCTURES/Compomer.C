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

#include <OpenMS/DATASTRUCTURES/Compomer.h>

using namespace std;

namespace OpenMS
{
	/// Sort compomer by (in order of importance): net-charge, mass, probability
	OPENMS_DLLAPI bool operator< (const Compomer &c1, const Compomer &c2)
	{  
		// how to sort Compomers:
		// first by net_charge
		if (c1.net_charge_ < c2.net_charge_) return true;
		else if (c1.net_charge_ > c2.net_charge_) return false;
		else
		{	
			// then my mass
			if (c1.mass_ < c2.mass_) return true;
			else if (c1.mass_ > c2.mass_) return false;
			else
			{
				// then by log probability (most probable compomers first!)
				return c1.log_p_ > c2.log_p_;
			}
		}
	} 
	
	/// Print the contents of a Compomer to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Compomer & cmp)
	{
		os << "Compomer: ";
		os << "Da " << cmp.mass_ << "; q_net " << cmp.net_charge_  << "; logP " << cmp.log_p_ << "[[ ";
		os << cmp.getAdductsAsString();
		os << " ]]\n";
		return os;
	}
	
	
	bool operator==(const Compomer& a, const  Compomer& b)
	{
		return (a.cmp_ == b.cmp_ 
						&& a.net_charge_ == b.net_charge_
						&& a.mass_ == b.mass_
						&& a.pos_charges_ == b.pos_charges_
						&& a.neg_charges_ == b.neg_charges_
						&& a.log_p_ == b.log_p_
						&& a.id_ == b.id_);
				
	}
}


