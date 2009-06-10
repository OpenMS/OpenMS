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

#include <OpenMS/DATASTRUCTURES/Adduct.h>

using namespace std;

namespace OpenMS
{
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

