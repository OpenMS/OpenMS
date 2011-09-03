// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Element.h>

namespace OpenMS {

namespace ims {

/**
 * @note Value for electron mass is taken from 
 * @link www.mcelwee.net/html/table_of_physical_constants.html
 */
const Element::mass_type Element::ELECTRON_MASS_IN_U = 0.00054858;

Element& Element::operator =(const Element& element)
{
	// if one doesn't assign object to itself,
	// assign all object elements to the elements of the given object
  if (this != &element)
  {
		name = element.name;		
		sequence = element.sequence;
		isotopes = element.isotopes;
	}
	return *this;
}


bool Element::operator ==(const Element& element) const
{
	return ( this == &element ||
          (name == element.name &&
           sequence == element.sequence &&
           isotopes == element.isotopes));
}


bool Element::operator !=(const Element& element) const
{
	return !this->operator==(element);
}


std::ostream& operator <<(std::ostream& os, const Element& element)
{
	os << "name:\t" << element.getName() << "\nsequence:\t" << element.getSequence()
	   << "\nisotope distribution:\n" << element.getIsotopeDistribution() << '\n';
	return os;
}

} // namespace ims
} // namespace OpenMS
