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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(Weights, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Weights* ptr = 0;
Weights* null_ptr = 0;
START_SECTION(Weights())
{
	ptr = new Weights();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Weights())
{
	delete ptr;
}
END_SECTION

START_SECTION((Weights(const alphabet_masses_type &masses, alphabet_mass_type prec)))
{
  // TODO
}
END_SECTION

START_SECTION((Weights(const Weights &other)))
{
  // TODO
}
END_SECTION

START_SECTION((Weights& operator=(const Weights &weights_)))
{
  // TODO
}
END_SECTION

START_SECTION((size_type size() const ))
{
  // TODO
}
END_SECTION

START_SECTION((weight_type getWeight(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecision(alphabet_mass_type precision_)))
{
  // TODO
}
END_SECTION

START_SECTION((alphabet_mass_type getPrecision() const ))
{
  // TODO
}
END_SECTION

START_SECTION((weight_type operator[](size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((weight_type back() const ))
{
  // TODO
}
END_SECTION

START_SECTION((alphabet_mass_type getAlphabetMass(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((alphabet_mass_type getParentMass(const std::vector< unsigned int > &decomposition) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void swap(size_type index1, size_type index2)))
{
  // TODO
}
END_SECTION

START_SECTION((bool divideByGCD()))
{
  // TODO
}
END_SECTION

START_SECTION((alphabet_mass_type getMinRoundingError() const ))
{
  // TODO
}
END_SECTION

START_SECTION((alphabet_mass_type getMaxRoundingError() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



