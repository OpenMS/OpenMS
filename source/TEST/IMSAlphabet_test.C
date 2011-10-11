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
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSAlphabet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSAlphabet* ptr = 0;
IMSAlphabet* null_ptr = 0;
START_SECTION(IMSAlphabet())
{
	ptr = new IMSAlphabet();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSAlphabet())
{
	delete ptr;
}
END_SECTION

START_SECTION((const element_type& getElement(const name_type &name) const ))
{
  // TODO
}
END_SECTION

START_SECTION((const name_type& getName(size_type index) const ))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getMass(const name_type &name) const ))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getMass(size_type index) const ))
{
  // TODO
}
END_SECTION

START_SECTION((masses_type getMasses(size_type isotope_index=0) const ))
{
  // TODO
}
END_SECTION

START_SECTION((masses_type getAverageMasses() const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasName(const name_type &name) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void push_back(const name_type &name, mass_type value)))
{
  // TODO
}
END_SECTION

START_SECTION((void push_back(const element_type &element)))
{
  // TODO
}
END_SECTION

START_SECTION((void clear()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void sortByNames()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void sortByValues()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void load(const std::string &fname)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void load(const std::string &fname, IMSAlphabetParser<> *parser)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~IMSAlphabet()))
{
  // TODO
}
END_SECTION

START_SECTION((IMSAlphabet(const container &elements)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSAlphabet(const IMSAlphabet &alphabet)))
{
  // TODO
}
END_SECTION

START_SECTION((size_type size() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const element_type& getElement(size_type index) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setElement(const name_type &name, mass_type mass, bool forced=false)))
{
  // TODO
}
END_SECTION

START_SECTION((bool erase(const name_type &name)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



