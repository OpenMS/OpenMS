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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSIsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSIsotopeDistribution* ptr = 0;
IMSIsotopeDistribution* null_ptr = 0;
START_SECTION(IMSIsotopeDistribution())
{
	ptr = new IMSIsotopeDistribution();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSIsotopeDistribution())
{
	delete ptr;
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(nominal_mass_type nominalMass=0)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(mass_type mass)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(const peaks_container &peaks, nominal_mass_type nominalMass=0)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((size_type size() const ))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator=(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IMSIsotopeDistribution &distribution) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IMSIsotopeDistribution &distribution) const ))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator*=(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator*=(unsigned int pow)))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getMass(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((abundance_type getAbundance(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getAverageMass() const ))
{
  // TODO
}
END_SECTION

START_SECTION((nominal_mass_type getNominalMass() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setNominalMass(nominal_mass_type nominalMass)))
{
  // TODO
}
END_SECTION

START_SECTION((masses_container getMasses() const ))
{
  // TODO
}
END_SECTION

START_SECTION((abundances_container getAbundances() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void normalize()))
{
  // TODO
}
END_SECTION

START_SECTION((bool empty() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ims::IMSIsotopeDistribution::Peak] Peak(mass_type mass=0.0, abundance_type abundance=0.0)))
{
  // TODO
}
END_SECTION

START_SECTION(([ims::IMSIsotopeDistribution::Peak] bool operator==(const Peak &peak) const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



