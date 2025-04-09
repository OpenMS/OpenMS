// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSIsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSIsotopeDistribution* ptr = nullptr;
IMSIsotopeDistribution* null_ptr = nullptr;
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



