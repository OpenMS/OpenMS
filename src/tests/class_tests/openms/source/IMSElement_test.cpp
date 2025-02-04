// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/ResidueDB.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSElement, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSElement* ptr = nullptr;
IMSElement* null_ptr = nullptr;

IMSIsotopeDistribution::peaks_container peaks;
peaks.push_back(IMSIsotopeDistribution::peak_type(0.0078250319,.999885));
peaks.push_back(IMSIsotopeDistribution::peak_type(0.01410178,.000115));
peaks.push_back(IMSIsotopeDistribution::peak_type(0.01604927,.0));

IMSIsotopeDistribution iso(peaks, 1);

START_SECTION(IMSElement())
{
	ptr = new IMSElement();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSElement())
{
	delete ptr;
}
END_SECTION


IMSElement * hydrogen = nullptr;

START_SECTION((IMSElement(const name_type &name, const isotopes_type &isotopes)))
{
  hydrogen = new IMSElement("H",iso);
  TEST_NOT_EQUAL(hydrogen, null_ptr)
  TEST_STRING_EQUAL(hydrogen->getName(), "H")
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), iso)
}
END_SECTION

START_SECTION((IMSElement(const IMSElement &element)))
{
  IMSElement hydrogen_copy(*hydrogen);
  TEST_EQUAL(hydrogen->getAverageMass(), hydrogen_copy.getAverageMass())
  TEST_EQUAL(hydrogen->getIonMass(), hydrogen_copy.getIonMass())
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), hydrogen_copy.getIsotopeDistribution())
  TEST_EQUAL(hydrogen->getMass(), hydrogen_copy.getMass())
  TEST_EQUAL(hydrogen->getName(), hydrogen_copy.getName())
  TEST_EQUAL(hydrogen->getNominalMass(), hydrogen_copy.getNominalMass())
  TEST_EQUAL(hydrogen->getSequence(), hydrogen_copy.getSequence())
}
END_SECTION

START_SECTION((IMSElement(const name_type &name, mass_type mass)))
{
  double oxygen_mass = 15.9994;
  IMSElement element("O", oxygen_mass);
  IMSIsotopeDistribution oxygen(oxygen_mass);

  TEST_EQUAL(element.getName(), "O")
  TEST_EQUAL(element.getNominalMass(), 0)
  TEST_EQUAL(element.getMass(), oxygen_mass)
  TEST_EQUAL(element.getIsotopeDistribution(), oxygen)
}
END_SECTION

START_SECTION((IMSElement(const name_type &name, nominal_mass_type nominal_mass=0)))
{
  IMSElement::nominal_mass_type nominal_mass = 16;
  IMSElement element("O", nominal_mass);
  IMSIsotopeDistribution oxygen(nominal_mass);

  TEST_EQUAL(element.getName(), "O")
  TEST_EQUAL(element.getNominalMass(), nominal_mass)
  TEST_EQUAL(element.getIsotopeDistribution(), oxygen)
}
END_SECTION

START_SECTION((const name_type& getName() const ))
{
  TEST_STRING_SIMILAR(hydrogen->getName(), "H")
}
END_SECTION

START_SECTION((void setName(const name_type &name)))
{
  hydrogen->setName("D");
  TEST_STRING_SIMILAR(hydrogen->getName(), "D")
  hydrogen->setName("H");
  TEST_STRING_SIMILAR(hydrogen->getName(), "H")
}
END_SECTION

START_SECTION((const name_type& getSequence() const ))
{
  TEST_STRING_SIMILAR(hydrogen->getSequence(), "H")
}
END_SECTION

START_SECTION((void setSequence(const name_type &sequence)))
{
  hydrogen->setSequence("H2");
  TEST_STRING_SIMILAR(hydrogen->getSequence(), "H2")
  hydrogen->setSequence("H");
  TEST_STRING_SIMILAR(hydrogen->getSequence(), "H")
}
END_SECTION

START_SECTION((nominal_mass_type getNominalMass() const ))
{
  TEST_EQUAL(hydrogen->getNominalMass(), iso.getNominalMass())
}
END_SECTION

START_SECTION((mass_type getMass(size_type index=0) const ))
{
  TEST_EQUAL(hydrogen->getMass(), 1.0078250319)
  TEST_EQUAL(hydrogen->getMass(0), 1.0078250319)
  TEST_EQUAL(hydrogen->getMass(1), 2.01410178)
  TEST_EQUAL(hydrogen->getMass(2), 3.01604927)
}
END_SECTION

START_SECTION((mass_type getAverageMass() const ))
{
  TEST_EQUAL(hydrogen->getAverageMass(), iso.getAverageMass())
}
END_SECTION

START_SECTION((mass_type getIonMass(int electrons_number=1) const ))
{
  double expected_ion_mass = hydrogen->getMass() - IMSElement::ELECTRON_MASS_IN_U;
  TEST_EQUAL(hydrogen->getIonMass(), expected_ion_mass)
  TEST_EQUAL(hydrogen->getIonMass(1), expected_ion_mass)

  expected_ion_mass = hydrogen->getMass() - 2 * IMSElement::ELECTRON_MASS_IN_U;
  TEST_EQUAL(hydrogen->getIonMass(2), expected_ion_mass)
}
END_SECTION

START_SECTION((const IMSIsotopeDistribution& getIsotopeDistribution() const ))
{
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), iso)
}
END_SECTION

START_SECTION((void setIsotopeDistribution(const IMSIsotopeDistribution &isotopes)))
{
  IMSIsotopeDistribution::peaks_container peaks_copy(peaks);
  peaks_copy.push_back(IMSIsotopeDistribution::peak_type(0.03604927,.0));
  IMSIsotopeDistribution modified_iso(peaks_copy, 1);
  hydrogen->setIsotopeDistribution(modified_iso);
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), modified_iso)
  hydrogen->setIsotopeDistribution(iso);
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), iso)
}
END_SECTION

START_SECTION((IMSElement& operator=(const IMSElement &element)))
{
  IMSElement hydrogen_copy;
  hydrogen_copy = *hydrogen;

  TEST_EQUAL(*hydrogen == hydrogen_copy, true)
  TEST_EQUAL(hydrogen->getAverageMass(), hydrogen_copy.getAverageMass())
  TEST_EQUAL(hydrogen->getIonMass(), hydrogen_copy.getIonMass())
  TEST_EQUAL(hydrogen->getIsotopeDistribution(), hydrogen_copy.getIsotopeDistribution())
  TEST_EQUAL(hydrogen->getMass(), hydrogen_copy.getMass())
  TEST_EQUAL(hydrogen->getName(), hydrogen_copy.getName())
  TEST_EQUAL(hydrogen->getNominalMass(), hydrogen_copy.getNominalMass())
  TEST_EQUAL(hydrogen->getSequence(), hydrogen_copy.getSequence())
}
END_SECTION

START_SECTION((bool operator==(const IMSElement &element) const ))
{
  IMSElement also_hydrogen("H",iso);

  TEST_EQUAL(*hydrogen == also_hydrogen, true)
  also_hydrogen.setName("D");
  TEST_EQUAL(*hydrogen == also_hydrogen, false)
  also_hydrogen.setName("H");
  also_hydrogen.setSequence("D");
  TEST_EQUAL(*hydrogen == also_hydrogen, false)
  also_hydrogen.setSequence("H");

  IMSIsotopeDistribution::peaks_container peaks_copy(peaks);
  peaks_copy.push_back(IMSIsotopeDistribution::peak_type(0.03604927,.0));
  IMSIsotopeDistribution modified_iso(peaks_copy, 1);

  also_hydrogen.setIsotopeDistribution(modified_iso);
  TEST_EQUAL(*hydrogen == also_hydrogen, false)

  IMSElement not_hydrogen;
  TEST_EQUAL(*hydrogen == not_hydrogen, false)
}
END_SECTION

START_SECTION((bool operator!=(const IMSElement &element) const ))
{
  IMSElement also_hydrogen("H",iso);
  TEST_EQUAL(*hydrogen != also_hydrogen, false)
  also_hydrogen.setName("D");
  TEST_EQUAL(*hydrogen != also_hydrogen, true)
  also_hydrogen.setName("H");
  also_hydrogen.setSequence("D");
  TEST_EQUAL(*hydrogen != also_hydrogen, true)
  also_hydrogen.setSequence("H");

  IMSIsotopeDistribution::peaks_container peaks_copy(peaks);
  peaks_copy.push_back(IMSIsotopeDistribution::peak_type(0.03604927,.0));
  IMSIsotopeDistribution modified_iso(peaks_copy, 1);

  also_hydrogen.setIsotopeDistribution(modified_iso);
  TEST_EQUAL(*hydrogen != also_hydrogen, true)


  IMSElement not_hydrogen;
  TEST_EQUAL(*hydrogen != not_hydrogen, true)
}
END_SECTION

delete hydrogen;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



