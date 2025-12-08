// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Isotope.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Isotope, "$Id$")

/////////////////////////////////////////////////////////////

Isotope* iso_ptr = nullptr;
Isotope* iso_nullPointer = nullptr;

START_SECTION(Isotope())
  iso_ptr = new Isotope;
  TEST_NOT_EQUAL(iso_ptr, iso_nullPointer)
END_SECTION

START_SECTION(~Isotope())
  delete iso_ptr;
END_SECTION

START_SECTION((Isotope(const std::string& name, const std::string& symbol, unsigned int atomic_number, unsigned int mass_number, double atomic_mass, double abundance, double half_life, DecayMode decay_mode)))
  iso_ptr = new Isotope("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  TEST_NOT_EQUAL(iso_ptr, iso_nullPointer)
  TEST_EQUAL(iso_ptr->getName(), "(14)Carbon")
  TEST_EQUAL(iso_ptr->getSymbol(), "(14)C")
  TEST_EQUAL(iso_ptr->getAtomicNumber(), 6)
  TEST_EQUAL(iso_ptr->getMassNumber(), 14)
  TEST_REAL_SIMILAR(iso_ptr->getMonoWeight(), 14.003241989)
  TEST_REAL_SIMILAR(iso_ptr->getAbundance(), 0.0)
  TEST_REAL_SIMILAR(iso_ptr->getHalfLife(), 1.807e11)
  TEST_EQUAL(iso_ptr->getDecayMode(), Isotope::DecayMode::BETA_MINUS)
  TEST_EQUAL(iso_ptr->isStable(), false)
  delete iso_ptr;
END_SECTION

START_SECTION(Isotope(const Isotope& isotope))
  Isotope original("(12)Carbon", "(12)C", 6, 12, 12.0, 0.9893, -1.0, Isotope::DecayMode::STABLE);
  Isotope copy(original);
  TEST_EQUAL(copy == original, true)
END_SECTION

iso_ptr = new Isotope("(235)Uranium", "(235)U", 92, 235, 235.043928, 0.007204, 2.22e16, Isotope::DecayMode::ALPHA);

START_SECTION(void setMassNumber(unsigned int mass_number))
  iso_ptr->setMassNumber(238);
  NOT_TESTABLE
END_SECTION

START_SECTION(unsigned int getMassNumber() const)
  TEST_EQUAL(iso_ptr->getMassNumber(), 238)
END_SECTION

START_SECTION(unsigned int getNeutrons() const)
  // Neutrons = mass_number - atomic_number = 238 - 92 = 146
  TEST_EQUAL(iso_ptr->getNeutrons(), 146)
END_SECTION

START_SECTION(void setAbundance(double abundance))
  iso_ptr->setAbundance(0.5);
  NOT_TESTABLE
END_SECTION

START_SECTION(double getAbundance() const)
  TEST_REAL_SIMILAR(iso_ptr->getAbundance(), 0.5)
END_SECTION

START_SECTION(void setHalfLife(double half_life))
  iso_ptr->setHalfLife(1.0e10);
  NOT_TESTABLE
END_SECTION

START_SECTION(double getHalfLife() const)
  TEST_REAL_SIMILAR(iso_ptr->getHalfLife(), 1.0e10)
END_SECTION

START_SECTION(void setDecayMode(DecayMode decay_mode))
  iso_ptr->setDecayMode(Isotope::DecayMode::BETA_MINUS);
  NOT_TESTABLE
END_SECTION

START_SECTION(DecayMode getDecayMode() const)
  TEST_EQUAL(iso_ptr->getDecayMode(), Isotope::DecayMode::BETA_MINUS)
END_SECTION

START_SECTION(bool isStable() const)
  // With positive half-life, not stable
  TEST_EQUAL(iso_ptr->isStable(), false)

  // With negative half-life, stable
  iso_ptr->setHalfLife(-1.0);
  TEST_EQUAL(iso_ptr->isStable(), true)
END_SECTION

START_SECTION(virtual bool isIsotope() const)
  TEST_EQUAL(iso_ptr->isIsotope(), true)
END_SECTION

START_SECTION(Isotope& operator=(const Isotope& isotope))
  Isotope original("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  Isotope copy;
  copy = original;
  TEST_EQUAL(copy == original, true)
END_SECTION

START_SECTION(bool operator==(const Isotope& isotope) const)
  Isotope iso1("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  Isotope iso2("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  TEST_EQUAL(iso1 == iso2, true)

  iso2.setMassNumber(15);
  TEST_EQUAL(iso1 == iso2, false)
END_SECTION

START_SECTION(bool operator!=(const Isotope& isotope) const)
  Isotope iso1("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  Isotope iso2("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  TEST_EQUAL(iso1 != iso2, false)

  iso2.setMassNumber(15);
  TEST_EQUAL(iso1 != iso2, true)
END_SECTION

START_SECTION(bool operator<(const Isotope& isotope) const)
  Isotope iso1("(12)Carbon", "(12)C", 6, 12, 12.0, 0.9893, -1.0, Isotope::DecayMode::STABLE);
  Isotope iso2("(14)Carbon", "(14)C", 6, 14, 14.003241989, 0.0, 1.807e11, Isotope::DecayMode::BETA_MINUS);
  TEST_EQUAL(iso1 < iso2, true)
  TEST_EQUAL(iso2 < iso1, false)
END_SECTION

START_SECTION(static std::string decayModeToString(DecayMode mode))
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::STABLE), "STABLE")
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::ALPHA), "ALPHA")
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::BETA_PLUS), "BETA_PLUS")
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::BETA_MINUS), "BETA_MINUS")
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::PROTON), "PROTON")
  TEST_EQUAL(Isotope::decayModeToString(Isotope::DecayMode::UNKNOWN), "UNKNOWN")
END_SECTION

START_SECTION(static DecayMode stringToDecayMode(const std::string& mode_str))
  TEST_EQUAL(Isotope::stringToDecayMode("STABLE"), Isotope::DecayMode::STABLE)
  TEST_EQUAL(Isotope::stringToDecayMode("ALPHA"), Isotope::DecayMode::ALPHA)
  TEST_EQUAL(Isotope::stringToDecayMode("BETA_PLUS"), Isotope::DecayMode::BETA_PLUS)
  TEST_EQUAL(Isotope::stringToDecayMode("BETA_MINUS"), Isotope::DecayMode::BETA_MINUS)
  TEST_EQUAL(Isotope::stringToDecayMode("PROTON"), Isotope::DecayMode::PROTON)
  TEST_EQUAL(Isotope::stringToDecayMode("UNKNOWN"), Isotope::DecayMode::UNKNOWN)
  TEST_EQUAL(Isotope::stringToDecayMode("invalid"), Isotope::DecayMode::UNKNOWN)
END_SECTION

delete iso_ptr;

// Test ElementDB isotope retrieval
START_SECTION([EXTRA] ElementDB::getIsotope)
{
  ElementDB* db = ElementDB::getInstance();

  // Test Carbon-14
  const Isotope* c14 = db->getIsotope("C", 14);
  TEST_NOT_EQUAL(c14, iso_nullPointer)
  if (c14 != nullptr)
  {
    TEST_EQUAL(c14->getMassNumber(), 14)
    TEST_EQUAL(c14->getAtomicNumber(), 6)
    TEST_EQUAL(c14->isStable(), false)
    TEST_EQUAL(c14->getDecayMode(), Isotope::DecayMode::BETA_MINUS)
  }

  // Test Carbon-12 (stable)
  const Isotope* c12 = db->getIsotope("C", 12);
  TEST_NOT_EQUAL(c12, iso_nullPointer)
  if (c12 != nullptr)
  {
    TEST_EQUAL(c12->getMassNumber(), 12)
    TEST_EQUAL(c12->isStable(), true)
    TEST_EQUAL(c12->getDecayMode(), Isotope::DecayMode::STABLE)
  }

  // Test Uranium-235
  const Isotope* u235 = db->getIsotope("U", 235);
  TEST_NOT_EQUAL(u235, iso_nullPointer)
  if (u235 != nullptr)
  {
    TEST_EQUAL(u235->getMassNumber(), 235)
    TEST_EQUAL(u235->getAtomicNumber(), 92)
    TEST_EQUAL(u235->isStable(), false)
    TEST_EQUAL(u235->getDecayMode(), Isotope::DecayMode::ALPHA)
  }

  // Test Uranium-238
  const Isotope* u238 = db->getIsotope("U", 238);
  TEST_NOT_EQUAL(u238, iso_nullPointer)
  if (u238 != nullptr)
  {
    TEST_EQUAL(u238->getMassNumber(), 238)
    TEST_EQUAL(u238->getDecayMode(), Isotope::DecayMode::ALPHA)
  }

  // Test retrieval by atomic number
  const Isotope* c14_by_an = db->getIsotope(6, 14);
  TEST_EQUAL(c14_by_an, c14)

  // Test non-existent isotope
  const Isotope* non_existent = db->getIsotope("C", 100);
  TEST_EQUAL(non_existent, iso_nullPointer)
}
END_SECTION

START_SECTION([EXTRA] ElementDB::hasIsotope)
{
  ElementDB* db = ElementDB::getInstance();

  TEST_EQUAL(db->hasIsotope("C", 12), true)
  TEST_EQUAL(db->hasIsotope("C", 13), true)
  TEST_EQUAL(db->hasIsotope("C", 14), true)
  TEST_EQUAL(db->hasIsotope("C", 100), false)

  TEST_EQUAL(db->hasIsotope("U", 234), true)
  TEST_EQUAL(db->hasIsotope("U", 235), true)
  TEST_EQUAL(db->hasIsotope("U", 238), true)

  // Test by atomic number
  TEST_EQUAL(db->hasIsotope(6, 14), true)
  TEST_EQUAL(db->hasIsotope(92, 235), true)
}
END_SECTION

START_SECTION([EXTRA] ElementDB::getIsotopeSymbols)
{
  ElementDB* db = ElementDB::getInstance();

  const auto& iso_symbols = db->getIsotopeSymbols();
  TEST_EQUAL(iso_symbols.count("(14)C") > 0, true)
  TEST_EQUAL(iso_symbols.count("(235)U") > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
