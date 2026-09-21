// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/NUXL/NuXLModificationsGenerator.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(NuXLModificationsGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NuXLModificationsGenerator* ptr = nullptr;
NuXLModificationsGenerator* null_ptr = nullptr;
START_SECTION(NuXLModificationsGenerator())
{
	ptr = new NuXLModificationsGenerator();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~NuXLModificationsGenerator())
{
	delete ptr;
}
END_SECTION

START_SECTION((static NuXLModificationMassesResult initModificationMassesNA(StringList target_nucleotides, StringList mappings, StringList restrictions, StringList modifications, std::string sequence_restriction, bool cysteine_adduct, Int max_length=4)))
{
  // TODO
}
END_SECTION

START_SECTION((static const ResidueModification* registerPrecursorAdduct(const std::string& nucleotide_composition, const EmpiricalFormula& adduct_formula, const Residue& residue)))
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  const EmpiricalFormula adduct("C9H11N2O8P1"); // U-H2O1
  const double adduct_mass = adduct.getMonoWeight();

  // a clean residue gets a definition named after the adduct, with the adduct's chemistry
  {
    AASequence seq = AASequence::fromString("AEADNLDDKK");
    const double unmodified_mass = seq.getMonoWeight();
    const EmpiricalFormula unmodified_formula = seq.getFormula();
    const Size before = db->getNumberOfModifications();

    const ResidueModification* mod = NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, seq[8]);
    TEST_TRUE(mod != nullptr)
    ABORT_IF(mod == nullptr)
    TEST_STRING_EQUAL(mod->getId(), "NuXL:U-H2O1")
    TEST_STRING_EQUAL(mod->getFullId(), "NuXL:U-H2O1 (K)")
    TEST_EQUAL(mod->getOrigin(), 'K')
    TEST_EQUAL(mod->getTermSpecificity(), ResidueModification::ANYWHERE)
    TEST_EQUAL(mod->getProvenance(), ResidueModification::DEFINED)
    TEST_STRING_EQUAL(mod->getDiffFormula().toString(), adduct.toString())
    TEST_REAL_SIMILAR(mod->getDiffMonoMass(), 306.0253048409)
    TEST_REAL_SIMILAR(mod->getDiffAverageMass(), adduct.getAverageWeight())
    TEST_EQUAL(db->getNumberOfModifications(), before + 1)
    TEST_TRUE(db->hasDefinedModification("NuXL:U-H2O1"))

    // registering the same adduct on the same residue again is a silent no-op
    std::ostringstream captured;
    OPENMS_LOG_WARN.insert(captured);
    const ResidueModification* again = NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, seq[8]);
    OPENMS_LOG_WARN.remove(captured);
    TEST_TRUE(again == mod)
    TEST_EQUAL(db->getNumberOfModifications(), before + 1)
    TEST_TRUE(captured.str().empty())

    // the same adduct on another residue is a separate entry
    AASequence tyr = AASequence::fromString("PEPTYDE");
    const ResidueModification* on_y = NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, tyr[4]);
    TEST_TRUE(on_y != nullptr)
    ABORT_IF(on_y == nullptr)
    TEST_TRUE(on_y != mod)
    TEST_STRING_EQUAL(on_y->getFullId(), "NuXL:U-H2O1 (Y)")
    TEST_EQUAL(db->getNumberOfModifications(), before + 2)

    // applied, it is written by name, re-parses to the same sequence and carries mass and formula
    seq.setModification(8, mod);
    TEST_STRING_EQUAL(seq.toString(), "AEADNLDDK(NuXL:U-H2O1)K")
    TEST_TRUE(AASequence::fromString(seq.toString()) == seq)
    TEST_REAL_SIMILAR(seq.getMonoWeight(), unmodified_mass + adduct_mass)
    TEST_STRING_EQUAL(seq.getFormula().toString(), (unmodified_formula + adduct).toString())
  }

  // an already modified residue: both modifications fold into one definition
  {
    AASequence seq = AASequence::fromString("PEPM(Oxidation)TIDE");
    const double oxidised_mass = seq.getMonoWeight();
    const EmpiricalFormula oxidised_formula = seq.getFormula();
    const double oxidation_mass = seq[3].getModification()->getDiffMonoMass();
    const Size before = db->getNumberOfModifications();

    std::ostringstream captured;
    OPENMS_LOG_WARN.insert(captured);
    const ResidueModification* combined = NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, seq[3]);
    OPENMS_LOG_WARN.remove(captured);
    TEST_TRUE(combined != nullptr)
    ABORT_IF(combined == nullptr)
    TEST_STRING_EQUAL(combined->getId(), "NuXL:U-H2O1~Oxidation")
    TEST_STRING_EQUAL(combined->getFullId(), "NuXL:U-H2O1~Oxidation (M)")
    TEST_EQUAL(combined->getOrigin(), 'M')
    TEST_STRING_EQUAL(combined->getDiffFormula().toString(), EmpiricalFormula("C9H11N2O9P1").toString())
    TEST_REAL_SIMILAR(combined->getDiffMonoMass(), adduct_mass + oxidation_mass)
    TEST_EQUAL(db->getNumberOfModifications(), before + 1)
    TEST_TRUE(captured.str().empty())

    // applying it replaces the oxidation, which the combined definition already contains
    seq.setModification(3, combined);
    TEST_STRING_EQUAL(seq.toString(), "PEPM(NuXL:U-H2O1~Oxidation)TIDE")
    TEST_TRUE(AASequence::fromString(seq.toString()) == seq)
    TEST_REAL_SIMILAR(seq.getMonoWeight(), oxidised_mass + adduct_mass)
    TEST_STRING_EQUAL(seq.getFormula().toString(), (oxidised_formula + adduct).toString())
  }

  // an existing modification whose formula contradicts its mass: the combination is defined by mass only
  {
    ResidueModification bogus;
    bogus.setId("TestNuXL:BogusFormula");
    bogus.setOrigin('M');
    bogus.setTermSpecificity(ResidueModification::ANYWHERE);
    bogus.setFullId();
    bogus.setDiffFormula(EmpiricalFormula("O"));
    bogus.setDiffMonoMass(100.0);
    bogus.setDiffAverageMass(100.0);
    db->registerDefinition(bogus);
    AASequence seq = AASequence::fromString("PEPM(TestNuXL:BogusFormula)TIDE");

    std::ostringstream captured;
    OPENMS_LOG_WARN.insert(captured);
    const ResidueModification* combined = NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, seq[3]);
    OPENMS_LOG_WARN.remove(captured);
    TEST_TRUE(combined != nullptr)
    ABORT_IF(combined == nullptr)
    TEST_STRING_EQUAL(combined->getId(), "NuXL:U-H2O1~TestNuXL:BogusFormula")
    TEST_TRUE(combined->getDiffFormula().isEmpty())
    TEST_REAL_SIMILAR(combined->getDiffMonoMass(), adduct_mass + 100.0)
    TEST_TRUE(captured.str().find("defined by mass only") != std::string::npos)
  }

  // no definition can be derived: anonymous mass shift, or a definition without a formula
  {
    AASequence anonymous = AASequence::fromString("PEPM[+12.3456]TIDE");
    TEST_TRUE(anonymous[3].isModified())
    TEST_TRUE(anonymous[3].getModification()->getId().empty())
    TEST_TRUE(NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, anonymous[3]) == nullptr)

    ResidueModification formula_less;
    formula_less.setId("TestNuXL:NoFormula");
    formula_less.setOrigin('M');
    formula_less.setTermSpecificity(ResidueModification::ANYWHERE);
    formula_less.setFullId();
    formula_less.setDiffMonoMass(42.0);
    db->registerDefinition(formula_less);
    AASequence named = AASequence::fromString("PEPM(TestNuXL:NoFormula)TIDE");
    TEST_TRUE(named[3].getModification()->getDiffFormula().isEmpty())
    TEST_TRUE(NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct, named[3]) == nullptr)
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



