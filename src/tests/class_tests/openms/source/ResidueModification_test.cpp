// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

// Modification tests
ResidueModification* ptr = nullptr;
ResidueModification* nullPointer = nullptr;
START_SECTION(ResidueModification())
  ptr = new ResidueModification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ResidueModification())
	delete ptr;
END_SECTION

ModificationsDB* mod_DB = ModificationsDB::getInstance();
ResidueDB* res_DB = ResidueDB::getInstance();

ptr = new ResidueModification();

START_SECTION(ResidueModification(const ResidueModification& modification))
  ResidueModification m(*ptr);
	TEST_EQUAL(m == *ptr, true)
END_SECTION

START_SECTION(ResidueModification& operator=(const ResidueModification& modification))
	ResidueModification m;
	m = *ptr;
	TEST_EQUAL(m == *ptr, true)
END_SECTION

START_SECTION(bool ResidueModification::operator<(const ResidueModification& rhs) const)
  ModificationsDB* ptr = ModificationsDB::getInstance();
  const ResidueModification* cm = ptr->getModification("Carboxymethyl (C)");
  const ResidueModification* pm = ptr->getModification("Phospho (S)");
  TEST_EQUAL(*cm < *pm, true);
END_SECTION

START_SECTION(void setId(const String& id))
	ptr->setId("blubb_new_id");
	TEST_STRING_EQUAL(ptr->getId(), "blubb_new_id")
END_SECTION

START_SECTION(const String& getId() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setFullName(const String& full_name))
	ptr->setFullName("blubb_new_full_name");
	TEST_STRING_EQUAL(ptr->getFullName(), "blubb_new_full_name")
END_SECTION

START_SECTION(const String& getFullName() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setName(const String& name))
	ptr->setName("blubb_new_name");
	TEST_STRING_EQUAL(ptr->getName(), "blubb_new_name")
END_SECTION

START_SECTION(const String& getName() const)
	NOT_TESTABLE
END_SECTION

START_SECTION((void setNeutralLossDiffFormula(const EmpiricalFormula& loss)))
  vector<EmpiricalFormula> loss_formulas;
  loss_formulas.push_back(EmpiricalFormula("H2O2"));
	ptr->setNeutralLossDiffFormulas(loss_formulas);
	TEST_EQUAL(ptr->getNeutralLossDiffFormulas()[0] == EmpiricalFormula("H2O2"), true)
END_SECTION

START_SECTION(const EmpiricalFormula& getNeutralLossDiffFormula() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setNeutralLossMonoMass(double mono_mass))
  vector<double> neutral_loss_mono_masses;
  neutral_loss_mono_masses.push_back(123.345678);
	ptr->setNeutralLossMonoMasses(neutral_loss_mono_masses);
	TEST_REAL_SIMILAR(ptr->getNeutralLossMonoMasses()[0], 123.345678);
END_SECTION

START_SECTION((double getNeutralLossMonoMass() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setNeutralLossAverageMass(double average_mass)))
  vector<double> neutral_loss_avg_masses;
  neutral_loss_avg_masses.push_back(23.345678);
	ptr->setNeutralLossAverageMasses(neutral_loss_avg_masses);
	TEST_REAL_SIMILAR(ptr->getNeutralLossAverageMasses()[0], 23.345678)
END_SECTION

START_SECTION(double getNeutralLossAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION((bool hasNeutralLoss() const))
	TEST_EQUAL(ptr->hasNeutralLoss(), true)
	ResidueModification mod;
	TEST_EQUAL(mod.hasNeutralLoss(), false)
  vector<EmpiricalFormula> loss_formulas;
  loss_formulas.push_back(EmpiricalFormula("H2O"));
	mod.setNeutralLossDiffFormulas(loss_formulas);
	TEST_EQUAL(mod.hasNeutralLoss(), true)
END_SECTION

START_SECTION((void setFullId(const String& full_id)))
	ptr->setFullId("blubb_new_fullid");
	TEST_STRING_EQUAL(ptr->getFullId(), "blubb_new_fullid")
END_SECTION

START_SECTION((const String& getFullId() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setUniModRecordId(const Int& id)))
	ptr->setUniModRecordId(42);
	TEST_EQUAL(ptr->getUniModRecordId(), 42)
END_SECTION

START_SECTION((const String& getUniModRecordId() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((const String& getUniModAccession() const))
	ptr->setUniModRecordId(42);
	TEST_STRING_EQUAL(ptr->getUniModAccession(), "UniMod:42")
END_SECTION

START_SECTION((void setPSIMODAccession(const String& id)))
	ptr->setPSIMODAccession("blubb_new_PSIMODAccession");
	TEST_STRING_EQUAL(ptr->getPSIMODAccession(), "blubb_new_PSIMODAccession")
END_SECTION

START_SECTION((const String& getPSIMODAccession() const))
	NOT_TESTABLE
END_SECTION


START_SECTION(void setTermSpecificity(TermSpecificity term_spec))
	ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
END_SECTION

START_SECTION(void setTermSpecificity(const String& name))
	ptr->setTermSpecificity("C-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity("N-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
	ptr->setTermSpecificity("none");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
END_SECTION

START_SECTION(TermSpecificity getTermSpecificity() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(String getTermSpecificityName(TermSpecificity=NUMBER_OF_TERM_SPECIFICITY) const)
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "C-term")
	ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "N-term")
	ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "none")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::C_TERM), "C-term")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::N_TERM), "N-term")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::ANYWHERE), "none")
END_SECTION

START_SECTION(void setOrigin(char origin))
	ptr->setOrigin('A');
	TEST_EQUAL(ptr->getOrigin(), 'A')
END_SECTION

START_SECTION(char getOrigin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setSourceClassification(SourceClassification classification))
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
END_SECTION

START_SECTION(void setSourceClassification(const String& classification))
	ptr->setSourceClassification("Artifact");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification("Natural");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification("Hypothetical");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
END_SECTION

START_SECTION(SourceClassification getSourceClassification() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(String getSourceClassificationName(SourceClassification classification=NUMBER_OF_SOURCE_CLASSIFICATIONS) const)
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Artefact")
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Natural")
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Hypothetical")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::ARTIFACT), "Artefact")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::NATURAL), "Natural")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::HYPOTHETICAL), "Hypothetical")
END_SECTION

START_SECTION(void setAverageMass(double mass))
	ptr->setAverageMass(2.0);
	TEST_REAL_SIMILAR(ptr->getAverageMass(), 2.0)
END_SECTION

START_SECTION(double getAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setMonoMass(double mass))
	ptr->setMonoMass(3.0);
	TEST_REAL_SIMILAR(ptr->getMonoMass(), 3.0)
END_SECTION

START_SECTION(double getMonoMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDiffAverageMass(double mass))
	ptr->setDiffAverageMass(4.0);
	TEST_REAL_SIMILAR(ptr->getDiffAverageMass(), 4.0)
END_SECTION

START_SECTION(double getDiffAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDiffMonoMass(double mass))
	ptr->setDiffMonoMass(5.0);
	TEST_REAL_SIMILAR(ptr->getDiffMonoMass(), 5.0)
END_SECTION

START_SECTION(double getDiffMonoMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setFormula(const String& composition))
	ptr->setFormula("blubb_new_formula");
	TEST_STRING_EQUAL(ptr->getFormula(), "blubb_new_formula")
END_SECTION

START_SECTION(const String& getFormula() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDiffFormula(const EmpiricalFormula& diff_formula))
	EmpiricalFormula ef("C3H4S-3");
	ptr->setDiffFormula(ef);
	TEST_EQUAL(ptr->getDiffFormula() == ef, true)
END_SECTION

START_SECTION(const EmpiricalFormula& getDiffFormula() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setSynonyms(const std::set<String>& synonyms))
	set<String> synonyms;
	synonyms.insert("blubb_syn1");
	synonyms.insert("blubb_syn2");
	ptr->setSynonyms(synonyms);
	TEST_EQUAL(ptr->getSynonyms() == synonyms, true)
END_SECTION

START_SECTION(void addSynonym(const String& synonym))
	ptr->addSynonym("blubb_syn3");
	TEST_EQUAL(ptr->getSynonyms().size(), 3)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(bool operator==(const ResidueModification& modification) const)
	ResidueModification mod1, mod2;
	mod1.setId("Id");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setId("Id");
	TEST_TRUE(mod1 == mod2)

	mod1.setFullName("FullName");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setFullName("FullName");
	TEST_TRUE(mod1 == mod2)

	mod1.setName("Name");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setName("Name");
	TEST_TRUE(mod1 == mod2)

	mod1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setTermSpecificity(ResidueModification::N_TERM);
	TEST_TRUE(mod1 == mod2)

	mod1.setOrigin('C');
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setOrigin('C');
	TEST_TRUE(mod1 == mod2)

	mod1.setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setSourceClassification(ResidueModification::NATURAL);
	TEST_TRUE(mod1 == mod2)

	mod1.setAverageMass(0.123);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setAverageMass(0.123);
	TEST_TRUE(mod1 == mod2)

	mod1.setMonoMass(1.23);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setMonoMass(1.23);
	TEST_TRUE(mod1 == mod2)

	mod1.setDiffAverageMass(2.34);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffAverageMass(2.34);
	TEST_TRUE(mod1 == mod2)

	mod1.setDiffMonoMass(3.45);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffMonoMass(3.45);
	TEST_TRUE(mod1 == mod2)

	mod1.setFormula("C 3 H 4");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setFormula("C 3 H 4");
	TEST_TRUE(mod1 == mod2)

	mod1.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
	TEST_TRUE(mod1 == mod2)

	mod1.addSynonym("new_syn");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.addSynonym("new_syn");
	TEST_TRUE(mod1 == mod2)
END_SECTION

START_SECTION(bool operator!=(const ResidueModification& modification) const)
	ResidueModification mod1, mod2;
  mod1.setId("Id");
  TEST_FALSE(mod1 == mod2)
  mod2.setId("Id");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setFullName("FullName");
  TEST_FALSE(mod1 == mod2)
  mod2.setFullName("FullName");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setName("Name");
  TEST_FALSE(mod1 == mod2)
  mod2.setName("Name");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setTermSpecificity(ResidueModification::N_TERM);
  TEST_FALSE(mod1 == mod2)
  mod2.setTermSpecificity(ResidueModification::N_TERM);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setOrigin('C');
  TEST_FALSE(mod1 == mod2)
  mod2.setOrigin('C');
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setSourceClassification(ResidueModification::NATURAL);
  TEST_FALSE(mod1 == mod2)
  mod2.setSourceClassification(ResidueModification::NATURAL);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setAverageMass(0.123);
  TEST_FALSE(mod1 == mod2)
  mod2.setAverageMass(0.123);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setMonoMass(1.23);
  TEST_FALSE(mod1 == mod2)
  mod2.setMonoMass(1.23);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffAverageMass(2.34);
  TEST_FALSE(mod1 == mod2)
  mod2.setDiffAverageMass(2.34);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffMonoMass(3.45);
  TEST_FALSE(mod1 == mod2)
  mod2.setDiffMonoMass(3.45);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setFormula("C 3 H 4");
  TEST_FALSE(mod1 == mod2)
  mod2.setFormula("C 3 H 4");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
  TEST_FALSE(mod1 == mod2)
  mod2.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
  TEST_EQUAL(mod1 != mod2, false)

  mod1.addSynonym("new_syn");
  TEST_FALSE(mod1 == mod2)
  mod2.addSynonym("new_syn");
  TEST_EQUAL(mod1 != mod2, false)
END_SECTION

const ResidueModification* combined_mod;

START_SECTION(static const ResidueModification* combineMods(const ResidueModification* base,
	const std::set<const ResidueModification*>& addons,
	bool allow_unknown_masses = false,
	const Residue* residue = nullptr))


  const ResidueModification* base = mod_DB->getModification("Phospho (S)");
	std::set<const ResidueModification*> addons;
	addons.insert(mod_DB->getModification("Label:15N(1) (S)"));
	// boring case: 1 base + 1 addon
	auto res = ResidueModification::combineMods(base, addons, false, res_DB->getResidue('S'));
	TEST_EQUAL(res == nullptr, false);
	TEST_EQUAL(res->getOrigin(), 'S');
	TEST_EQUAL(res->getTermSpecificity(), ResidueModification::ANYWHERE);
	TEST_REAL_SIMILAR(res->getDiffMonoMass(), base->getDiffMonoMass() + (*addons.begin())->getDiffMonoMass());
	combined_mod = res; // useful for downstream tests

	// test empty addons
	TEST_EQUAL(base, ResidueModification::combineMods(base, {}, false, res_DB->getResidue('S')));

	// test empty base + 1 addon
	res = ResidueModification::combineMods(nullptr, addons, false, res_DB->getResidue('S'));
	TEST_EQUAL(res, *addons.begin());

	// 1 base, 2 addons (1 is invalid)
	addons.insert(mod_DB->getModification("Label:15N(1) (T)"));
	TEST_EXCEPTION(Exception::Precondition, ResidueModification::combineMods(base, addons, false, res_DB->getResidue('S'));)

	// both empty
	TEST_EQUAL(ResidueModification::combineMods(nullptr, {}, true) == nullptr, true)
END_SECTION

START_SECTION(String toString() const)
	const ResidueModification* base = mod_DB->getModification("Phospho (S)");
	TEST_EQUAL(base->toString(), "S(Phospho)")
  TEST_EQUAL(combined_mod->toString(), "S[+80.963365999999994]")

END_SECTION

START_SECTION(static String getDiffMonoMassString(const double diff_mono_mass))
	TEST_EQUAL(ResidueModification::getDiffMonoMassString(16), "+16.0");
	TEST_EQUAL(ResidueModification::getDiffMonoMassString(-16), "-16.0");
END_SECTION

/// return a string of the form '[+&gt;mass&lt;] (the '+' might be a '-', if mass is negative).
START_SECTION(static String getDiffMonoMassWithBracket(const double diff_mono_mass))
	TEST_EQUAL(ResidueModification::getDiffMonoMassWithBracket(16), "[+16.0]");
	TEST_EQUAL(ResidueModification::getDiffMonoMassWithBracket(-16), "[-16.0]");
END_SECTION

/// return a string of the form '[&gt;mass&lt;]
START_SECTION(static String getMonoMassWithBracket(const double mono_mass))
	TEST_EQUAL(ResidueModification::getDiffMonoMassWithBracket(16), "[+16.0]");
	TEST_EXCEPTION(Exception::InvalidValue, ResidueModification::getMonoMassWithBracket(-16));
END_SECTION

delete ptr;

END_TEST
