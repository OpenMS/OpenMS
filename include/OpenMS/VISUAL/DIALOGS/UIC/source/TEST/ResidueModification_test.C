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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id: ResidueModification_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////

// Modification tests
ResidueModification* ptr = 0;
START_SECTION(ResidueModification())
  ptr = new ResidueModification();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~ResidueModification())
	delete ptr;
END_SECTION

ptr = new ResidueModification();

START_SECTION(ResidueModification(const ResidueModification& modification))
  ResidueModification m(*ptr);
	TEST_EQUAL(m == *ptr, true)
END_SECTION

START_SECTION(ResidueModification& operator = (const ResidueModification& modification))
	ResidueModification m;
	m = *ptr;
	TEST_EQUAL(m == *ptr, true)
END_SECTION

START_SECTION(void setId(const String &id))
	ptr->setId("blubb_new_id");
	TEST_STRING_EQUAL(ptr->getId(), "blubb_new_id")
END_SECTION

START_SECTION(const String& getId() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setFullName(const String &full_name))
	ptr->setFullName("blubb_new_full_name");
	TEST_STRING_EQUAL(ptr->getFullName(), "blubb_new_full_name")
END_SECTION

START_SECTION(const String& getFullName() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setName(const String &name))
	ptr->setName("blubb_new_name");
	TEST_STRING_EQUAL(ptr->getName(), "blubb_new_name")
END_SECTION

START_SECTION(const String& getName() const)
	NOT_TESTABLE
END_SECTION

START_SECTION((void setNeutralLossDiffFormula(const EmpiricalFormula& loss)))
	ptr->setNeutralLossDiffFormula(EmpiricalFormula("H2O2"));
	TEST_EQUAL(ptr->getNeutralLossDiffFormula() == EmpiricalFormula("H2O2"), true)
END_SECTION

START_SECTION(const EmpiricalFormula& getNeutralLossDiffFormula() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setNeutralLossMonoMass(DoubleReal mono_mass))
	ptr->setNeutralLossMonoMass(123.345678);
	TEST_REAL_SIMILAR(ptr->getNeutralLossMonoMass(), 123.345678);
END_SECTION

START_SECTION((DoubleReal getNeutralLossMonoMass() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setNeutralLossAverageMass(DoubleReal average_mass)))
	ptr->setNeutralLossAverageMass(23.345678);
	TEST_REAL_SIMILAR(ptr->getNeutralLossAverageMass(), 23.345678)
END_SECTION

START_SECTION(DoubleReal getNeutralLossAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION((bool hasNeutralLoss() const))
	TEST_EQUAL(ptr->hasNeutralLoss(), true)
	ResidueModification mod;
	TEST_EQUAL(mod.hasNeutralLoss(), false)
	mod.setNeutralLossDiffFormula(EmpiricalFormula("H2O"));
	TEST_EQUAL(mod.hasNeutralLoss(), true)
END_SECTION

START_SECTION((void setFullId(const String& full_id)))
	ptr->setFullId("blubb_new_fullid");
	TEST_STRING_EQUAL(ptr->getFullId(), "blubb_new_fullid")
END_SECTION

START_SECTION((const String& getFullId() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setUniModAccession(const String &id)))
	ptr->setUniModAccession("blubb_new_UniModAccession");
	TEST_STRING_EQUAL(ptr->getUniModAccession(), "blubb_new_UniModAccession")
END_SECTION

START_SECTION((const String& getUniModAccession() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setPSIMODAccession(const String& id)))
	ptr->setPSIMODAccession("blubb_new_PSIMODAccession");
	TEST_STRING_EQUAL(ptr->getPSIMODAccession(), "blubb_new_PSIMODAccession")
END_SECTION

START_SECTION((const String& getPSIMODAccession() const))
	NOT_TESTABLE
END_SECTION


START_SECTION(void setTermSpecificity(Term_Specificity term_spec))
	ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
END_SECTION

START_SECTION(void setTermSpecificity(const String &name))
	ptr->setTermSpecificity("C-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity("N-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
	ptr->setTermSpecificity("none");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
END_SECTION

START_SECTION(Term_Specificity getTermSpecificity() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(String getTermSpecificityName(Term_Specificity=NUMBER_OF_TERM_SPECIFICITY) const)
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

START_SECTION(void setOrigin(const String &origin))
	ptr->setOrigin("blubb_new_origin");
	TEST_STRING_EQUAL(ptr->getOrigin(), "blubb_new_origin")
END_SECTION

START_SECTION(const String& getOrigin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setSourceClassification(Source_Classification classification))
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
END_SECTION

START_SECTION(void setSourceClassification(const String &classification))
	ptr->setSourceClassification("Artifact");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification("Natural");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification("Hypothetical");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
END_SECTION

START_SECTION(Source_Classification getSourceClassification() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(String getSourceClassificationName(Source_Classification classification=NUMBER_OF_SOURCE_CLASSIFICATIONS) const)
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Artifact")
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Natural")
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Hypothetical")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::ARTIFACT), "Artifact")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::NATURAL), "Natural")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::HYPOTHETICAL), "Hypothetical")
END_SECTION

START_SECTION(void setAverageMass(DoubleReal mass))
	ptr->setAverageMass(2.0);
	TEST_REAL_SIMILAR(ptr->getAverageMass(), 2.0)
END_SECTION

START_SECTION(DoubleReal getAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setMonoMass(DoubleReal mass))
	ptr->setMonoMass(3.0);
	TEST_REAL_SIMILAR(ptr->getMonoMass(), 3.0)
END_SECTION

START_SECTION(DoubleReal getMonoMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDiffAverageMass(DoubleReal mass))
	ptr->setDiffAverageMass(4.0);
	TEST_REAL_SIMILAR(ptr->getDiffAverageMass(), 4.0)
END_SECTION

START_SECTION(DoubleReal getDiffAverageMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDiffMonoMass(DoubleReal mass))
	ptr->setDiffMonoMass(5.0);
	TEST_REAL_SIMILAR(ptr->getDiffMonoMass(), 5.0)
END_SECTION

START_SECTION(DoubleReal getDiffMonoMass() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setFormula(const String &composition))
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

START_SECTION(void setSynonyms(const std::set< String > &synonyms))
	set<String> synonyms;
	synonyms.insert("blubb_syn1");
	synonyms.insert("blubb_syn2");
	ptr->setSynonyms(synonyms);
	TEST_EQUAL(ptr->getSynonyms() == synonyms, true)
END_SECTION

START_SECTION(void addSynonym(const String &synonym))
	ptr->addSynonym("blubb_syn3");
	TEST_EQUAL(ptr->getSynonyms().size(), 3)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(bool operator==(const ResidueModification &modification) const)
	ResidueModification mod1, mod2;
	mod1.setId("Id");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setId("Id");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setFullName("FullName");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setFullName("FullName");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setName("Name");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setName("Name");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setOrigin("C");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setOrigin("C");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setAverageMass(0.123);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setAverageMass(0.123);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setMonoMass(1.23);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setMonoMass(1.23);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setDiffAverageMass(2.34);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffAverageMass(2.34);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setDiffMonoMass(3.45);
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffMonoMass(3.45);
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setFormula("C 3 H 4");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setFormula("C 3 H 4");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
	TEST_EQUAL(mod1 == mod2, true)

	mod1.addSynonym("new_syn");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.addSynonym("new_syn");
	TEST_EQUAL(mod1 == mod2, true)
END_SECTION

START_SECTION(bool operator!=(const ResidueModification &modification) const)
	ResidueModification mod1, mod2;
  mod1.setId("Id");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setId("Id");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setFullName("FullName");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setFullName("FullName");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setName("Name");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setName("Name");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setTermSpecificity(ResidueModification::N_TERM);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setTermSpecificity(ResidueModification::N_TERM);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setOrigin("C");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setOrigin("C");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setSourceClassification(ResidueModification::NATURAL);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setSourceClassification(ResidueModification::NATURAL);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setAverageMass(0.123);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setAverageMass(0.123);
  TEST_EQUAL(mod1 != mod2, false)


  mod1.setMonoMass(1.23);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setMonoMass(1.23);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffAverageMass(2.34);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setDiffAverageMass(2.34);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffMonoMass(3.45);
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setDiffMonoMass(3.45);
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setFormula("C 3 H 4");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setFormula("C 3 H 4");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setDiffFormula(EmpiricalFormula("C0H-2N0O0"));
  TEST_EQUAL(mod1 != mod2, false)

  mod1.addSynonym("new_syn");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.addSynonym("new_syn");
  TEST_EQUAL(mod1 != mod2, false)
END_SECTION


END_TEST
