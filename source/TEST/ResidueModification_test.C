// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

// Modification tests
ResidueModification* ptr = 0;
CHECK(ResidueModification())
  ptr = new ResidueModification();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ResidueModification())
	delete ptr;
RESULT

ptr = new ResidueModification();

CHECK(ResidueModification(const ResidueModification& modification))
  ResidueModification m(*ptr);
	TEST_EQUAL(m == *ptr, true)
RESULT

CHECK(ResidueModification& operator = (const ResidueModification& modification))
	ResidueModification m;
	m = *ptr;
	TEST_EQUAL(m == *ptr, true)
RESULT

CHECK(void setId(const String &id))
	ptr->setId("blubb_new_id");
	TEST_STRING_EQUAL(ptr->getId(), "blubb_new_id")
RESULT

CHECK(const String& getId() const)
	NOT_TESTABLE
RESULT

CHECK(void setFullName(const String &full_name))
	ptr->setFullName("blubb_new_full_name");
	TEST_STRING_EQUAL(ptr->getFullName(), "blubb_new_full_name")
RESULT

CHECK(const String& getFullName() const)
	NOT_TESTABLE
RESULT

CHECK(void setName(const String &name))
	ptr->setName("blubb_new_name");
	TEST_STRING_EQUAL(ptr->getName(), "blubb_new_name")
RESULT

CHECK(const String& getName() const)
	NOT_TESTABLE
RESULT

CHECK(void setTermSpecificity(Term_Specificity term_spec))
	ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
RESULT

CHECK(void setTermSpecificity(const String &name))
	ptr->setTermSpecificity("C-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM)
	ptr->setTermSpecificity("N-term");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
	ptr->setTermSpecificity("none");
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE)
RESULT

CHECK(Term_Specificity getTermSpecificity() const)
	NOT_TESTABLE
RESULT

CHECK(String getTermSpecificityName(Term_Specificity=NUMBER_OF_TERM_SPECIFICITY) const)
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "C-term")
	ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "N-term")
	ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(), "none")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::C_TERM), "C-term")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::N_TERM), "N-term")
	TEST_STRING_EQUAL(ptr->getTermSpecificityName(ResidueModification::ANYWHERE), "none")
RESULT

CHECK(void setOrigin(const String &origin))
	ptr->setOrigin("blubb_new_origin");
	TEST_STRING_EQUAL(ptr->getOrigin(), "blubb_new_origin")
RESULT

CHECK(const String& getOrigin() const)
	NOT_TESTABLE
RESULT

CHECK(void setSourceClassification(Source_Classification classification))
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
RESULT

CHECK(void setSourceClassification(const String &classification))
	ptr->setSourceClassification("Artifact");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::ARTIFACT)
	ptr->setSourceClassification("Natural");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::NATURAL)
	ptr->setSourceClassification("Hypothetical");
	TEST_EQUAL(ptr->getSourceClassification(), ResidueModification::HYPOTHETICAL)
RESULT

CHECK(Source_Classification getSourceClassification() const)
	NOT_TESTABLE
RESULT

CHECK(String getSourceClassificationName(Source_Classification classification=NUMBER_OF_SOURCE_CLASSIFICATIONS) const)
	ptr->setSourceClassification(ResidueModification::ARTIFACT);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Artifact")
	ptr->setSourceClassification(ResidueModification::NATURAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Natural")
	ptr->setSourceClassification(ResidueModification::HYPOTHETICAL);
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(), "Hypothetical")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::ARTIFACT), "Artifact")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::NATURAL), "Natural")
	TEST_STRING_EQUAL(ptr->getSourceClassificationName(ResidueModification::HYPOTHETICAL), "Hypothetical")
RESULT

CHECK(void setAverageMass(double mass))
	ptr->setAverageMass(2.0);
	TEST_REAL_EQUAL(ptr->getAverageMass(), 2.0)
RESULT

CHECK(double getAverageMass() const)
	NOT_TESTABLE
RESULT

CHECK(void setMonoMass(double mass))
	ptr->setMonoMass(3.0);
	TEST_REAL_EQUAL(ptr->getMonoMass(), 3.0)
RESULT

CHECK(double getMonoMass() const)
	NOT_TESTABLE
RESULT

CHECK(void setDiffAverageMass(double mass))
	ptr->setDiffAverageMass(4.0);
	TEST_REAL_EQUAL(ptr->getDiffAverageMass(), 4.0)
RESULT

CHECK(double getDiffAverageMass() const)
	NOT_TESTABLE
RESULT

CHECK(void setDiffMonoMass(double mass))
	ptr->setDiffMonoMass(5.0);
	TEST_REAL_EQUAL(ptr->getDiffMonoMass(), 5.0)
RESULT

CHECK(double getDiffMonoMass() const)
	NOT_TESTABLE
RESULT

CHECK(void setFormula(const String &composition))
	ptr->setFormula("blubb_new_formula");
	TEST_STRING_EQUAL(ptr->getFormula(), "blubb_new_formula")
RESULT

CHECK(const String& getFormula() const)
	NOT_TESTABLE
RESULT

CHECK(void setDiffFormula(const String &diff_formula))
	ptr->setDiffFormula("blubb_new_diff_formula");
	TEST_STRING_EQUAL(ptr->getDiffFormula(), "blubb_new_diff_formula")
RESULT

CHECK(const String& getDiffFormula() const)
	NOT_TESTABLE
RESULT

CHECK(void setSynonyms(const std::set< String > &synonyms))
	set<String> synonyms;
	synonyms.insert("blubb_syn1");
	synonyms.insert("blubb_syn2");
	ptr->setSynonyms(synonyms);
	TEST_EQUAL(ptr->getSynonyms() == synonyms, true)
RESULT

CHECK(void addSynonym(const String &synonym))
	ptr->addSynonym("blubb_syn3");
	TEST_EQUAL(ptr->getSynonyms().size(), 3)
RESULT

CHECK(const std::set<String>& getSynonyms() const)
	NOT_TESTABLE
RESULT

CHECK(bool operator==(const ResidueModification &modification) const)
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

	mod1.setDiffFormula("C 0 H -2 N 0 O 0");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.setDiffFormula("C 0 H -2 N 0 O 0");
	TEST_EQUAL(mod1 == mod2, true)

	mod1.addSynonym("new_syn");
	TEST_EQUAL(mod1 == mod2, false)
	mod2.addSynonym("new_syn");
	TEST_EQUAL(mod1 == mod2, true)
RESULT

CHECK(bool operator!=(const ResidueModification &modification) const)
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

  mod1.setDiffFormula("C 0 H -2 N 0 O 0");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.setDiffFormula("C 0 H -2 N 0 O 0");
  TEST_EQUAL(mod1 != mod2, false)

  mod1.addSynonym("new_syn");
  TEST_EQUAL(mod1 != mod2, true)
  mod2.addSynonym("new_syn");
  TEST_EQUAL(mod1 != mod2, false)
RESULT


END_TEST
