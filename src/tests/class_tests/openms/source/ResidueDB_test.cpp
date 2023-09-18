// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

ResidueDB* ptr = nullptr;
ResidueDB* nullPointer = nullptr;
START_SECTION(ResidueDB* getInstance())
	ptr = ResidueDB::getInstance();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~ResidueDB())
	NOT_TESTABLE
END_SECTION

START_SECTION((const Residue* getResidue(const String& name) const))
  TEST_EQUAL(ptr->getResidue("C")->getOneLetterCode(), "C")
END_SECTION

START_SECTION((const Residue* getResidue(const unsigned char& one_letter_code) const))
  TEST_EQUAL(ptr->getResidue('C')->getOneLetterCode(), "C")
END_SECTION

START_SECTION((bool hasResidue(const String& name) const))
  TEST_EQUAL(ptr->hasResidue("BLUBB"), false)
	TEST_EQUAL(ptr->hasResidue("Lys"), true)
	TEST_EQUAL(ptr->hasResidue("K"), true)
END_SECTION

START_SECTION(bool hasResidue(const Residue* residue) const)
	TEST_EXCEPTION(Exception::InvalidValue, ptr->hasResidue(ptr->getResidue("BLUBB")))
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("Lys")), true)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("K")), true)
END_SECTION

START_SECTION(Size getNumberOfResidues() const)
	TEST_EQUAL(ptr->getNumberOfResidues() >= 20 , true);
END_SECTION

START_SECTION(const Residue* getModifiedResidue(const String& name))
	const Residue* mod_res = ptr->getModifiedResidue("Oxidation (M)"); // ox methionine
	TEST_STRING_EQUAL(mod_res->getOneLetterCode(), "M")
	TEST_STRING_EQUAL(mod_res->getModificationName(), "Oxidation")
END_SECTION

START_SECTION(const Residue* getModifiedResidue(const Residue* residue, const String& name))
	const Residue* mod_res = ptr->getModifiedResidue(ptr->getResidue("M"), "Oxidation (M)");
	TEST_STRING_EQUAL(mod_res->getOneLetterCode(), "M")
	TEST_STRING_EQUAL(mod_res->getModificationName(), "Oxidation")

	const Residue* term_mod_res = ptr->getModifiedResidue(ptr->getResidue("C"), "Pyro-carbamidomethyl (N-term C)");
	TEST_STRING_EQUAL(term_mod_res->getOneLetterCode(), "C")
	TEST_STRING_EQUAL(term_mod_res->getModificationName(), "Pyro-carbamidomethyl")
END_SECTION

START_SECTION((const std::set<const Residue*> getResidues(const String& residue_set="All") const))
	set<const Residue*> residues = ptr->getResidues("All");
	TEST_EQUAL(residues.size() >= 21, true)
	residues = ptr->getResidues("Natural20");
	TEST_EQUAL(residues.size(), 20)
	residues = ptr->getResidues("Natural19WithoutL");
	TEST_EQUAL(residues.size(), 19)
END_SECTION

START_SECTION((const std::set<String>& getResidueSets() const))
	set<String> res_sets = ResidueDB::getInstance()->getResidueSets();
	TEST_EQUAL(res_sets.find("All") != res_sets.end(), true)
	TEST_EQUAL(res_sets.find("Natural20") != res_sets.end(), true)
	TEST_EQUAL(res_sets.find("Natural19WithoutL") != res_sets.end(), true)
	TEST_EQUAL(res_sets.find("Natural19WithoutI") != res_sets.end(), true)
END_SECTION


START_SECTION(void setResidues(const String& filename))
	NOT_TESTABLE // this method is hard to test, just provided for convenience
END_SECTION

START_SECTION(Size getNumberOfModifiedResidues() const)
	TEST_EQUAL(ptr->getNumberOfModifiedResidues(), 1)
	const Residue* mod_res = nullptr;
    const Residue* mod_res_nullPointer = nullptr;
	mod_res = ptr->getModifiedResidue("Carbamidomethyl (C)");
    TEST_NOT_EQUAL(mod_res, mod_res_nullPointer)
	TEST_EQUAL(ptr->getNumberOfModifiedResidues(), 2)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
