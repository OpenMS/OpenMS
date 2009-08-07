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

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

ResidueDB* ptr = 0;
START_SECTION(ResidueDB* getInstance())
	ptr = ResidueDB::getInstance();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(virtual ~ResidueDB())
	NOT_TESTABLE
END_SECTION

START_SECTION((const Residue* getResidue(const String &name) const))
  TEST_EQUAL(ptr->getResidue("C")->getOneLetterCode(), "C")
END_SECTION

START_SECTION((bool hasResidue(const String &name) const))
  TEST_EQUAL(ptr->hasResidue("BLUBB"), false)
	TEST_EQUAL(ptr->hasResidue("LYS"), true)
	TEST_EQUAL(ptr->hasResidue("K"), true)
END_SECTION

START_SECTION(bool hasResidue(const Residue *residue) const)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("BLUBB")), false)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("LYS")), true)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("K")), true)
END_SECTION

START_SECTION(Size getNumberOfResidues() const)
	TEST_EQUAL(ptr->getNumberOfResidues() >= 20 , true);
END_SECTION

START_SECTION(const Residue* getModifiedResidue(const String &name))
	const Residue* mod_res = ptr->getModifiedResidue("MOD:00720"); // ox methionine
	TEST_STRING_EQUAL(mod_res->getOneLetterCode(), "M")
	TEST_STRING_EQUAL(mod_res->getModification(), "MOD:00720")
END_SECTION

START_SECTION(const Residue* getModifiedResidue(const Residue *residue, const String &name))
	const Residue* mod_res = ptr->getModifiedResidue(ptr->getResidue("M"), "MOD:00720");
	TEST_STRING_EQUAL(mod_res->getOneLetterCode(), "M")
	TEST_STRING_EQUAL(mod_res->getModification(), "MOD:00720")
END_SECTION

START_SECTION(const std::set<const Residue*> getResidues(AminoAcidSet aa_set = "All") const)
	set<const Residue*> residues = ptr->getResidues("All");
	TEST_EQUAL(residues.size() >= 21, true)
	residues = ptr->getResidues("Natural20");
	TEST_EQUAL(residues.size(), 20)
	residues = ptr->getResidues("Natural19WithoutL");
	TEST_EQUAL(residues.size(), 19)
END_SECTION

START_SECTION(void setResidues(const String &filename))
	NOT_TESTABLE // this method is hard to test, just provided for convenience
END_SECTION
    
START_SECTION(void addResidue(const Residue &residue))
	TEST_EQUAL(ptr->hasResidue("UGU"), false)
	TEST_EQUAL(ptr->hasResidue("$"), false)
	Residue res;
	res.setShortName("$");
	res.setOneLetterCode("$");
	res.setThreeLetterCode("UGU");
	res.setName("MyLittleUGUResidue");
	res.setFormula(EmpiricalFormula("C3H4O4"));
	ptr->addResidue(res);
	TEST_EQUAL(ptr->hasResidue("UGU"), true)
	TEST_EQUAL(ptr->hasResidue("$"), true)
END_SECTION

START_SECTION(ResidueIterator beginResidue())
	ResidueDB::ResidueIterator it = ptr->beginResidue();
	Size count(0);
	while (it != ptr->endResidue())
	{
		++it;
		++count;
	}

	TEST_EQUAL(count >= 22, true)
END_SECTION
  
START_SECTION(ResidueIterator endResidue())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(ResidueConstIterator beginResidue() const)
	const ResidueDB* const_ptr = ptr;
	ResidueDB::ResidueConstIterator it = const_ptr->beginResidue();
	Size count(0);
	while (it != const_ptr->endResidue())
	{
		++it;
		++count;
	}
	TEST_EQUAL(count >= 22, true)
END_SECTION

START_SECTION(ResidueConstIterator endResidue() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(Size getNumberOfModifiedResidues() const)
	TEST_EQUAL(ptr->getNumberOfModifiedResidues(), 1)
	const Residue* mod_res = 0;
	mod_res = ptr->getModifiedResidue("MOD:01214");
	TEST_NOT_EQUAL(mod_res, 0)
	TEST_EQUAL(ptr->getNumberOfModifiedResidues(), 2)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
