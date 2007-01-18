// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

Residue* e_ptr = 0;
CHECK(Residue())
	e_ptr = new Residue();
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~Residue())
	delete e_ptr;
RESULT


///////// Modification tests /////////////////

Residue::Modification* m_ptr = 0;
CHECK(Modification())
  m_ptr = new Residue::Modification();
	TEST_NOT_EQUAL(m_ptr, 0)
RESULT

CHECK(~Modification())
	delete m_ptr;
RESULT

m_ptr = new Residue::Modification();

CHECK(Modification(const Modification& modification))
  Residue::Modification m(*m_ptr);
	TEST_EQUAL(m == *m_ptr, true)
RESULT

CHECK(Modification& operator = (const Modification& modification))
	Residue::Modification m;
	m = *m_ptr;
	TEST_EQUAL(m == *m_ptr, true)
RESULT

CHECK(void setName(const String& name))
	m_ptr->setName("the_name");
RESULT

CHECK(const String& getName() const)
	TEST_EQUAL("the_name", m_ptr->getName())
RESULT

CHECK(void setShortName(const String& name))
	m_ptr->setShortName("short_name");
RESULT

CHECK(const String& getShortName() const)
	TEST_EQUAL("short_name", m_ptr->getShortName())
RESULT

CHECK(void setNamePrefix(const String& name_prefix))
	m_ptr->setNamePrefix("the_name_prefix");
RESULT

CHECK(const String& getNamePrefix() const)
	TEST_EQUAL("the_name_prefix", m_ptr->getNamePrefix())
RESULT

set<String> mod_synonyms;
mod_synonyms.insert("syn1");
mod_synonyms.insert("syn2");
CHECK(void setSynonyms(const std::set<String>& synonyms))
	m_ptr->setSynonyms(mod_synonyms);
RESULT

CHECK(void addSynonym(const String& synonym))
	m_ptr->addSynonym("syn3");
RESULT

mod_synonyms.insert("syn3");

CHECK(const std::set<String>& getSynonyms() const)
	TEST_EQUAL(mod_synonyms == m_ptr->getSynonyms(), true)
RESULT

EmpiricalFormula mod_add_formula("CH3");
CHECK(void setAddFormula(const EmpiricalFormula& formula))
	m_ptr->setAddFormula(mod_add_formula);
RESULT

CHECK(const EmpiricalFormula& getAddFormula() const)
	TEST_EQUAL(mod_add_formula == m_ptr->getAddFormula(), true)
RESULT

CHECK(void setAddAverageWeight(Real weight))
	m_ptr->setAddAverageWeight(0.12345);
RESULT

CHECK(Real getAddAverageWeight() const)
	TEST_REAL_EQUAL(m_ptr->getAddAverageWeight(), 0.12345)
RESULT

CHECK(void setAddMonoWeight(Real weight))
	m_ptr->setAddMonoWeight(0.54321);
RESULT

CHECK(Real getAddMonoWeight() const)
	TEST_REAL_EQUAL(m_ptr->getAddMonoWeight(), 0.54321)
RESULT

EmpiricalFormula mod_del_formula("CH4");
CHECK(void setDelFormula(const EmpiricalFormula& formula))
	m_ptr->setDelFormula(mod_del_formula);
RESULT

CHECK(const EmpiricalFormula& getDelFormula() const)
	TEST_EQUAL(m_ptr->getDelFormula() == mod_del_formula, true)
RESULT

CHECK(void setDelAverageWeight(Real weight))
	m_ptr->setDelAverageWeight(0.1234);
RESULT

CHECK(Real getDelAverageWeight() const)
	TEST_REAL_EQUAL(m_ptr->getDelAverageWeight(), 0.1234)
RESULT

CHECK(void setDelMonoWeight(Real weight))
	m_ptr->setDelMonoWeight(0.4321);
RESULT

CHECK(Real getDelMonoWeight() const)
	TEST_REAL_EQUAL(m_ptr->getDelMonoWeight(), 0.4321)
RESULT

set<Residue*> mod_valid_residues;
Residue * mod_valid_res1 = new Residue();
Residue * mod_valid_res2 = new Residue();
mod_valid_residues.insert(mod_valid_res1);

CHECK(void setValidResidues(const std::set<Residue*>& valid_residues))
	m_ptr->setValidResidues(mod_valid_residues);
RESULT

CHECK(void addValidResidue(Residue* valid_residue))
	m_ptr->addValidResidue(mod_valid_res2);
RESULT

mod_valid_residues.insert(mod_valid_res2);

CHECK(const std::set<Residue*>& getValidResidues() const)
	TEST_EQUAL(m_ptr->getValidResidues() == mod_valid_residues, true)
RESULT

CHECK(bool operator == (const Modification& modification) const)
	Residue::Modification m(*m_ptr);
	TEST_EQUAL(m == *m_ptr, true)
	m.addSynonym("syn4");
	TEST_EQUAL(m == *m_ptr, false)
RESULT

CHECK(bool operator != (const Modification& modification) const)
	Residue::Modification m(*m_ptr);
	TEST_EQUAL(m != *m_ptr, false)
	m.addValidResidue(new Residue);
	TEST_EQUAL(m != *m_ptr, true)
RESULT

END_TEST
