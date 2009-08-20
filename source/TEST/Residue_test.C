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

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

Residue* e_ptr = 0;
START_SECTION((Residue()))
	e_ptr = new Residue();
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((virtual ~Residue()))
	delete e_ptr;
END_SECTION

ResidueDB* db = ResidueDB::getInstance();
e_ptr = new Residue(*db->getResidue("LYS"));

EmpiricalFormula h2o("H2O");

START_SECTION((static const EmpiricalFormula& getInternalToFull()))
	TEST_EQUAL(e_ptr->getInternalToFull(), h2o)
END_SECTION


START_SECTION((static DoubleReal getInternalToFullAverageWeight()))
	TEST_EQUAL(e_ptr->getInternalToFullAverageWeight(), h2o.getAverageWeight())
END_SECTION

START_SECTION((static DoubleReal getInternalToFullMonoWeight()))
	TOLERANCE_ABSOLUTE(0.001)
	TEST_REAL_SIMILAR(e_ptr->getInternalToFullMonoWeight(), 18.0106)
END_SECTION

START_SECTION((static const EmpiricalFormula& getNTerminalToFull()))
	TEST_EQUAL(e_ptr->getNTerminalToFull() == EmpiricalFormula("OH"), true)
END_SECTION

START_SECTION((static DoubleReal getNTerminalToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getNTerminalToFullAverageWeight(), 17.0073)
END_SECTION

START_SECTION((static DoubleReal getNTerminalToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getNTerminalToFullMonoWeight(), 17.0027)
END_SECTION

START_SECTION((static const EmpiricalFormula& getCTerminalToFull()))
	TEST_EQUAL(e_ptr->getCTerminalToFull() == EmpiricalFormula("H"), true)
END_SECTION

START_SECTION((static DoubleReal getCTerminalToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getCTerminalToFullAverageWeight(), 1.00794)
END_SECTION

START_SECTION((static DoubleReal getCTerminalToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getCTerminalToFullMonoWeight(), 1.00783)
END_SECTION

START_SECTION((static const EmpiricalFormula& getBIonToFull()))
	TEST_EQUAL(e_ptr->getBIonToFull(), EmpiricalFormula("OH"))
END_SECTION

START_SECTION((static DoubleReal getBIonToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getBIonToFullAverageWeight(), 17.0073)
END_SECTION

START_SECTION((static DoubleReal getBIonToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getBIonToFullMonoWeight(), 17.0027)
END_SECTION

START_SECTION((static const EmpiricalFormula& getAIonToFull()))
	TEST_EQUAL(e_ptr->getAIonToFull(), EmpiricalFormula("O2HC"))
END_SECTION

START_SECTION((static DoubleReal getAIonToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getAIonToFullAverageWeight(), 45.0174)
END_SECTION

START_SECTION((static DoubleReal getAIonToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getAIonToFullMonoWeight(), 44.9976)
END_SECTION

START_SECTION((static const EmpiricalFormula& getYIonToFull()))
	TEST_EQUAL(e_ptr->getYIonToFull(), EmpiricalFormula(""))
END_SECTION

START_SECTION((static DoubleReal getYIonToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getYIonToFullAverageWeight(), 0)
END_SECTION

START_SECTION((static DoubleReal getYIonToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getYIonToFullMonoWeight(), 0)
END_SECTION

START_SECTION((static const EmpiricalFormula& getCIonToFull()))
	TEST_EQUAL(e_ptr->getCIonToFull(), EmpiricalFormula(""))
END_SECTION

START_SECTION((static DoubleReal getCIonToFullAverageWeight()))
	TEST_EQUAL(e_ptr->getCIonToFullAverageWeight(), 0)
END_SECTION

START_SECTION((static DoubleReal getCIonToFullMonoWeight()))
	TEST_EQUAL(e_ptr->getCIonToFullMonoWeight(), 0)
END_SECTION

START_SECTION((static const EmpiricalFormula& getXIonToFull()))
	TEST_EQUAL(e_ptr->getXIonToFull(), EmpiricalFormula("OCH"))
END_SECTION

START_SECTION((static DoubleReal getXIonToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getXIonToFullAverageWeight(), 29.018)
END_SECTION

START_SECTION((static DoubleReal getXIonToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getXIonToFullMonoWeight(), 29.0027)
END_SECTION

START_SECTION((static const EmpiricalFormula& getZIonToFull()))
	TEST_EQUAL(e_ptr->getZIonToFull(), EmpiricalFormula("NH2"))
END_SECTION

START_SECTION((static DoubleReal getZIonToFullAverageWeight()))
	TEST_REAL_SIMILAR(e_ptr->getZIonToFullAverageWeight(), 16.0226)
END_SECTION

START_SECTION((static DoubleReal getZIonToFullMonoWeight()))
	TEST_REAL_SIMILAR(e_ptr->getZIonToFullMonoWeight(), 16.0187)
END_SECTION

START_SECTION(Residue(const Residue &residue))
	Residue copy(*e_ptr);
	TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(Residue(const String &name, const String &three_letter_code, const String &one_letter_code, const EmpiricalFormula &formula))
	Residue copy(e_ptr->getName(), e_ptr->getThreeLetterCode(), e_ptr->getOneLetterCode(), e_ptr->getFormula());
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getThreeLetterCode(), e_ptr->getThreeLetterCode())
	TEST_EQUAL(copy.getOneLetterCode(), e_ptr->getOneLetterCode())
	TEST_EQUAL(copy.getFormula(), e_ptr->getFormula())
END_SECTION

START_SECTION(Residue& operator=(const Residue &residue))
	Residue copy;
	copy = *e_ptr;
	TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(void setName(const String &name))
	Residue copy(*e_ptr);
	e_ptr->setName("BLUBB");
	TEST_NOT_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), "BLUBB")
END_SECTION

START_SECTION(void setShortName(const String &short_name))
	Residue copy(*e_ptr);
	e_ptr->setShortName("BB");
	TEST_NOT_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(const String& getShortName() const)
	TEST_EQUAL(e_ptr->getShortName(), "BB")
END_SECTION

START_SECTION(void setSynonyms(const std::set< String > &synonyms))
	Residue copy(*e_ptr);
	set<String> syn;
	syn.insert("BLI");
	syn.insert("BLA");
	e_ptr->setSynonyms(syn);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION 

START_SECTION(void addSynonym(const String &synonym))
	Residue copy(*e_ptr);
	e_ptr->addSynonym("BLUFF");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
	TEST_EQUAL(e_ptr->getSynonyms().size(), 3)
END_SECTION

START_SECTION(void setThreeLetterCode(const String &three_letter_code))
	Residue copy(*e_ptr);
	e_ptr->setThreeLetterCode("BLA");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getThreeLetterCode() const)
	TEST_EQUAL(e_ptr->getThreeLetterCode(), "BLA")
END_SECTION

START_SECTION(void setOneLetterCode(const String &one_letter_code))
	Residue copy(*e_ptr);
	e_ptr->setOneLetterCode("B");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getOneLetterCode() const)
	TEST_EQUAL(e_ptr->getOneLetterCode(), "B")
END_SECTION

START_SECTION(void addLossFormula(const EmpiricalFormula&))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy)
	e_ptr->addLossFormula(EmpiricalFormula("H2O"));
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(void setLossFormulas(const std::vector<EmpiricalFormula> &))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy)
	vector<EmpiricalFormula> losses;
	losses.push_back(EmpiricalFormula("H2O"));
	e_ptr->setLossFormulas(losses);
	TEST_NOT_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getLossFormulas() const)
	vector<EmpiricalFormula> losses;
	losses.push_back(EmpiricalFormula("H2O"));
	TEST_EQUAL(e_ptr->getLossFormulas() == losses, true)
END_SECTION

START_SECTION(void setLossNames(const std::vector<String> &name))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy)
	vector<String> names;
	names.push_back("Waesserchen");
	e_ptr->setLossNames(names);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::vector<String>& getLossNames() const)
	vector<String> names;
	names.push_back("Waesserchen");
	TEST_EQUAL(e_ptr->getLossNames() == names, true)
END_SECTION

START_SECTION(void addLossName(const String& name))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy)
	copy.addLossName("Waesserchen2");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(void setNTermLossFormulas(const std::vector< EmpiricalFormula > &))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy)
	vector<EmpiricalFormula> losses;
	losses.push_back(EmpiricalFormula("H3O"));
	e_ptr->setNTermLossFormulas(losses);
	TEST_NOT_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getNTermLossFormulas() const)
	vector<EmpiricalFormula> losses;
	losses.push_back(EmpiricalFormula("H3O"));
	TEST_EQUAL(e_ptr->getNTermLossFormulas() == losses, true)
END_SECTION

START_SECTION(void addNTermLossFormula(const EmpiricalFormula&))
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  e_ptr->addNTermLossFormula(EmpiricalFormula("H4O"));
  TEST_NOT_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION(void setNTermLossNames(const std::vector< String > &name))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy);
	vector<String> names;
	names.push_back("Nwaesserchen");
	e_ptr->setNTermLossNames(names);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::vector<String>& getNTermLossNames() const)
	vector<String> names;
	names.push_back("Nwaesserchen");
	TEST_EQUAL(e_ptr->getNTermLossNames() == names, true)
END_SECTION

START_SECTION(void addNTermLossName(const String &name))
	Residue copy(*e_ptr);
	TEST_EQUAL(*e_ptr, copy);
	e_ptr->addNTermLossName("Nwaesserchen2");
	TEST_NOT_EQUAL(*e_ptr == copy, true);
END_SECTION

START_SECTION(bool hasNTermNeutralLosses() const)
	Residue copy(*e_ptr);
	TEST_EQUAL(copy.hasNTermNeutralLosses(), true)
	copy.setNTermLossFormulas(vector<EmpiricalFormula>());
	copy.setNTermLossNames(vector<String>());
	TEST_EQUAL(copy.hasNTermNeutralLosses(), false)
END_SECTION

START_SECTION(void setFormula(const EmpiricalFormula &formula, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setFormula(EmpiricalFormula("C2H6O"));
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(EmpiricalFormula getFormula(ResidueType res_type=Full) const)
	TEST_EQUAL(e_ptr->getFormula(), EmpiricalFormula("C2H6O"))
END_SECTION

START_SECTION(void setAverageWeight(DoubleReal weight, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setAverageWeight(123.4);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getAverageWeight(ResidueType res_type=Full) const)
	TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), 123.4)
END_SECTION
    
START_SECTION(void setMonoWeight(DoubleReal weight, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setMonoWeight(1234.5);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getMonoWeight(ResidueType res_type=Full) const)
	TEST_REAL_SIMILAR(e_ptr->getMonoWeight(), 1234.5)
END_SECTION
 

START_SECTION(void setModification(const String& name))
	e_ptr->setOneLetterCode("M"); // we need M for this mod
	e_ptr->setModification("Oxidation");
	TEST_EQUAL(e_ptr->getModification(), "Oxidation")
	e_ptr->setOneLetterCode("B");
END_SECTION


START_SECTION(const String& getModification() const)
	NOT_TESTABLE
END_SECTION
 
START_SECTION(void setLowMassIons(const std::vector< EmpiricalFormula > &low_mass_ions))
	Residue copy(*e_ptr);
	vector<EmpiricalFormula> ions;
	ions.push_back(EmpiricalFormula("NH3"));
	ions.push_back(EmpiricalFormula("PO4"));
	e_ptr->setLowMassIons(ions);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getLowMassIons() const)
	TEST_EQUAL(e_ptr->getLowMassIons()[0], EmpiricalFormula("NH3"))
END_SECTION

START_SECTION(bool hasNeutralLoss() const)
	Residue res;
	TEST_EQUAL(res.hasNeutralLoss(), false)
	res.addLossFormula(EmpiricalFormula("H2O"));
	TEST_EQUAL(res.hasNeutralLoss(), true)
END_SECTION

START_SECTION(bool operator==(const Residue &residue) const)
	Residue r;
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setName("other_name");
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setShortName("other_short_name");
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	set<String> syns;
	syns.insert("new_syn");
	r.setSynonyms(syns);
	TEST_EQUAL(r == *e_ptr, false)

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setThreeLetterCode("new_3lc");
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setOneLetterCode("new_1lc");
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.addLossFormula(EmpiricalFormula("C1H3"));
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.addLossName("new_loss_name");
	TEST_EQUAL(r == *e_ptr, false)

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setFormula(EmpiricalFormula("C16H18N3O5"));
	TEST_EQUAL(r == *e_ptr, false)

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setAverageWeight(12345.678);
	TEST_EQUAL(r == *e_ptr, false)

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setMonoWeight(12345.6789);
	TEST_EQUAL(r == *e_ptr, false)
				
	e_ptr->setOneLetterCode("M");
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setModification("Oxidation");
	TEST_EQUAL(r == *e_ptr, false)
	e_ptr->setOneLetterCode("B");

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true);
	vector<EmpiricalFormula> low_mass_ions;
	low_mass_ions.push_back(EmpiricalFormula("H"));
	r.setLowMassIons(low_mass_ions);
	TEST_EQUAL(r == *e_ptr, false);

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setPka(123456.789);
	TEST_EQUAL(r == *e_ptr, false);

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true);
	r.setPkb(1234567.89);
	TEST_EQUAL(r == *e_ptr, false);

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true);
	r.setPkc(12345678.9);
	TEST_EQUAL(r == *e_ptr, false)

	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setSideChainBasicity(111.2345);
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setBackboneBasicityLeft(1112.345);
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setBackboneBasicityRight(11123.45);
	TEST_EQUAL(r == *e_ptr, false)
END_SECTION
    
START_SECTION(bool operator!=(const Residue &residue) const)
  Residue r;
  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setName("other_name");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setShortName("other_short_name");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  set<String> syns;
  syns.insert("new_syn");
  r.setSynonyms(syns);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setThreeLetterCode("new_3lc");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setOneLetterCode("new_1lc");
  TEST_EQUAL(r != *e_ptr, true)

	r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.addLossFormula(EmpiricalFormula("C1H3"));
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.addLossName("new_loss_name");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setFormula(EmpiricalFormula("C16H18N3O5"));
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setAverageWeight(12345.678);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setMonoWeight(12345.6789);
  TEST_EQUAL(r != *e_ptr, true)


	r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false);
  vector<EmpiricalFormula> low_mass_ions;
  low_mass_ions.push_back(EmpiricalFormula("H"));
  r.setLowMassIons(low_mass_ions);
  TEST_EQUAL(r != *e_ptr, true);

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setPka(123456.789);
  TEST_EQUAL(r != *e_ptr, true);

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false);
  r.setPkb(1234567.89);
  TEST_EQUAL(r != *e_ptr, true);

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false);
  r.setPkc(12345678.9);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setSideChainBasicity(111.2345);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setBackboneBasicityLeft(1112.345);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setBackboneBasicityRight(11123.45);
  TEST_EQUAL(r != *e_ptr, true)
				
END_SECTION

START_SECTION(bool operator==(char one_letter_code) const)
	TEST_EQUAL(*e_ptr == 'B', true)
END_SECTION
    
START_SECTION(bool operator!=(char one_letter_code) const)
	TEST_EQUAL(*e_ptr != 'C', true)
END_SECTION

START_SECTION(void setPka(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPka(345.5);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getPka() const)
	TEST_REAL_SIMILAR(e_ptr->getPka(), 345.5)
END_SECTION

START_SECTION(void setPkb(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPkb(675.8);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getPkb() const)
	TEST_REAL_SIMILAR(e_ptr->getPkb(), 675.8)
END_SECTION

START_SECTION(void setPkc(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPkc(9329.0);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getPkc() const)
	TEST_REAL_SIMILAR(e_ptr->getPkc(), 9329.0)
END_SECTION

START_SECTION(DoubleReal getPiValue() const)
	TEST_REAL_SIMILAR(e_ptr->getPiValue(), 4837.25)
END_SECTION

START_SECTION(void setSideChainBasicity(DoubleReal gb_sc))
	Residue copy(*e_ptr);
	e_ptr->setSideChainBasicity(654.3);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getSideChainBasicity() const)
	TEST_REAL_SIMILAR(e_ptr->getSideChainBasicity(), 654.3)
END_SECTION


START_SECTION(void setBackboneBasicityLeft(DoubleReal gb_bb_l))
	Residue copy(*e_ptr);
	e_ptr->setBackboneBasicityLeft(123.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(DoubleReal getBackboneBasicityLeft() const)
	TEST_REAL_SIMILAR(e_ptr->getBackboneBasicityLeft(), 123.6)
END_SECTION


START_SECTION(void setBackboneBasicityRight(DoubleReal gb_bb_r))
	Residue copy(*e_ptr);
	e_ptr->setBackboneBasicityRight(12345.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION


START_SECTION(DoubleReal getBackboneBasicityRight() const)
	TEST_REAL_SIMILAR(e_ptr->getBackboneBasicityRight(), 12345.6)
END_SECTION

START_SECTION(bool isModified() const)
	Residue res;
	res.setOneLetterCode("M"); // we need M for this mod
	TEST_EQUAL(res.isModified(), false)
	res.setModification("Oxidation");
	TEST_EQUAL(res.isModified(), true)
END_SECTION


START_SECTION((void setResidueSets(const std::set< String > &residues_sets)))
	set<String> res_sets;
	res_sets.insert("rs1");
	res_sets.insert("rs2");
	e_ptr->setResidueSets(res_sets);
	TEST_EQUAL(res_sets == e_ptr->getResidueSets(), true)
END_SECTION

START_SECTION((void addResidueSet(const String &residue_sets)))
	e_ptr->addResidueSet("rs3");
	TEST_EQUAL(e_ptr->getResidueSets().size(), 3)
END_SECTION

START_SECTION((const std::set<String>& getResidueSets() const))
	set<String> res_sets;
	res_sets.insert("rs1");
	res_sets.insert("rs2");
	res_sets.insert("rs3");
	TEST_EQUAL(e_ptr->getResidueSets() == res_sets, true)
END_SECTION

START_SECTION((bool isInResidueSet(const String &residue_set)))
	TEST_EQUAL(e_ptr->isInResidueSet("rs1"), true)
	TEST_EQUAL(e_ptr->isInResidueSet("rs3"), true)
	TEST_EQUAL(e_ptr->isInResidueSet("rs4"), false)
END_SECTION

END_TEST

