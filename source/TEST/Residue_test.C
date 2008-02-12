// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

Residue* e_ptr = 0;
CHECK((Residue()))
	e_ptr = new Residue();
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((virtual ~Residue()))
	delete e_ptr;
RESULT

ResidueDB db;
e_ptr = new Residue(*db.getResidue("LYS"));

EmpiricalFormula h2o("H2O");

CHECK((static const EmpiricalFormula& getInternalToFull()))
	TEST_EQUAL(e_ptr->getInternalToFull(), h2o)
RESULT


CHECK((static DoubleReal getInternalToFullAverageWeight()))
	TEST_EQUAL(e_ptr->getInternalToFullAverageWeight(), h2o.getAverageWeight())
RESULT

CHECK((static DoubleReal getInternalToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getNTerminalToFull()))

RESULT

CHECK((static DoubleReal getNTerminalToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getNTerminalToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getCTerminalToFull()))

RESULT

CHECK((static DoubleReal getCTerminalToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getCTerminalToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getBIonToFull()))

RESULT

CHECK((static DoubleReal getBIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getBIonToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getAIonToFull()))

RESULT

CHECK((static DoubleReal getAIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getAIonToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getYIonToFull()))

RESULT

CHECK((static DoubleReal getYIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getYIonToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getCIonToFull()))

RESULT

CHECK((static DoubleReal getCIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getCIonToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getXIonToFull()))

RESULT

CHECK((static DoubleReal getXIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getXIonToFullMonoWeight()))

RESULT

CHECK((static const EmpiricalFormula& getZIonToFull()))

RESULT

CHECK((static DoubleReal getZIonToFullAverageWeight()))

RESULT

CHECK((static DoubleReal getZIonToFullMonoWeight()))

RESULT

CHECK(Residue(const Residue &residue))
	Residue copy(*e_ptr);
	TEST_EQUAL(copy, *e_ptr)
RESULT

CHECK(Residue(const String &name, const String &three_letter_code, const String &one_letter_code, const EmpiricalFormula &formula, const EmpiricalFormula &neutral_loss))
	Residue copy(e_ptr->getName(), e_ptr->getThreeLetterCode(), e_ptr->getOneLetterCode(), e_ptr->getFormula(), e_ptr->getLossFormula());
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getThreeLetterCode(), e_ptr->getThreeLetterCode())
	TEST_EQUAL(copy.getOneLetterCode(), e_ptr->getOneLetterCode())
	TEST_EQUAL(copy.getFormula(), e_ptr->getFormula())
	TEST_EQUAL(copy.getLossFormula(), e_ptr->getLossFormula())
RESULT

CHECK(Residue& operator=(const Residue &residue))
	Residue copy;
	copy = *e_ptr;
	TEST_EQUAL(copy, *e_ptr)
RESULT

CHECK(void setName(const String &name))
	Residue copy(*e_ptr);
	e_ptr->setName("BLUBB");
	TEST_NOT_EQUAL(copy, *e_ptr)
RESULT

CHECK(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), "BLUBB")
RESULT

CHECK(void setShortName(const String &short_name))
	Residue copy(*e_ptr);
	e_ptr->setShortName("BB");
	TEST_NOT_EQUAL(copy, *e_ptr)
RESULT

CHECK(const String& getShortName() const)
	TEST_EQUAL(e_ptr->getShortName(), "BB")
RESULT

CHECK(void setSynonyms(const std::set< String > &synonyms))
	Residue copy(*e_ptr);
	set<String> syn;
	syn.insert("BLI");
	syn.insert("BLA");
	e_ptr->setSynonyms(syn);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT 

CHECK(void addSynonym(const String &synonym))
	Residue copy(*e_ptr);
	e_ptr->addSynonym("BLUFF");
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const std::set<String>& getSynonyms() const)
	TEST_EQUAL(e_ptr->getSynonyms().size(), 3)
RESULT

CHECK(void setThreeLetterCode(const String &three_letter_code))
	Residue copy(*e_ptr);
	e_ptr->setThreeLetterCode("BLA");
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const String& getThreeLetterCode() const)
	TEST_EQUAL(e_ptr->getThreeLetterCode(), "BLA")
RESULT

CHECK(void setOneLetterCode(const String &one_letter_code))
	Residue copy(*e_ptr);
	e_ptr->setOneLetterCode("B");
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const String& getOneLetterCode() const)
	TEST_EQUAL(e_ptr->getOneLetterCode(), "B")
RESULT

CHECK(void setLossFormula(const EmpiricalFormula &))
	Residue copy(*e_ptr);
	e_ptr->setLossFormula(EmpiricalFormula("H2O"));
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const EmpiricalFormula& getLossFormula() const)
	TEST_EQUAL(e_ptr->getLossFormula(), EmpiricalFormula("H2O"))
RESULT

CHECK(void setLossAverageWeight(DoubleReal weight))
	Residue copy(*e_ptr);
	e_ptr->setLossAverageWeight(18.5);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getLossAverageWeight() const)
	TEST_REAL_EQUAL(e_ptr->getLossAverageWeight(), 18.5)
RESULT

CHECK(void setLossMonoWeight(DoubleReal weight))
	Residue copy(*e_ptr);
	e_ptr->setLossMonoWeight(18.6);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getLossMonoWeight() const)
	TEST_EQUAL(e_ptr->getLossMonoWeight(), 18.6)
RESULT

CHECK(void setLossName(const String &name))
	Residue copy(*e_ptr);
	e_ptr->setLossName("Waesserchen");
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const String& getLossName() const)
	TEST_EQUAL(e_ptr->getLossName(), "Waesserchen")
RESULT

CHECK(void setFormula(const EmpiricalFormula &formula, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setFormula(EmpiricalFormula("C2H6O"));
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(EmpiricalFormula getFormula(ResidueType res_type=Full) const)
	TEST_EQUAL(e_ptr->getFormula(), EmpiricalFormula("C2H6O"))
RESULT

CHECK(void setAverageWeight(DoubleReal weight, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setAverageWeight(123.4);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getAverageWeight(ResidueType res_type=Full) const)
	TEST_REAL_EQUAL(e_ptr->getAverageWeight(), 123.4)
RESULT
    
CHECK(void setMonoWeight(DoubleReal weight, ResidueType res_type=Full))
	Residue copy(*e_ptr);
	e_ptr->setMonoWeight(1234.5);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getMonoWeight(ResidueType res_type=Full) const)
	TEST_REAL_EQUAL(e_ptr->getMonoWeight(), 1234.5)
RESULT
 
CHECK(void setModification(ResidueModification *modification))
	Residue copy(*e_ptr);
	ResidueModification* mod = new ResidueModification();
	mod->setName("DA_MOD");
	e_ptr->setModification(mod);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const ResidueModification* getModification() const)
	TEST_EQUAL(e_ptr->getModification()->getName(), "DA_MOD")
RESULT
 
CHECK(void setUnmodifiedName(const String &name))
	Residue copy(*e_ptr);
	e_ptr->setUnmodifiedName("NATURAL");
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT
    
CHECK(const String& getUnmodifiedName() const)
	TEST_EQUAL(e_ptr->getUnmodifiedName(), "NATURAL")
RESULT
    
CHECK(void setLowMassIons(const std::vector< EmpiricalFormula > &low_mass_ions))
	Residue copy(*e_ptr);
	vector<EmpiricalFormula> ions;
	ions.push_back(EmpiricalFormula("NH3"));
	ions.push_back(EmpiricalFormula("PO4"));
	e_ptr->setLowMassIons(ions);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(const std::vector<EmpiricalFormula>& getLowMassIons() const)
	TEST_EQUAL(e_ptr->getLowMassIons()[0], EmpiricalFormula("NH3"))
RESULT

CHECK(bool hasNeutralLoss() const)
	Residue res;
	TEST_EQUAL(res.hasNeutralLoss(), false)
	res.setLossFormula(EmpiricalFormula("H2O"));
	TEST_EQUAL(res.hasNeutralLoss(), true)
RESULT

CHECK(bool operator==(const Residue &residue) const)
	// implicitly testet above
RESULT
    
CHECK(bool operator!=(const Residue &residue) const)
	// implicitly tested above
RESULT

CHECK(bool operator==(char one_letter_code) const)
	TEST_EQUAL(*e_ptr == 'B', true)
RESULT
    
CHECK(bool operator!=(char one_letter_code) const)
	TEST_EQUAL(*e_ptr != 'C', true)
RESULT

CHECK(void setPka(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPka(345.5);
	TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getPka() const)
	TEST_REAL_EQUAL(e_ptr->getPka(), 345.5)
RESULT

CHECK(void setPkb(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPkb(675.8);
  TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getPkb() const)
	TEST_REAL_EQUAL(e_ptr->getPkb(), 675.8)
RESULT

CHECK(void setPkc(DoubleReal value))
	Residue copy(*e_ptr);
	e_ptr->setPkc(9329.0);
  TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getPkc() const)
	TEST_REAL_EQUAL(e_ptr->getPkc(), 9329.0)
RESULT

CHECK(DoubleReal getPiValue() const)
	TEST_REAL_EQUAL(e_ptr->getPiValue(), 4837.25)
RESULT

CHECK(void setSideChainBasicity(DoubleReal gb_sc))
	Residue copy(*e_ptr);
	e_ptr->setSideChainBasicity(654.3);
  TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getSideChainBasicity() const)
	TEST_REAL_EQUAL(e_ptr->getSideChainBasicity(), 654.3)
RESULT


CHECK(void setBackboneBasicityLeft(DoubleReal gb_bb_l))
	Residue copy(*e_ptr);
	e_ptr->setBackboneBasicityLeft(123.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
RESULT

CHECK(DoubleReal getBackboneBasicityLeft() const)
	TEST_REAL_EQUAL(e_ptr->getBackboneBasicityLeft(), 123.6)
RESULT


CHECK(void setBackboneBasicityRight(DoubleReal gb_bb_r))
	Residue copy(*e_ptr);
	e_ptr->setBackboneBasicityRight(12345.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
RESULT


CHECK(DoubleReal getBackboneBasicityRight() const)
	TEST_REAL_EQUAL(e_ptr->getBackboneBasicityRight(), 12345.6)
RESULT


CHECK(bool isModified() const)
	ResidueModification mod;
	Residue res;
	TEST_EQUAL(res.isModified(), false)
	res.setModification(&mod);
	TEST_EQUAL(res.isModified(), true)
RESULT


END_TEST

