// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

Residue* e_ptr = nullptr;
Residue* e_nullPointer = nullptr;
START_SECTION((Residue()))
{
  e_ptr = new Residue();
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
}
END_SECTION

START_SECTION((virtual ~Residue()))
{
  delete e_ptr;
}
END_SECTION

ResidueDB* db = ResidueDB::getInstance();
e_ptr = new Residue(*db->getResidue("Lys"));

EmpiricalFormula h2o("H2O");

START_SECTION((static const EmpiricalFormula& getInternalToFull()))
{
  TEST_EQUAL(e_ptr->getInternalToFull(), h2o)
}
END_SECTION

TOLERANCE_ABSOLUTE(0.001)

START_SECTION((static const EmpiricalFormula& getInternalToNTerm()))
  TEST_EQUAL(e_ptr->getInternalToNTerm(), EmpiricalFormula("H"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToCTerm()))
  TEST_EQUAL(e_ptr->getInternalToCTerm(), EmpiricalFormula("OH"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToAIon()))
  TEST_EQUAL(e_ptr->getInternalToAIon(), EmpiricalFormula("")-EmpiricalFormula("CO"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToBIon()))
  TEST_EQUAL(e_ptr->getInternalToBIon(), EmpiricalFormula(""))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToCIon()))
  TEST_EQUAL(e_ptr->getInternalToCIon(), EmpiricalFormula("NH3"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToXIon()))
  TEST_EQUAL(e_ptr->getInternalToXIon(), EmpiricalFormula("CO2"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToYIon()))
  TEST_EQUAL(e_ptr->getInternalToYIon(), EmpiricalFormula("H2O"))
END_SECTION

START_SECTION((static const EmpiricalFormula& getInternalToZIon()))
  TEST_EQUAL(e_ptr->getInternalToZIon(), EmpiricalFormula("OH") - EmpiricalFormula("NH2"))
END_SECTION

START_SECTION(Residue(const Residue &residue))
  Residue copy(*e_ptr);
  TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(Residue(const String &name, const String &three_letter_code, const String &one_letter_code, const EmpiricalFormula &formula))
{
  Residue copy(e_ptr->getName(), e_ptr->getThreeLetterCode(), e_ptr->getOneLetterCode(), e_ptr->getFormula());
  TEST_EQUAL(copy.getName(), e_ptr->getName())
  TEST_EQUAL(copy.getThreeLetterCode(), e_ptr->getThreeLetterCode())
  TEST_EQUAL(copy.getOneLetterCode(), e_ptr->getOneLetterCode())
  TEST_EQUAL(copy.getFormula(), e_ptr->getFormula())
}
END_SECTION

START_SECTION(Residue& operator=(const Residue &residue))
{
  Residue copy;
  copy = *e_ptr;
  TEST_EQUAL(copy, *e_ptr)
}
END_SECTION

START_SECTION(void setName(const String &name))
{
  Residue copy(*e_ptr);
  e_ptr->setName("BLUBB");
  TEST_NOT_EQUAL(copy, *e_ptr)
}
END_SECTION

START_SECTION(const String& getName() const)
{
  TEST_EQUAL(e_ptr->getName(), "BLUBB")
}
END_SECTION

START_SECTION(void setSynonyms(const std::set< String > &synonyms))
{
  Residue copy(*e_ptr);
  set<String> syn;
  syn.insert("BLI");
  syn.insert("BLA");
  e_ptr->setSynonyms(syn);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(void addSynonym(const String &synonym))
{
  Residue copy(*e_ptr);
  e_ptr->addSynonym("BLUFF");
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
{
  TEST_EQUAL(e_ptr->getSynonyms().size(), 3)
}
END_SECTION

START_SECTION(void setThreeLetterCode(const String &three_letter_code))
{
  Residue copy(*e_ptr);
  e_ptr->setThreeLetterCode("BLA");
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const String& getThreeLetterCode() const)
{
  TEST_EQUAL(e_ptr->getThreeLetterCode(), "BLA")
}
END_SECTION

START_SECTION(void setOneLetterCode(const String &one_letter_code))
{
  Residue copy(*e_ptr);
  e_ptr->setOneLetterCode("B");
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const String& getOneLetterCode() const)
{
  TEST_EQUAL(e_ptr->getOneLetterCode(), "B")
}
END_SECTION

START_SECTION(void addLossFormula(const EmpiricalFormula&))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  e_ptr->addLossFormula(EmpiricalFormula("H2O"));
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(void setLossFormulas(const std::vector<EmpiricalFormula> &))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  std::vector<EmpiricalFormula> losses;
  losses.push_back(EmpiricalFormula("H2O"));
  e_ptr->setLossFormulas(losses);
  TEST_NOT_EQUAL(*e_ptr == copy, true)
}
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getLossFormulas() const)
{
  std::vector<EmpiricalFormula> losses;
  losses.push_back(EmpiricalFormula("H2O"));
  TEST_EQUAL(e_ptr->getLossFormulas() == losses, true)
}
END_SECTION

START_SECTION(void setLossNames(const std::vector<String> &name))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  vector<String> names;
  names.push_back("Waesserchen");
  e_ptr->setLossNames(names);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const std::vector<String>& getLossNames() const)
{
  vector<String> names;
  names.push_back("Waesserchen");
  TEST_EQUAL(e_ptr->getLossNames() == names, true)
}
END_SECTION

START_SECTION(void addLossName(const String& name))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  copy.addLossName("Waesserchen2");
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(void setNTermLossFormulas(const std::vector< EmpiricalFormula > &))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  vector<EmpiricalFormula> losses;
  losses.push_back(EmpiricalFormula("H3O"));
  e_ptr->setNTermLossFormulas(losses);
  TEST_NOT_EQUAL(*e_ptr == copy, true)
}
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getNTermLossFormulas() const)
{
  vector<EmpiricalFormula> losses;
  losses.push_back(EmpiricalFormula("H3O"));
  TEST_EQUAL(e_ptr->getNTermLossFormulas() == losses, true)
}
END_SECTION

START_SECTION(void addNTermLossFormula(const EmpiricalFormula&))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy)
  e_ptr->addNTermLossFormula(EmpiricalFormula("H4O"));
  TEST_NOT_EQUAL(*e_ptr == copy, true)
}
END_SECTION

START_SECTION(void setNTermLossNames(const std::vector< String > &name))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy);
  vector<String> names;
  names.push_back("Nwaesserchen");
  e_ptr->setNTermLossNames(names);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const std::vector<String>& getNTermLossNames() const)
{
  vector<String> names;
  names.push_back("Nwaesserchen");
  TEST_EQUAL(e_ptr->getNTermLossNames() == names, true)
}
END_SECTION

START_SECTION(void addNTermLossName(const String &name))
{
  Residue copy(*e_ptr);
  TEST_EQUAL(*e_ptr, copy);
  e_ptr->addNTermLossName("Nwaesserchen2");
  TEST_NOT_EQUAL(*e_ptr == copy, true);
}
END_SECTION

START_SECTION(bool hasNTermNeutralLosses() const)
{
  Residue copy(*e_ptr);
  TEST_EQUAL(copy.hasNTermNeutralLosses(), true)
  copy.setNTermLossFormulas(vector<EmpiricalFormula>());
  copy.setNTermLossNames(vector<String>());
  TEST_EQUAL(copy.hasNTermNeutralLosses(), false)
}
END_SECTION

START_SECTION(void setFormula(const EmpiricalFormula &formula))
{
  Residue copy(*e_ptr);
  e_ptr->setFormula(EmpiricalFormula("C2H6O"));
  TEST_NOT_EQUAL(*e_ptr, copy)
  TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), 46.06852164251619);
  TEST_REAL_SIMILAR(e_ptr->getMonoWeight(), 46.0418651914);
}
END_SECTION

START_SECTION(EmpiricalFormula getFormula(ResidueType res_type=Full) const)
  TEST_EQUAL(e_ptr->getFormula(), EmpiricalFormula("C2H6O"))
END_SECTION

START_SECTION(void setAverageWeight(double weight))
{
  Residue copy(*e_ptr);
  e_ptr->setAverageWeight(123.4);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(double getAverageWeight(ResidueType res_type=Full) const)
  TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), 123.4)
END_SECTION

START_SECTION(void setMonoWeight(double weight))
  Residue copy(*e_ptr);
  e_ptr->setMonoWeight(1234.5);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(double getMonoWeight(ResidueType res_type=Full) const)
  TEST_REAL_SIMILAR(e_ptr->getMonoWeight(), 1234.5)
END_SECTION


START_SECTION(void setModification(const String& name))
{
  e_ptr->setOneLetterCode("M"); // we need M for this mod
  TEST_EQUAL(e_ptr->getModificationName(), "")
  TEST_EQUAL(e_ptr->getModification() == nullptr, true)
  e_ptr->setModification("Oxidation");
  TEST_EQUAL(e_ptr->getModificationName(), "Oxidation")
  TEST_EQUAL(e_ptr->getModification()->getFullId(), "Oxidation (M)")
  e_ptr->setOneLetterCode("C"); // TODO Why do we allow setting of the OneLetterCode?? This should be const after construction!
  // TODO Shouldn't the modification be re-checked or deleted if the OneLetterCode changes???
}
END_SECTION

START_SECTION(void setModificationByDiffMonoMass(double diffMonoMass))
{
  e_ptr->setModificationByDiffMonoMass(-1000000);
  TEST_EQUAL(e_ptr->getModificationName(), "")
  TEST_EQUAL(e_ptr->getModification()->getFullId(), "C[-1.0e06]")
  e_ptr->setOneLetterCode("M"); // we need M for the next mod to be findable
  e_ptr->setModificationByDiffMonoMass(15.9949);
  TEST_EQUAL(e_ptr->getModificationName(), "Oxidation")
  TEST_EQUAL(e_ptr->getModification()->getFullId(), "Oxidation (M)")
  e_ptr->setOneLetterCode("B");
}
END_SECTION

START_SECTION(String Residue::toString() const)
{
  auto rr(*db->getResidue("Met"));
  TEST_EQUAL(rr.toString(), "M");
  TEST_EQUAL(rr.getModification() == nullptr, true)
  rr.setModification("Oxidation");
  TEST_EQUAL(rr.getModificationName(), "Oxidation")
  TEST_EQUAL(rr.toString(), "M(Oxidation)");
  const ResidueModification* mod = ResidueModification::createUnknownFromMassString("123", 123.0, false, ResidueModification::ANYWHERE, &rr);
  rr.setModification(mod);
  TEST_EQUAL(rr.toString(), "M[123]");
}
END_SECTION

START_SECTION(const String& getModificationName() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(const ResidueModification* getModification() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setLowMassIons(const std::vector< EmpiricalFormula > &low_mass_ions))
{
  Residue copy(*e_ptr);
  vector<EmpiricalFormula> ions;
  ions.push_back(EmpiricalFormula("NH3"));
  ions.push_back(EmpiricalFormula("PO4"));
  e_ptr->setLowMassIons(ions);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(const std::vector<EmpiricalFormula>& getLowMassIons() const)
  TEST_EQUAL(e_ptr->getLowMassIons()[0], EmpiricalFormula("NH3"))
END_SECTION

START_SECTION(bool hasNeutralLoss() const)
{
  Residue res;
  TEST_EQUAL(res.hasNeutralLoss(), false)
  res.addLossFormula(EmpiricalFormula("H2O"));
  TEST_EQUAL(res.hasNeutralLoss(), true)
}
END_SECTION

START_SECTION(bool operator==(const Residue &residue) const)
{
  Residue r;
  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setName("other_name");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)

  set<String> syns;
  syns.insert("new_syn");
  r.setSynonyms(syns);
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setThreeLetterCode("3lc");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setOneLetterCode("1");
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
}
END_SECTION

START_SECTION(bool operator!=(const Residue &residue) const)
{
  Residue r;
  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setName("other_name");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)

  set<String> syns;
  syns.insert("new_syn");
  r.setSynonyms(syns);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setThreeLetterCode("3lc");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setOneLetterCode("1");
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
}
END_SECTION

START_SECTION(bool operator==(char one_letter_code) const)
  TEST_EQUAL(*e_ptr == 'B', true)
END_SECTION

START_SECTION(bool operator!=(char one_letter_code) const)
  TEST_EQUAL(*e_ptr != 'C', true)
END_SECTION

START_SECTION(void setPka(double value))
  Residue copy(*e_ptr);
  e_ptr->setPka(345.5);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(double getPka() const)
  TEST_REAL_SIMILAR(e_ptr->getPka(), 345.5)
END_SECTION

START_SECTION(void setPkb(double value))
{
  Residue copy(*e_ptr);
  e_ptr->setPkb(675.8);
  TEST_NOT_EQUAL(*e_ptr, copy)
}
END_SECTION

START_SECTION(double getPkb() const)
  TEST_REAL_SIMILAR(e_ptr->getPkb(), 675.8)
END_SECTION

START_SECTION(void setPkc(double value))
  Residue copy(*e_ptr);
  e_ptr->setPkc(9329.0);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(double getPkc() const)
  TEST_REAL_SIMILAR(e_ptr->getPkc(), 9329.0)
END_SECTION

START_SECTION(double getPiValue() const)
  TEST_REAL_SIMILAR(db->getResidue("A")->getPiValue(), 6.11)
END_SECTION

START_SECTION(void setSideChainBasicity(double gb_sc))
  Residue copy(*e_ptr);
  e_ptr->setSideChainBasicity(654.3);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(double getSideChainBasicity() const)
  TEST_REAL_SIMILAR(e_ptr->getSideChainBasicity(), 654.3)
END_SECTION


START_SECTION(void setBackboneBasicityLeft(double gb_bb_l))
  Residue copy(*e_ptr);
  e_ptr->setBackboneBasicityLeft(123.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(double getBackboneBasicityLeft() const)
  TEST_REAL_SIMILAR(e_ptr->getBackboneBasicityLeft(), 123.6)
END_SECTION


START_SECTION(void setBackboneBasicityRight(double gb_bb_r))
  Residue copy(*e_ptr);
  e_ptr->setBackboneBasicityRight(12345.6);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION


START_SECTION(double getBackboneBasicityRight() const)
  TEST_REAL_SIMILAR(e_ptr->getBackboneBasicityRight(), 12345.6)
END_SECTION

START_SECTION(bool isModified() const)
{
  Residue res;
  res.setOneLetterCode("M"); // we need M for this mod
  TEST_EQUAL(res.isModified(), false)
  res.setModification("Oxidation");
  TEST_EQUAL(res.isModified(), true)
}
END_SECTION

START_SECTION((void setResidueSets(const std::set< String > &residues_sets)))
{
  set<String> res_sets;
  res_sets.insert("rs1");
  res_sets.insert("rs2");
  e_ptr->setResidueSets(res_sets);
  TEST_EQUAL(res_sets == e_ptr->getResidueSets(), true)
}
END_SECTION

START_SECTION((void addResidueSet(const String &residue_sets)))
  e_ptr->addResidueSet("rs3");
  TEST_EQUAL(e_ptr->getResidueSets().size(), 3)
END_SECTION

START_SECTION((const std::set<String>& getResidueSets() const))
{
  set<String> res_sets;
  res_sets.insert("rs1");
  res_sets.insert("rs2");
  res_sets.insert("rs3");
  TEST_EQUAL(e_ptr->getResidueSets() == res_sets, true)
}
END_SECTION

START_SECTION((bool isInResidueSet(const String &residue_set)))
{
  TEST_EQUAL(e_ptr->isInResidueSet("rs1"), true)
  TEST_EQUAL(e_ptr->isInResidueSet("rs3"), true)
  TEST_EQUAL(e_ptr->isInResidueSet("rs4"), false)
}
END_SECTION

START_SECTION((static String getResidueTypeName(const ResidueType res_type)))
{
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Full), "full")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Internal), "internal")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::NTerminal), "N-terminal")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::CTerminal), "C-terminal")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::AIon), "a-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::BIon), "b-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::CIon), "c-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::XIon), "x-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::YIon), "y-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::ZIon), "z-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Zp1Ion), "z+1-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Zp2Ion), "z+2-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Precursor), "precursor-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::BIonMinusH20), "b-H2O-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::YIonMinusH20), "y-H2O-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::BIonMinusNH3), "b-NH3-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::YIonMinusNH3), "y-NH3-ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::NonIdentified), "Non-identified ion")
  TEST_STRING_EQUAL(Residue::getResidueTypeName(Residue::Unannotated), "unannotated")


  TEST_EQUAL(Residue::SizeOfResidueType, Residue::names_of_residuetype.size());
  // test that no entries in array are empty (the initializer of std::array<> can have less values than the std::array<>)
  for (int i = 0; i < Residue::names_of_residuetype.size(); ++i)
  {
    TEST_FALSE(Residue::getResidueTypeName(static_cast<Residue::ResidueType>(i)).empty());
  }
}
END_SECTION

delete e_ptr;

END_TEST

