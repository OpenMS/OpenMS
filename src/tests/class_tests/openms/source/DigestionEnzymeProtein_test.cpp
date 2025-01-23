// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DigestionEnzymeProtein, "$Id$")

/////////////////////////////////////////////////////////////

DigestionEnzymeProtein* e_ptr = nullptr;
DigestionEnzymeProtein* e_null = nullptr;

START_SECTION((DigestionEnzymeProtein()))
  e_ptr = new DigestionEnzymeProtein();
  TEST_NOT_EQUAL(e_ptr, e_null)
END_SECTION

START_SECTION((virtual ~DigestionEnzymeProtein()))
  delete e_ptr;
END_SECTION

ProteaseDB* db = ProteaseDB::getInstance();
e_ptr = new DigestionEnzymeProtein(*db->getEnzyme("Trypsin"));

String RKP("(?<=[RKP])(?!P)");

START_SECTION(DigestionEnzymeProtein(const DigestionEnzymeProtein& enzyme))
  DigestionEnzymeProtein copy(*e_ptr);
  TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(DigestionEnzymeProtein(const String& name,
                                     const String& cleavage_regex,
                                     const std::set<String> & synonyms,
                                     String regex_description,
                                     EmpiricalFormula n_term_gain,
                                     EmpiricalFormula c_term_gain,
                                     String psi_id,
                                     String xtandem_id,
                                     Int comet_id,
                                     Int msgf_id,
                                     Int omssa_id))
  DigestionEnzymeProtein copy(e_ptr->getName(), e_ptr->getRegEx(), e_ptr->getSynonyms(), e_ptr->getRegExDescription(), e_ptr->getNTermGain(), e_ptr->getCTermGain(), e_ptr->getPSIID(), e_ptr->getXTandemID(), e_ptr->getCometID(),  e_ptr->getMSGFID(), e_ptr->getOMSSAID());
  TEST_EQUAL(copy.getName(), e_ptr->getName())
  TEST_EQUAL(copy.getRegEx(), e_ptr->getRegEx())
  TEST_EQUAL(copy.getRegExDescription(), e_ptr->getRegExDescription())
  TEST_EQUAL(copy.getNTermGain(), e_ptr->getNTermGain())
  TEST_EQUAL(copy.getCTermGain(), e_ptr->getCTermGain())
  TEST_EQUAL(copy.getPSIID(), e_ptr->getPSIID())
  TEST_EQUAL(copy.getXTandemID(), e_ptr->getXTandemID())
  TEST_EQUAL(copy.getCometID(), e_ptr->getCometID())
  TEST_EQUAL(copy.getMSGFID(), e_ptr->getMSGFID())
  TEST_EQUAL(copy.getOMSSAID(), e_ptr->getOMSSAID())
END_SECTION

START_SECTION(DigestionEnzymeProtein& operator=(const DigestionEnzymeProtein& enzyme))
  DigestionEnzymeProtein copy("","");
  copy = *e_ptr;
  TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(void setName(const String& name))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setName("PepsinA");
  TEST_NOT_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(const String& getName() const)
  TEST_EQUAL(e_ptr->getName(), "PepsinA")
END_SECTION

START_SECTION(void setSynonyms(const std::set<String>& synonyms))
  DigestionEnzymeProtein copy(*e_ptr);
  set<String> syn;
  syn.insert("BLI");
  syn.insert("BLA");
  e_ptr->setSynonyms(syn);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(void addSynonym(const String& synonym))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->addSynonym("Tryp");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
  TEST_EQUAL(e_ptr->getSynonyms().size(), 3)
END_SECTION

START_SECTION(void setRegEx(const String& cleavage_regex))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setRegEx(RKP);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getRegEx() const)
  TEST_EQUAL(e_ptr->getRegEx(), RKP)
END_SECTION

START_SECTION(void setRegExDescription(String value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setRegExDescription("cutting after R K unless followed by P");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getRegExDescription() const)
  TEST_EQUAL(e_ptr->getRegExDescription(), "cutting after R K unless followed by P")
END_SECTION

START_SECTION(void setNTermGain(EmpiricalFormula value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setNTermGain(EmpiricalFormula("H2"));
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(EmpiricalFormula getNTermGain() const)
  TEST_EQUAL(e_ptr->getNTermGain(), EmpiricalFormula("H2"))
END_SECTION

START_SECTION(void setCTermGain(EmpiricalFormula value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setCTermGain(EmpiricalFormula("OH2"));
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(EmpiricalFormula getCTermGain() const)
  TEST_EQUAL(e_ptr->getCTermGain(), EmpiricalFormula("OH2"))
END_SECTION

START_SECTION(void setPSIID(String value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setPSIID("MS:000");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getPSIID() const)
	TEST_EQUAL(e_ptr->getPSIID(), "MS:000")
END_SECTION

START_SECTION(void setXTandemID(String value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setXTandemID("[]|[]");
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getXTandemID() const)
  TEST_EQUAL(e_ptr->getXTandemID(), "[]|[]")
END_SECTION

START_SECTION(void setOMSSAID(UInt value))
  DigestionEnzymeProtein copy(*e_ptr);
  e_ptr->setOMSSAID(2);
  TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(UInt getOMSSAID() const)
  TEST_EQUAL(e_ptr->getOMSSAID(), 2)
END_SECTION

START_SECTION(bool operator==(const DigestionEnzymeProtein& enzyme) const)
  DigestionEnzymeProtein r("","");
  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setName("other_name");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setRegEx("?<=[P]");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  set<String> syns;
  syns.insert("new_syn");
  r.setSynonyms(syns);
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setRegExDescription("new description");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setNTermGain(EmpiricalFormula("H2O"));
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setCTermGain(EmpiricalFormula("H6O"));
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setPSIID("new id");
  TEST_EQUAL(r == *e_ptr, false)

  r = *e_ptr;
  TEST_EQUAL(r == *e_ptr, true)
  r.setOMSSAID(-2);
  TEST_EQUAL(r == *e_ptr, false)
END_SECTION

START_SECTION(bool operator!=(const DigestionEnzymeProtein& enzyme) const)
  DigestionEnzymeProtein r("","");
  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setName("other_name");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setRegEx("?<=[P]");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  set<String> syns;
  syns.insert("new_syn");
  r.setSynonyms(syns);
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setRegExDescription("new description");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setNTermGain(EmpiricalFormula("H2O"));
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setCTermGain(EmpiricalFormula("O"));
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setPSIID("new id");
  TEST_EQUAL(r != *e_ptr, true)

  r = *e_ptr;
  TEST_EQUAL(r != *e_ptr, false)
  r.setOMSSAID(4);
  TEST_EQUAL(r != *e_ptr, true)
END_SECTION

START_SECTION(bool operator==(String cleavage_regex) const)
  TEST_EQUAL(*e_ptr == RKP, true)
END_SECTION

START_SECTION(bool operator!=(String cleavage_regex) const)
  TEST_EQUAL(*e_ptr != "?<=[P]", true)
END_SECTION

delete e_ptr;

END_TEST

