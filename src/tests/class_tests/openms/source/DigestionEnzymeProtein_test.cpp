// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
                                     UInt omssa_id))
  DigestionEnzymeProtein copy(e_ptr->getName(), e_ptr->getRegEx(), e_ptr->getSynonyms(), e_ptr->getRegExDescription(), e_ptr->getNTermGain(), e_ptr->getCTermGain(), e_ptr->getPSIID(), e_ptr->getXTandemID(), e_ptr->getOMSSAID());
  TEST_EQUAL(copy.getName(), e_ptr->getName())
  TEST_EQUAL(copy.getRegEx(), e_ptr->getRegEx())
  TEST_EQUAL(copy.getRegExDescription(), e_ptr->getRegExDescription())
  TEST_EQUAL(copy.getNTermGain(), e_ptr->getNTermGain())
  TEST_EQUAL(copy.getCTermGain(), e_ptr->getCTermGain())
  TEST_EQUAL(copy.getPSIID(), e_ptr->getPSIID())
  TEST_EQUAL(copy.getXTandemID(), e_ptr->getXTandemID())
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

END_TEST

