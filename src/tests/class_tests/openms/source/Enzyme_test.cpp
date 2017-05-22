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

#include <OpenMS/CHEMISTRY/Enzyme.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Enzyme, "$Id$")

/////////////////////////////////////////////////////////////

Enzyme* e_ptr = 0;
/*
Enzyme* e_nullPointer = 0;

START_SECTION((Enzyme()))
	e_ptr = new Enzyme();
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION
*/

START_SECTION((virtual ~Enzyme()))
	delete e_ptr;
END_SECTION

EnzymesDB* db = EnzymesDB::getInstance();
e_ptr = new Enzyme(*db->getEnzyme("Trypsin"));

String RKP("(?<=[RKP])(?!P)");

START_SECTION(Enzyme(const Enzyme &enzyme))
	Enzyme copy(*e_ptr);
	TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(Enzyme(const String &name, 
                       const String &cleavage_regex,
                       const std::set<String> & synonyms,
                       String regex_description,
                       EmpiricalFormula n_term_gain,
                       EmpiricalFormula c_term_gain,
                       String psi_id,
                       String xtandem_id,
                       UInt omssa_id))
	Enzyme copy(e_ptr->getName(), e_ptr->getRegEx(), e_ptr->getSynonyms(), e_ptr->getRegExDescription(), e_ptr->getNTermGain(), e_ptr->getCTermGain(), e_ptr->getPSIid(), e_ptr->getXTANDEMid(), e_ptr->getOMSSAid());
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getRegEx(), e_ptr->getRegEx())
        TEST_EQUAL(copy.getRegExDescription(), e_ptr->getRegExDescription())
        TEST_EQUAL(copy.getNTermGain(), e_ptr->getNTermGain())
        TEST_EQUAL(copy.getCTermGain(), e_ptr->getCTermGain())
        TEST_EQUAL(copy.getPSIid(), e_ptr->getPSIid())
        TEST_EQUAL(copy.getXTANDEMid(), e_ptr->getXTANDEMid())
        TEST_EQUAL(copy.getOMSSAid(), e_ptr->getOMSSAid())
END_SECTION

START_SECTION(Enzyme& operator=(const Enzyme &enzyme))
	Enzyme copy("","");
	copy = *e_ptr;
	TEST_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(void setName(const String &name))
	Enzyme copy(*e_ptr);
	e_ptr->setName("PepsinA");
	TEST_NOT_EQUAL(copy, *e_ptr)
END_SECTION

START_SECTION(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), "PepsinA")
END_SECTION

START_SECTION(void setSynonyms(const std::set< String > &synonyms))
	Enzyme copy(*e_ptr);
	set<String> syn;
	syn.insert("BLI");
	syn.insert("BLA");
	e_ptr->setSynonyms(syn);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION 

START_SECTION(void addSynonym(const String &synonym))
	Enzyme copy(*e_ptr);
	e_ptr->addSynonym("Tryp");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const std::set<String>& getSynonyms() const)
	TEST_EQUAL(e_ptr->getSynonyms().size(), 3)
END_SECTION

START_SECTION(void setRegEx(const String & cleavage_regex))
	Enzyme copy(*e_ptr);
	e_ptr->setRegEx(RKP);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(const String& getRegEx() const)
	TEST_EQUAL(e_ptr->getRegEx(), RKP)
END_SECTION

START_SECTION(void setRegExDescription(String value))
	Enzyme copy(*e_ptr);
	e_ptr->setRegExDescription("cutting after R K unless followed by P");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getRegExDescription() const)
	TEST_EQUAL(e_ptr->getRegExDescription(), "cutting after R K unless followed by P")
END_SECTION

START_SECTION(void setNTermGain(EmpiricalFormula value))
	Enzyme copy(*e_ptr);
	e_ptr->setNTermGain(EmpiricalFormula("H2"));
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(EmpiricalFormula getNTermGain() const)
	TEST_EQUAL(e_ptr->getNTermGain(), EmpiricalFormula("H2"))
END_SECTION

START_SECTION(void setCTermGain(EmpiricalFormula value))
	Enzyme copy(*e_ptr);
	e_ptr->setCTermGain(EmpiricalFormula("OH2"));
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(EmpiricalFormula getCTermGain() const)
	TEST_EQUAL(e_ptr->getCTermGain(), EmpiricalFormula("OH2"))
END_SECTION

START_SECTION(void setPSIid(String value))
	Enzyme copy(*e_ptr);
	e_ptr->setPSIid("MS:000");
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getPSIid() const)
	TEST_EQUAL(e_ptr->getPSIid(), "MS:000")
END_SECTION

START_SECTION(void setXTANDEMid(String value))
        Enzyme copy(*e_ptr);
        e_ptr->setXTANDEMid("[]|[]");
        TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(String getXTANDEMid() const)
        TEST_EQUAL(e_ptr->getXTANDEMid(), "[]|[]")
END_SECTION

START_SECTION(void setOMSSAid(UInt value))
	Enzyme copy(*e_ptr);
	e_ptr->setOMSSAid(2);
	TEST_NOT_EQUAL(*e_ptr, copy)
END_SECTION

START_SECTION(UInt getOMSSAid() const)
	TEST_EQUAL(e_ptr->getOMSSAid(), 2)
END_SECTION

START_SECTION(bool operator==(const Enzyme &enzyme) const)
	Enzyme r("","");
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
	r.setPSIid("new id");
	TEST_EQUAL(r == *e_ptr, false)
	
	r = *e_ptr;
	TEST_EQUAL(r == *e_ptr, true)
	r.setOMSSAid(-2);
	TEST_EQUAL(r == *e_ptr, false)
END_SECTION
    
START_SECTION(bool operator!=(const Enzyme &enzyme) const)
    Enzyme r("","");
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
	r.setPSIid("new id");
	TEST_EQUAL(r != *e_ptr, true)
	
	r = *e_ptr;
	TEST_EQUAL(r != *e_ptr, false)
	r.setOMSSAid(4);
	TEST_EQUAL(r != *e_ptr, true)
END_SECTION

START_SECTION(bool operator==(String cleavage_regex) const)
	TEST_EQUAL(*e_ptr == RKP, true)
END_SECTION
    
START_SECTION(bool operator!=(String cleavage_regex) const)
	TEST_EQUAL(*e_ptr != "?<=[P]", true)
END_SECTION

END_TEST

