// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/Enzyme.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(EnzymesDB, "$Id$")

/////////////////////////////////////////////////////////////

EnzymesDB* ptr = 0;
EnzymesDB* nullPointer = 0;
String RKP("(?<=R)(?!P)");
START_SECTION(EnzymesDB* getInstance())
    ptr = EnzymesDB::getInstance();
    TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~EnzymesDB())
    NOT_TESTABLE
END_SECTION

START_SECTION((bool hasEnzyme(const String &name) const))
    TEST_EQUAL(ptr->hasEnzyme("Try"), false)
    TEST_EQUAL(ptr->hasEnzyme("Trypsin"), true)
END_SECTION

START_SECTION((const Enzyme* getEnzyme(const String &name) const))
    TEST_EQUAL(ptr->getEnzyme("Trypsin")->getName(), "Trypsin")
END_SECTION

START_SECTION((bool hasRegEx(const String & cleavage_regex) const))
    TEST_EQUAL(ptr->hasRegEx("(?<=[P])(?!P)"), false)
    TEST_EQUAL(ptr->hasRegEx(RKP), true)
END_SECTION

START_SECTION((const Enzyme* getEnzymeByRegEx(const String & cleavage_regex) const))
    TEST_EQUAL(ptr->getEnzymeByRegEx(RKP)->getName(), "Arg-C")
END_SECTION

START_SECTION(bool hasEnzyme(const Enzyme *enzyme) const)
    TEST_EQUAL(ptr->hasEnzyme(ptr->getEnzyme("Trypsin")), true)
END_SECTION

START_SECTION(void setEnzymes(const String &filename))
    NOT_TESTABLE // this method is hard to test, just provided for convenience
END_SECTION
    
START_SECTION(void addEnzyme(const Enzyme &enzyme))
    TEST_EQUAL(ptr->hasEnzyme("Try"), false)
    TEST_EQUAL(ptr->hasRegEx("(?<=[P])(?!P)"), false)
    Enzyme enzy = Enzyme("Try","(?<=[P])(?!P)");
    ptr->addEnzyme(enzy);
    TEST_EQUAL(ptr->hasEnzyme("Try"), true)
    TEST_EQUAL(ptr->hasRegEx("(?<=[P])(?!P)"), true)
END_SECTION

START_SECTION(EnzymeIterator beginEnzyme())
    EnzymesDB::EnzymeIterator it = ptr->beginEnzyme();
    Size count(0);
    while (it != ptr->endEnzyme())
    {
      ++it;
      ++count;
    }
    TEST_EQUAL(count >= 10, true)
END_SECTION
  
START_SECTION(EnzymeIterator endEnzyme())
    NOT_TESTABLE // tested above
END_SECTION

START_SECTION(EnzymeConstIterator beginEnzyme() const)
    const EnzymesDB* const_ptr = ptr;
    EnzymesDB::EnzymeConstIterator it = const_ptr->beginEnzyme();
    Size count(0);
    while (it != const_ptr->endEnzyme())
    {
      ++it;
      ++count;
    }
    TEST_EQUAL(count >= 10, true)
END_SECTION

START_SECTION(EnzymeConstIterator endEnzyme() const)
    NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void getAllNames(std::vector< String > &all_names)))
    vector<String> names;
    ptr->getAllNames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "Tryptryp") != names.end(), false)
    Size old_size=names.size();
    ptr->getAllNames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

START_SECTION((void getAllXTandemNames(std::vector< String > &all_names)))
    vector<String> names;
    ptr->getAllXTandemNames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "no cleavage") != names.end(), false)
    Size old_size=names.size();
    ptr->getAllXTandemNames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

START_SECTION((void getAllOMSSANames(std::vector< String > &all_names)))
    vector<String> names;
    ptr->getAllOMSSANames(names);
    TEST_EQUAL(find(names.begin(), names.end(), "Trypsin") != names.end(), true)
    TEST_EQUAL(find(names.begin(), names.end(), "leukocyte elastase") != names.end(), false)
    Size old_size=names.size();
    ptr->getAllOMSSANames(names);
    TEST_EQUAL(names.size(), old_size)
END_SECTION

END_TEST
