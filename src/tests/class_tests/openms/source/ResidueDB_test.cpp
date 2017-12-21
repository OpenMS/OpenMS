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
	TEST_EQUAL(ptr->hasResidue("LYS"), true)
	TEST_EQUAL(ptr->hasResidue("K"), true)
END_SECTION

START_SECTION(bool hasResidue(const Residue* residue) const)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("BLUBB")), false)
	TEST_EQUAL(ptr->hasResidue(ptr->getResidue("LYS")), true)
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

START_SECTION(void addResidue(const Residue& residue))
	TEST_EQUAL(ptr->hasResidue("UGU"), false)
	TEST_EQUAL(ptr->hasResidue("$"), false)
	TEST_EQUAL(ptr->hasResidue("$hortName"), false)
	TEST_EQUAL(ptr->hasResidue("MyLittleUGUResidue"), false)
	Residue res;
	res.setShortName("$hortName");
	res.setOneLetterCode("$");
	res.setThreeLetterCode("UGU");
	res.setName("MyLittleUGUResidue");
	res.setFormula(EmpiricalFormula("C3H4O4"));
	ptr->addResidue(res);
	TEST_EQUAL(ptr->hasResidue("UGU"), true)
	TEST_EQUAL(ptr->hasResidue("$"), true)
	TEST_EQUAL(ptr->hasResidue("$hortName"), true)
	TEST_EQUAL(ptr->hasResidue("MyLittleUGUResidue"), true)
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
	const Residue* mod_res = nullptr;
    const Residue* mod_res_nullPointer = nullptr;
	mod_res = ptr->getModifiedResidue("Carbamidomethyl (C)");
    TEST_NOT_EQUAL(mod_res, mod_res_nullPointer)
	TEST_EQUAL(ptr->getNumberOfModifiedResidues(), 2)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
