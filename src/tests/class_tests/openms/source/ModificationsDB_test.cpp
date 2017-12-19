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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <limits>
#include <algorithm>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct ResidueModificationOriginCmp
{
  bool operator() (const ResidueModification* a, const ResidueModification* b)
  {
    return a->getOrigin() < b->getOrigin();
  }
};

START_TEST(ModificationsDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(bool ModificationsDB::isInstantiated())
{
  bool instantiated = ModificationsDB::isInstantiated();
  TEST_EQUAL(instantiated, false);
}
END_SECTION

ModificationsDB* ptr = nullptr;
ModificationsDB* nullPointer = nullptr;

START_SECTION(ModificationsDB* getInstance())
{
	ptr = ModificationsDB::getInstance();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(bool ModificationsDB::isInstantiated())
{
  bool instantiated = ModificationsDB::isInstantiated();
  TEST_EQUAL(instantiated, true);
}
END_SECTION

START_SECTION(Size getNumberOfModifications() const)
	// range because data may change over time
	TEST_EQUAL(ptr->getNumberOfModifications() > 100, true);
END_SECTION

START_SECTION(const ResidueModification& getModification(Size index) const)
	TEST_EQUAL(ptr->getModification(0).getId().size() > 0, true)
END_SECTION

  START_SECTION((void searchModifications(std::set<const ResidueModification*>& mods, const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  set<const ResidueModification*> mods;
  ptr->searchModifications(mods, "Phosphorylation", "T", ResidueModification::ANYWHERE);
  TEST_EQUAL(mods.size(), 1);
  TEST_STRING_EQUAL((*mods.begin())->getFullId(), "Phospho (T)");
  // terminal mod:
  ptr->searchModifications(mods, "NIC", "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 1);

  ptr->searchModifications(mods, "Label:18O(1)");

  TEST_EQUAL(mods.size(), 4);
  ABORT_IF(mods.size() != 4);

  // Create a sorted set (sorted by getOrigin) instead of sorted by pointer
  // value -> this is more robust on different platforms.
  set<const ResidueModification*, ResidueModificationOriginCmp> mods_sorted;

  // Copy the mods into our sorted container
  set<const ResidueModification*>::const_iterator mod_it_ = mods.begin();
  for (; mod_it_ != mods.end(); mod_it_++)
  {
    mods_sorted.insert((*mod_it_));
  }

  set<const ResidueModification*, ResidueModificationOriginCmp>::const_iterator mod_it = mods_sorted.begin();

  TEST_EQUAL((*mod_it)->getOrigin(), 'S')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'T')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'Y')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)

  ptr->searchModifications(mods, "Label:18O(1)", "", ResidueModification::C_TERM);

  TEST_EQUAL(mods.size(), 1)
  ABORT_IF(mods.size() != 1)

  // Copy the mods into our sorted container
  mods_sorted.clear();
  for (mod_it_ = mods.begin(); mod_it_ != mods.end(); mod_it_++)
  {
    mods_sorted.insert((*mod_it_));
  }

  mod_it = mods_sorted.begin();
  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)

  // no match, thus mods should be empty
  ptr->searchModifications(mods, "Label:18O(1)", "", ResidueModification::N_TERM);

  TEST_EQUAL(mods.size(), 0)
}
END_SECTION


START_SECTION((void searchModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)))
{
  vector<String> mods;
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 0.1, "S");
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1, "S");
  TEST_EQUAL(mods.size(), 0);

  // terminal mod:
  ptr->searchModificationsByDiffMonoMass(mods, 42, 0.1, "", ResidueModification::N_TERM);
  set<String> uniq_mods;
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }
  TEST_EQUAL(mods.size(), 18);
  TEST_EQUAL(uniq_mods.size(), 18);
  TEST_EQUAL(uniq_mods.find("Acetyl (N-term)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 4200000, 0.1, "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 0);

  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 0.1);
  uniq_mods.clear();
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }

  TEST_EQUAL(uniq_mods.find("Phospho (S)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Phospho (T)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Phospho (Y)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Sulfo (S)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1);
  TEST_EQUAL(mods.size(), 0);

  // make sure the common ones are also found for integer masses (this is how
  // integer mass search is done)
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "S");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (S)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "T");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (T)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "Y");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (Y)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 16.0, 1.0, "M");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Oxidation (M)")
  ptr->searchModificationsByDiffMonoMass(mods, 1.0, 1.0, "N");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Deamidated (N)")
  ptr->searchModificationsByDiffMonoMass(mods, 1.0, 1.0, "Q");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Deamidated (Q)")
}
END_SECTION

START_SECTION((const ResidueModification& getModification(const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  TEST_EQUAL(ptr->getModification("Carboxymethyl (C)").getFullId(), "Carboxymethyl (C)");
  TEST_EQUAL(ptr->getModification("Carboxymethyl (C)").getId(), "Carboxymethyl");

  TEST_EQUAL(ptr->getModification("Phosphorylation", "S", ResidueModification::ANYWHERE).getId(), "Phospho");
  TEST_EQUAL(ptr->getModification("Phosphorylation", "S", ResidueModification::ANYWHERE).getFullId(), "Phospho (S)");

  // terminal mod:
  TEST_EQUAL(ptr->getModification("NIC", "", ResidueModification::N_TERM).getId(), "NIC");
  TEST_EQUAL(ptr->getModification("NIC", "", ResidueModification::N_TERM).getFullId(), "NIC (N-term)");
  TEST_EQUAL(ptr->getModification("Acetyl", "", ResidueModification::N_TERM).getFullId(), "Acetyl (N-term)");
}
END_SECTION

START_SECTION((Size findModificationIndex(const String& mod_name) const))
{
  Size index = numeric_limits<Size>::max();
  index = ptr->findModificationIndex("Phospho (T)");
  TEST_NOT_EQUAL(index, numeric_limits<Size>::max());
}
END_SECTION

START_SECTION(void readFromOBOFile(const String& filename))
	// implicitely tested above
	NOT_TESTABLE
END_SECTION

START_SECTION(void readFromUnimodXMLFile(const String& filename))
	// just provided for convenience at the moment
	NOT_TESTABLE
END_SECTION

START_SECTION((void getAllSearchModifications(std::vector<String>& modifications)))
{
  vector<String> mods;
  ptr->getAllSearchModifications(mods);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "NIC (N-term)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho") != mods.end(), false);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Dehydrated (N-term C)") != mods.end(), true);

  // repeat search .. return size should be the same
  Size old_size=mods.size();
  ptr->getAllSearchModifications(mods);
  TEST_EQUAL(mods.size(), old_size);
}
END_SECTION

START_SECTION((bool addModification(ResidueModification* modification)))
{
  TEST_EQUAL(ptr->has("Phospho (E)"), false);
  ResidueModification* modification = new ResidueModification();
  modification->setFullId("Phospho (E)");
  ptr->addModification(modification);
  TEST_EQUAL(ptr->has("Phospho (E)"), true);
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



