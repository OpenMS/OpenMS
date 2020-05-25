// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinition* ptr = nullptr;
ModificationDefinition* nullPointer = nullptr;
START_SECTION(ModificationDefinition())
{
  ptr = new ModificationDefinition();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~ModificationDefinition()))
{
  delete ptr;
}
END_SECTION

ptr = new ModificationDefinition();

START_SECTION((ModificationDefinition(const ModificationDefinition& rhs)))
{
  ModificationDefinition mod_def;
  mod_def.setFixedModification(true);
  ModificationDefinition copy(mod_def);
  TEST_EQUAL(mod_def.isFixedModification(), copy.isFixedModification())

  mod_def.setFixedModification(false);
  ModificationDefinition copy2(mod_def);
  TEST_EQUAL(mod_def.isFixedModification(), copy2.isFixedModification())
}
END_SECTION

START_SECTION((ModificationDefinition(const String& mod, bool fixed = true, UInt max_occur = 0)))
{
  ModificationDefinition mod1("Acetyl (N-term)");
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)");
  ModificationDefinition mod2("Oxidation (M)");
  TEST_EQUAL(mod2.getModificationName(), "Oxidation (M)");
  TEST_EQUAL(mod2.isFixedModification(), true);
  TEST_EQUAL(mod2.getMaxOccurrences(), 0);
  ModificationDefinition mod3("Carboxymethyl (C)", false, 2);
  TEST_EQUAL(mod3.getModificationName(), "Carboxymethyl (C)");
  TEST_EQUAL(mod3.isFixedModification(), false);
  TEST_EQUAL(mod3.getMaxOccurrences(), 2);
}
END_SECTION

START_SECTION((ModificationDefinition(const ResidueModification& mod, bool fixed = true, UInt max_occur = 0)))
{
  const ResidueModification res_mod1 = *ModificationsDB::getInstance()->getModification("Acetyl (N-term)");
  ModificationDefinition mod1(res_mod1);
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)");
  TEST_EQUAL(mod1.isFixedModification(), true);
  TEST_EQUAL(mod1.getMaxOccurrences(), 0);
  const ResidueModification res_mod2 = *ModificationsDB::getInstance()->getModification("Oxidation (M)");
  ModificationDefinition mod2(res_mod2, false, 2);
  TEST_EQUAL(mod2.isFixedModification(), false);
  TEST_EQUAL(mod2.getMaxOccurrences(), 2);
}
END_SECTION

START_SECTION((void setFixedModification(bool fixed)))
{
  ptr->setFixedModification(true);
  TEST_EQUAL(ptr->isFixedModification(), true)
  ptr->setFixedModification(false);
  TEST_EQUAL(ptr->isFixedModification(), false)
}
END_SECTION

START_SECTION((bool isFixedModification() const))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setMaxOccurrences(UInt num)))
{
  ptr->setMaxOccurrences(1);
  TEST_EQUAL(ptr->getMaxOccurrences(), 1)
  ptr->setMaxOccurrences(1000);
  TEST_EQUAL(ptr->getMaxOccurrences(), 1000)
}
END_SECTION

START_SECTION((UInt getMaxOccurrences() const))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((String getModificationName() const))
{
  ModificationDefinition mod1;
  mod1.setModification("Acetyl (N-term)");
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)")
  mod1.setModification("Oxidation (M)");
  TEST_EQUAL(mod1.getModificationName(), "Oxidation (M)")
}
END_SECTION

START_SECTION((String getModification() const))
{
  const ResidueModification* rm = ModificationsDB::getInstance()->getModification("Acetyl (N-term)");
  ModificationDefinition mod1;
  mod1.setModification(rm->getFullId());
  TEST_EQUAL(rm, &(mod1.getModification()));
}
END_SECTION

START_SECTION((void setModification(const String& modification)))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ModificationDefinition& operator=(const ModificationDefinition& element)))
{
  ModificationDefinition mod_def;
  mod_def.setFixedModification(true);
  *ptr = mod_def;
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())

  mod_def.setFixedModification(false);
  *ptr = mod_def;
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())
  
}
END_SECTION

START_SECTION((bool operator==(const ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  TEST_EQUAL(m1 == m2, true)
  m1.setFixedModification(false);
  TEST_EQUAL(m1 == m2, false)
  m1.setFixedModification(true);
  m1.setMaxOccurrences(15);
  TEST_EQUAL(m1 == m2, false)
  m1.setMaxOccurrences(0);
  m1.setModification("Oxidation (M)");
  TEST_EQUAL(m1 == m2, false)
  m2.setModification("Oxidation (M)");
  TEST_EQUAL(m1 == m2, true)
}
END_SECTION

START_SECTION((bool operator!=(const ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  TEST_EQUAL(m1 != m2, false)
  m1.setFixedModification(false);
  TEST_EQUAL(m1 != m2, true)
  m1.setFixedModification(true);
  m1.setMaxOccurrences(15);
  TEST_EQUAL(m1 != m2, true)
  m1.setMaxOccurrences(0);
  m1.setModification("Oxidation (M)");
  TEST_EQUAL(m1 != m2, true)
  m2.setModification("Oxidation (M)");
  TEST_EQUAL(m1 != m2, false)
}
END_SECTION

START_SECTION((bool operator<(const OpenMS::ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  m1.setModification("Oxidation (M)");
  m2.setModification("Carboxymethyl (C)");
  TEST_EQUAL(m1 < m2, false)
  TEST_EQUAL(m1 < m1, false)
  TEST_EQUAL(m2 < m1, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



