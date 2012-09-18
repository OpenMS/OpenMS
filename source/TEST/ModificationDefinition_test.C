// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinition* ptr = 0;
ModificationDefinition* nullPointer = 0;
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

START_SECTION((ModificationDefinition(const ModificationDefinition &rhs)))
{
  ModificationDefinition mod_def;
	mod_def.setTermSpecificity(ResidueModification::C_TERM);
	mod_def.setFixedModification(true);
	ModificationDefinition copy(mod_def);
	TEST_EQUAL(mod_def.getTermSpecificity(), copy.getTermSpecificity())
	TEST_EQUAL(mod_def.isFixedModification(), copy.isFixedModification())

	mod_def.setTermSpecificity(ResidueModification::ANYWHERE);
	mod_def.setFixedModification(false);
	ModificationDefinition copy2(mod_def);
	TEST_EQUAL(mod_def.getTermSpecificity(), copy2.getTermSpecificity())
	TEST_EQUAL(mod_def.isFixedModification(), copy2.isFixedModification())
}
END_SECTION

START_SECTION((ModificationDefinition(const String &mod)))
{
	ModificationDefinition mod1("Acetyl (N-term)");
	TEST_EQUAL(mod1.getModification(), "Acetyl (N-term)");
	ModificationDefinition mod2("Oxidation (M)");
	TEST_EQUAL(mod2.getModification(), "Oxidation (M)");
	ModificationDefinition mod3("Carboxymethyl (C)");
	TEST_EQUAL(mod3.getModification(), "Carboxymethyl (C)");	
}
END_SECTION

START_SECTION((void setTermSpecificity(ResidueModification::Term_Specificity pos)))
{
  ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE);
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM);
}
END_SECTION

START_SECTION((ResidueModification::Term_Specificity getTermSpecificity() const ))
{
  ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
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

START_SECTION((bool isFixedModification() const ))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setMaxOccurences(UInt num)))
{
  ptr->setMaxOccurences(1);
	TEST_EQUAL(ptr->getMaxOccurences(), 1)
	ptr->setMaxOccurences(1000);
	TEST_EQUAL(ptr->getMaxOccurences(), 1000)
}
END_SECTION

START_SECTION((UInt getMaxOccurences() const ))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((String getModification() const ))
{
	ModificationDefinition mod1;
	mod1.setModification("Acetyl (N-term)");
	TEST_EQUAL(mod1.getModification(), "Acetyl (N-term)")
  mod1.setModification("Oxidation (M)");
  TEST_EQUAL(mod1.getModification(), "Oxidation (M)")
}
END_SECTION

START_SECTION((void setModification(const String &modification)))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ModificationDefinition& operator=(const ModificationDefinition &element)))
{
  ModificationDefinition mod_def;
  mod_def.setTermSpecificity(ResidueModification::C_TERM);
  mod_def.setFixedModification(true);
	*ptr = mod_def;
  TEST_EQUAL(mod_def.getTermSpecificity(), ptr->getTermSpecificity())
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())

  mod_def.setTermSpecificity(ResidueModification::ANYWHERE);
  mod_def.setFixedModification(false);
  *ptr = mod_def;
  TEST_EQUAL(mod_def.getTermSpecificity(), ptr->getTermSpecificity())
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())
  
}
END_SECTION

START_SECTION((bool operator==(const ModificationDefinition &rhs) const ))
{
  ModificationDefinition m1, m2;
	TEST_EQUAL(m1 == m2, true)
	m1.setFixedModification(false);
	TEST_EQUAL(m1 == m2, false)
	m1.setFixedModification(true);
	m1.setMaxOccurences(15);
	TEST_EQUAL(m1 == m2, false)
	m1.setMaxOccurences(0);
	m1.setModification("Oxidation (M)");
	TEST_EQUAL(m1 == m2, false)
	m2.setModification("Oxidation (M)");
	TEST_EQUAL(m1 == m2, true)
	m1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(m1 == m2, false)
}
END_SECTION

START_SECTION((bool operator!=(const ModificationDefinition &rhs) const ))
{
  ModificationDefinition m1, m2;
  TEST_EQUAL(m1 != m2, false)
  m1.setFixedModification(false);
  TEST_EQUAL(m1 != m2, true)
  m1.setFixedModification(true);
  m1.setMaxOccurences(15);
  TEST_EQUAL(m1 != m2, true)
  m1.setMaxOccurences(0);
  m1.setModification("Oxidation (M)");
  TEST_EQUAL(m1 != m2, true)
  m2.setModification("Oxidation (M)");
  TEST_EQUAL(m1 != m2, false)
	m1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(m1 != m2, true)
}
END_SECTION

START_SECTION((bool operator<(const OpenMS::ModificationDefinition &) const ))
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



