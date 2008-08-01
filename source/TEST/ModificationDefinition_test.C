// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch$
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
CHECK(ModificationDefinition())
{
	ptr = new ModificationDefinition();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK((virtual ~ModificationDefinition()))
{
	delete ptr;
}
RESULT

ptr = new ModificationDefinition();

CHECK((ModificationDefinition(const ModificationDefinition &rhs)))
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
RESULT

CHECK((ModificationDefinition(const String &mod)))
{
	ModificationDefinition mod1("MOD:00323");
	TEST_EQUAL(mod1.getModification(), "MOD:00323");
	ModificationDefinition mod2("MOD:00719");
	TEST_EQUAL(mod2.getModification(), "MOD:00719");
	ModificationDefinition mod3("MOD:01061");
	TEST_EQUAL(mod3.getModification(), "MOD:01061");	
}
RESULT

CHECK((void setTermSpecificity(ResidueModification::Term_Specificity pos)))
{
  ptr->setTermSpecificity(ResidueModification::ANYWHERE);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::ANYWHERE);
	ptr->setTermSpecificity(ResidueModification::C_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::C_TERM);
}
RESULT

CHECK((ResidueModification::Term_Specificity getTermSpecificity() const ))
{
  ptr->setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(ptr->getTermSpecificity(), ResidueModification::N_TERM)
}
RESULT

CHECK((void setFixedModification(bool fixed)))
{
  ptr->setFixedModification(true);
	TEST_EQUAL(ptr->isFixedModification(), true)
	ptr->setFixedModification(false);
	TEST_EQUAL(ptr->isFixedModification(), false)
}
RESULT

CHECK((bool isFixedModification() const ))
{
	// tested above
  NOT_TESTABLE
}
RESULT

CHECK((void setMaxOccurences(UInt num)))
{
  ptr->setMaxOccurences(1);
	TEST_EQUAL(ptr->getMaxOccurences(), 1)
	ptr->setMaxOccurences(1000);
	TEST_EQUAL(ptr->getMaxOccurences(), 1000)
}
RESULT

CHECK((UInt getMaxOccurences() const ))
{
	// tested above
  NOT_TESTABLE
}
RESULT

CHECK((String getModification() const ))
{
	ModificationDefinition mod1;
	mod1.setModification("MOD:00323");
	TEST_EQUAL(mod1.getModification(), "MOD:00323")
  mod1.setModification("MOD:00719");
  TEST_EQUAL(mod1.getModification(), "MOD:00719")
}
RESULT

CHECK((void setModification(const String &modification)))
{
	// tested above
  NOT_TESTABLE
}
RESULT

CHECK((ModificationDefinition& operator=(const ModificationDefinition &element)))
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
RESULT

CHECK((bool operator==(const ModificationDefinition &rhs) const ))
{
  ModificationDefinition m1, m2;
	TEST_EQUAL(m1 == m2, true)
	m1.setFixedModification(false);
	TEST_EQUAL(m1 == m2, false)
	m1.setFixedModification(true);
	m1.setMaxOccurences(15);
	TEST_EQUAL(m1 == m2, false)
	m1.setMaxOccurences(0);
	m1.setModification("MOD:00720");
	TEST_EQUAL(m1 == m2, false)
	m2.setModification("MOD:00720");
	TEST_EQUAL(m1 == m2, true)
	m1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(m1 == m2, false)
}
RESULT

CHECK((bool operator!=(const ModificationDefinition &rhs) const ))
{
  ModificationDefinition m1, m2;
  TEST_EQUAL(m1 != m2, false)
  m1.setFixedModification(false);
  TEST_EQUAL(m1 != m2, true)
  m1.setFixedModification(true);
  m1.setMaxOccurences(15);
  TEST_EQUAL(m1 != m2, true)
  m1.setMaxOccurences(0);
  m1.setModification("MOD:00720");
  TEST_EQUAL(m1 != m2, true)
  m2.setModification("MOD:00720");
  TEST_EQUAL(m1 != m2, false)
	m1.setTermSpecificity(ResidueModification::N_TERM);
	TEST_EQUAL(m1 != m2, true)
}
RESULT

CHECK((bool operator<(const OpenMS::ModificationDefinition &) const ))
{
  ModificationDefinition m1, m2;
	m1.setModification("MOD:00720");
	m2.setModification("MOD:00719");
	TEST_EQUAL(m1 < m2, false)
	TEST_EQUAL(m1 < m1, false)
	TEST_EQUAL(m2 < m1, true)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



