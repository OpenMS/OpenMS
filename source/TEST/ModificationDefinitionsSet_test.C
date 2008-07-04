// -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinitionsSet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinitionsSet* ptr = 0;
CHECK(ModificationDefinitionsSet())
{
	ptr = new ModificationDefinitionsSet();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK((ModificationDefinitionsSet(const ModificationDefinitionsSet &rhs)))
{
  // TODO
}
RESULT

CHECK((ModificationDefinitionsSet(const String &fixed_modifications, const String &variable_modifications="")))
{
  // TODO
}
RESULT

CHECK((virtual ~ModificationDefinitionsSet()))
{
  // TODO
}
RESULT

CHECK((void setMaxModifications(UInt max_mod)))
{
  // TODO
}
RESULT

CHECK((UInt getMaxModifications() const ))
{
  // TODO
}
RESULT

CHECK((UInt getNumberOfModifications() const ))
{
  // TODO
}
RESULT

CHECK((UInt getNumberOfFixedModifications() const ))
{
  // TODO
}
RESULT

CHECK((UInt getNumberOfVariableModifications() const ))
{
  // TODO
}
RESULT

CHECK((void addModification(const ModificationDefinition& mod_def)))
{
  // TODO
}
RESULT

CHECK((void setModifications(const std::set<ModificationDefinition>& mod_defs)))
{
  // TODO
}
RESULT

CHECK((void setModifications(const String& fixed_modifications, const String& variable_modifications)))
{
  // TODO
}
RESULT

CHECK((std::set<ModificationDefinition> getModifications() const ))
{
  // TODO
}
RESULT

CHECK((std::set<String> getModificationNames() const ))
{
  // TODO
}
RESULT

CHECK((std::set<String> getFixedModificationNames() const ))
{
  // TODO
}
RESULT

CHECK((std::set<String> getVariableModificationNames() const ))
{
  // TODO
}
RESULT

CHECK((ModificationDefinitionsSet& operator=(const ModificationDefinitionsSet& element)))
{
  // TODO
}
RESULT

CHECK((bool operator==(const ModificationDefinitionsSet& rhs) const))
{
  // TODO
}
RESULT

CHECK((bool operator!=(const ModificationDefinitionsSet& rhs) const))
{
  // TODO
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



