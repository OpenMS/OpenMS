// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

ResidueDB* e_ptr = 0;
CHECK(ResidueDB())
	e_ptr = new ResidueDB;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(virtual ~ResidueDB())
	delete e_ptr;
RESULT

e_ptr = new ResidueDB;

CHECK((ResidueDB(const ResidueDB &residue_db)))
	TEST_EXCEPTION(Exception::NotImplemented, ResidueDB copy(*e_ptr))
RESULT

CHECK((ResidueDB(const String &res_filename, const String &mod_filename) throw (Exception::FileNotFound, Exception::ParseError)))
	TEST_EXCEPTION(Exception::FileNotFound, ResidueDB("FILE_DOES_NOT_EXIST", "FILE_DOES_NOT_EXIST"))
RESULT

CHECK((ResidueDB& operator=(const ResidueDB &aa)))
	ResidueDB copy;
	TEST_EXCEPTION(Exception::NotImplemented, copy = *e_ptr)
RESULT

CHECK((UInt getNumberOfResidues() const))
	TEST_EQUAL(e_ptr->getNumberOfResidues(), 28)
RESULT

CHECK((UInt getNumberOfResidueModifications() const))
	TEST_EQUAL(e_ptr->getNumberOfResidueModifications(), 4)
RESULT

CHECK((const ResidueModification* getModification(const String &name) const))
	TEST_NOT_EQUAL(e_ptr->getModification("phosphorylation"), 0)
	cerr << "." << endl;
	TEST_EQUAL(e_ptr->getModification("phosphorylation")->getName(), "phosphorylation")
RESULT

CHECK((std::set<const ResidueModification*> getModifications(const Residue *residue) const))
  TEST_EQUAL(e_ptr->getModifications(e_ptr->getResidue("C")).size(), 1)
RESULT

CHECK((std::set<const ResidueModification*> getModifications(const String &res_name) const))
 	 TEST_EQUAL(e_ptr->getModifications("C").size(), 1)
RESULT

CHECK((const std::set<const ResidueModification*>& getModifications() const))
  TEST_EQUAL(e_ptr->getModifications().size(), 4)
	TEST_EQUAL(e_ptr->getModifications().size(), e_ptr->getNumberOfResidueModifications())
RESULT

CHECK((const Residue* getResidue(const String &name) const))
  TEST_EQUAL(e_ptr->getResidue("C")->getOneLetterCode(), "C")
RESULT

CHECK((std::set<const Residue*> getResidues(const ResidueModification *modification) const))
  const ResidueModification* mod = e_ptr->getModification("phosphorylation");
	TEST_EQUAL(e_ptr->getResidues(mod).size(), 3)
RESULT

CHECK((std::set<const Residue*> getResidues(const String &mod_name) const))
  TEST_EQUAL(e_ptr->getResidues("phosphorylation").size(), 3)
RESULT

CHECK((const std::set<const Residue*>& getResidues() const))
  TEST_EQUAL(e_ptr->getResidues().size(), 28)
RESULT

CHECK((void setModifications(const String &filename) throw (Exception::FileNotFound, Exception::ParseError)))
	ResidueDB db;
  TEST_EXCEPTION(Exception::FileNotFound, db.setModifications("FILE_DOES_NOT_EXIST"))
RESULT

CHECK((void addResidueModification(ResidueModification modification)))
  ResidueModification new_mod;
	//e_ptr->addResidueModification(new_mod);
	//TEST_EQUAL(e_ptr->getNumberOfResidueModifications(), 4)
	TEST_EXCEPTION(Exception::NotImplemented, e_ptr->addResidueModification(new_mod))
RESULT

CHECK((void setResidues(const String &filename) throw (Exception::FileNotFound, Exception::ParseError)))
	ResidueDB db;
  TEST_EXCEPTION(Exception::FileNotFound, db.setResidues("FILE_DOES_NOT_EXIST"))
RESULT

CHECK((void addResidue(const Residue &residue)))
  Residue new_res;
	e_ptr->addResidue(new_res);
	TEST_EQUAL(e_ptr->getNumberOfResidues(), 29)
RESULT

CHECK((bool hasResidueModification(const String &name) const))
	TEST_EQUAL(e_ptr->hasResidueModification("phosphorylation"), true)
	TEST_EQUAL(e_ptr->hasResidueModification("not_a_real_modification"), false)
RESULT

CHECK((bool hasResidue(const String &name) const))
  TEST_EQUAL(e_ptr->hasResidue("BLUBB"), false)
	TEST_EQUAL(e_ptr->hasResidue("LYS"), true)
	TEST_EQUAL(e_ptr->hasResidue("K"), true)
RESULT

CHECK((ResidueIterator beginResidue()))
	UInt i = 0;
  for (ResidueDB::ResidueIterator it = e_ptr->beginResidue(); it != e_ptr->endResidue(); ++it, ++i)
	{
	}
	TEST_EQUAL(i, 29)
RESULT

CHECK((ResidueIterator endResidue()))
	// see above
RESULT

CHECK((ResidueConstIterator beginResidue() const))
	UInt i = 0;
	const ResidueDB db;
	for (ResidueDB::ResidueConstIterator it = db.beginResidue(); it != db.endResidue(); ++it, ++i)
	{
	}
	TEST_EQUAL(i, 28)
RESULT

CHECK((ResidueConstIterator endResidue() const))
	// see above
RESULT

CHECK((ResidueModificationIterator beginResidueModification()))
  UInt i = 0;
  for (ResidueDB::ResidueModificationIterator it = e_ptr->beginResidueModification(); it != e_ptr->endResidueModification(); ++it, ++i)
  {
  }
  TEST_EQUAL(i, 4)
RESULT

CHECK((ResidueModificationIterator endResidueModification()))
  // see above
RESULT

CHECK((ResidueModificationConstIterator beginResidueModification() const))
  UInt i = 0;
	const ResidueDB db;
  for (ResidueDB::ResidueModificationConstIterator it = db.beginResidueModification(); it != db.endResidueModification(); ++it, ++i)
  {
  }
  TEST_EQUAL(i, 4) 
RESULT

CHECK((ResidueModificationConstIterator endResidueModification() const))
  // see above
RESULT

CHECK((bool operator==(const ResidueDB &rhs) const))
	TEST_EXCEPTION(Exception::NotImplemented, *e_ptr == *e_ptr)
RESULT

CHECK((bool operator!=(const ResidueDB &rhs) const))
	TEST_EXCEPTION(Exception::NotImplemented, *e_ptr == *e_ptr)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
