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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>

///////////////////////////

START_TEST(ControlledVocabulary, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ControlledVocabulary* ptr = 0;
CHECK((ControlledVocabulary()))
	ptr = new ControlledVocabulary();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ControlledVocabulary()))
	delete ptr;
RESULT

CHECK(const String& name() const)
	TEST_EQUAL(ControlledVocabulary().name(),"")
RESULT

ControlledVocabulary cv;
CHECK(void loadFromOBO(const String &name, const String &filename) throw (Exception::FileNotFound, Exception::ParseError))
	cv.loadFromOBO("bla","data/ControlledVocabulary.obo");
	TEST_EQUAL(cv.name(),"bla")
RESULT

CHECK(bool exists(const String& id) const)
	TEST_EQUAL(cv.exists("OpenMS:1"),true)
	TEST_EQUAL(cv.exists("OpenMS:2"),true)
	TEST_EQUAL(cv.exists("OpenMS:3"),true)
	TEST_EQUAL(cv.exists("OpenMS:4"),true)
	TEST_EQUAL(cv.exists("OpenMS:5"),true)
	TEST_EQUAL(cv.exists("OpenMS:6"),true)
	TEST_EQUAL(cv.exists("OpenMS:7"),false)
RESULT

CHECK(const CVTerm& getTerm(const String& id) const throw (Exception::InvalidValue))
	const ControlledVocabulary::CVTerm* term=0;
	//Auto
	term = &(cv.getTerm("OpenMS:1"));
	TEST_EQUAL(term->id,"OpenMS:1")
	TEST_EQUAL(term->name,"Auto")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),0)
	//Ford
	term = &(cv.getTerm("OpenMS:2"));
	TEST_EQUAL(term->id,"OpenMS:2")
	TEST_EQUAL(term->name,"Ford")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(term->parents[0],"OpenMS:1")
	//Mercedes
	term = &(cv.getTerm("OpenMS:3"));
	TEST_EQUAL(term->id,"OpenMS:3")
	TEST_EQUAL(term->name,"Mercedes")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(term->parents[0],"OpenMS:1")
	//A-Klasse
	term = &(cv.getTerm("OpenMS:4"));
	TEST_EQUAL(term->id,"OpenMS:4")
	TEST_EQUAL(term->name,"A-Klasse")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(term->parents[0],"OpenMS:3")
	//Mustang
	term = &(cv.getTerm("OpenMS:5"));
	TEST_EQUAL(term->id,"OpenMS:5")
	TEST_EQUAL(term->name,"Mustang")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(term->parents[0],"OpenMS:2")
	//Ka
	term = &(cv.getTerm("OpenMS:6"));
	TEST_EQUAL(term->id,"OpenMS:6")
	TEST_EQUAL(term->name,"Ka")
	TEST_EQUAL(term->obsolete,true)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(term->parents[0],"OpenMS:2")
	
	TEST_EXCEPTION(Exception::InvalidValue , cv.getTerm("OpenMS:7"))
RESULT

CHECK(bool isChildOf(const String& child, const String& parent) const  throw (Exception::InvalidValue))
	TEST_EQUAL(cv.isChildOf("OpenMS:6","OpenMS:2"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:5","OpenMS:2"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:1"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:3","OpenMS:1"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:4","OpenMS:3"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:1","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:4","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:3"),false)
	TEST_EXCEPTION(Exception::InvalidValue, cv.isChildOf("OpenMS:7","OpenMS:3"))
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
