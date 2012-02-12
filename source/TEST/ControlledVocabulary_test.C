// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch  $
// $Authors: Marc Sturm, Andreas Bertsch $
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
ControlledVocabulary* nullPointer = 0;
START_SECTION((ControlledVocabulary()))
	ptr = new ControlledVocabulary();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ControlledVocabulary()))
	delete ptr;
END_SECTION

START_SECTION(const String& name() const)
	TEST_EQUAL(ControlledVocabulary().name(),"")
END_SECTION

ControlledVocabulary cv;
START_SECTION(void loadFromOBO(const String &name, const String &filename))
	cv.loadFromOBO("bla",OPENMS_GET_TEST_DATA_PATH("ControlledVocabulary.obo"));
	TEST_EQUAL(cv.name(),"bla")
END_SECTION

START_SECTION(bool exists(const String& id) const)
	TEST_EQUAL(cv.exists("OpenMS:1"),true)
	TEST_EQUAL(cv.exists("OpenMS:2"),true)
	TEST_EQUAL(cv.exists("OpenMS:3"),true)
	TEST_EQUAL(cv.exists("OpenMS:4"),true)
	TEST_EQUAL(cv.exists("OpenMS:5"),true)
	TEST_EQUAL(cv.exists("OpenMS:6"),true)
	TEST_EQUAL(cv.exists("OpenMS:7"),false)
END_SECTION

START_SECTION(const CVTerm& getTerm(const String& id) const)
	const ControlledVocabulary::CVTerm* term=0;
	//Auto
	term = &(cv.getTerm("OpenMS:1"));
	TEST_EQUAL(term->id,"OpenMS:1")
	TEST_EQUAL(term->name,"Auto")
	TEST_EQUAL(term->description,"Auto desc")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),0)
	TEST_EQUAL(term->unparsed.size(),0)
	TEST_EQUAL(term->synonyms.size(),2)
	TEST_STRING_EQUAL(term->synonyms[0],"Kutsche")
	TEST_STRING_EQUAL(term->synonyms[1],"Karre")
	//Ford
	term = &(cv.getTerm("OpenMS:2"));
	TEST_EQUAL(term->id,"OpenMS:2")
	TEST_EQUAL(term->name,"Ford")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:1")
	TEST_EQUAL(term->unparsed.size(),0)
	TEST_EQUAL(term->synonyms.size(),0)
	//Mercedes
	term = &(cv.getTerm("OpenMS:3"));
	TEST_EQUAL(term->id,"OpenMS:3")
	TEST_EQUAL(term->name,"Mercedes")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:1")
	TEST_STRING_EQUAL(term->synonyms[0],"Zedes")
	//A-Klasse
	term = &(cv.getTerm("OpenMS:4"));
	TEST_EQUAL(term->id,"OpenMS:4")
	TEST_EQUAL(term->name,"A-Klasse")
	TEST_EQUAL(term->description,"A-Klasse desc")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:3")
	TEST_EQUAL(term->unparsed.size(),3)
	TEST_EQUAL(term->unparsed[0],"xref: unparsed line 1")
	TEST_EQUAL(term->unparsed[1],"xref: unparsed line 2")
	TEST_EQUAL(term->unparsed[2],"xref: unparsed line 3")
	TEST_EQUAL(term->synonyms.size(),0)
	//Mustang
	term = &(cv.getTerm("OpenMS:5"));
	TEST_EQUAL(term->id,"OpenMS:5")
	TEST_EQUAL(term->name,"Mustang")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:2")
	TEST_EQUAL(term->unparsed.size(),0)
	TEST_EQUAL(term->synonyms.size(),0)
	//Ka
	term = &(cv.getTerm("OpenMS:6"));
	TEST_EQUAL(term->id,"OpenMS:6")
	TEST_EQUAL(term->name,"Ka")
	TEST_EQUAL(term->description,"Ka desc")
	TEST_EQUAL(term->obsolete,true)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:2")
	TEST_EQUAL(term->unparsed.size(),0)
	TEST_EQUAL(term->synonyms.size(),0)
	
	TEST_EXCEPTION(Exception::InvalidValue , cv.getTerm("OpenMS:7"))
END_SECTION

START_SECTION(bool isChildOf(const String& child, const String& parent) const)
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
END_SECTION

START_SECTION((const Map<String, CVTerm>& getTerms() const))
	Map<String, ControlledVocabulary::CVTerm> terms = cv.getTerms();
	TEST_EQUAL(terms.size(), 6)
	TEST_EQUAL(terms.has("OpenMS:1"), true)
	TEST_EQUAL(terms.has("OpenMS:2"), true)
	TEST_EQUAL(terms.has("OpenMS:3"), true)
	TEST_EQUAL(terms.has("OpenMS:4"), true)
	TEST_EQUAL(terms.has("OpenMS:5"), true)
	TEST_EQUAL(terms.has("OpenMS:6"), true)
	TEST_EQUAL(terms.has("OpenMS:7"), false)

END_SECTION

START_SECTION((void getAllChildTerms(std::set<String>& terms, const String& parent) const))
	set<String> terms;
	cv.getAllChildTerms(terms, "OpenMS:2");
	TEST_EQUAL(terms.size(), 2)
	TEST_EQUAL(terms.find("OpenMS:2") == terms.end(), true)
	TEST_EQUAL(terms.find("OpenMS:5") == terms.end(), false)
END_SECTION


ControlledVocabulary::CVTerm * cvterm = 0;
ControlledVocabulary::CVTerm * cvtermNullPointer = 0;

START_SECTION(([ControlledVocabulary::CVTerm] CVTerm()))
{
  cvterm = new ControlledVocabulary::CVTerm();
  TEST_NOT_EQUAL(cvterm, cvtermNullPointer)
}
END_SECTION


START_SECTION(([ControlledVocabulary::CVTerm] static String getXRefTypeName(XRefType type)))
{  
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_STRING), "xsd:string")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_INTEGER), "xsd:integer")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_DECIMAL), "xsd:decimal")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER), "xsd:negativeInteger")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER), "xsd:positiveInteger")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER), "xsd:nonNegativeInteger")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER), "xsd:nonPositiveInteger")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_BOOLEAN), "xsd:boolean")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_DATE), "xsd:date")
  TEST_STRING_EQUAL(ControlledVocabulary::CVTerm::getXRefTypeName(ControlledVocabulary::CVTerm::XSD_ANYURI), "xsd:anyURI")
}
END_SECTION

START_SECTION(([ControlledVocabulary::CVTerm] CVTerm(const CVTerm &rhs)))
{  
  ControlledVocabulary::CVTerm a;

  a.name = "test_cvterm";
  a.id = "test_id";
  a.parents.insert("test_parent");
  a.children.insert("test_children");
  a.obsolete = true;
  a.description = "test_description";
  a.synonyms = StringList::create("test,synonyms");
  a.unparsed = StringList::create("test,unparsed");
  a.xref_type = ControlledVocabulary::CVTerm::XSD_DECIMAL;
  a.xref_binary = StringList::create("test,xref_binary");
  a.units.insert("units");

  ControlledVocabulary::CVTerm b(a);

  TEST_STRING_EQUAL(b.name,a.name)
  TEST_STRING_EQUAL(b.id,a.id)
  TEST_EQUAL(b.parents == a.parents, true)
  TEST_EQUAL(b.children == a.children, true)
  TEST_EQUAL(b.obsolete, a.obsolete)
  TEST_STRING_EQUAL(b.description,a.description)
  TEST_EQUAL(b.synonyms, a.synonyms)
  TEST_EQUAL(b.unparsed, a.unparsed)
  TEST_EQUAL(b.xref_type == a.xref_type, true)
  TEST_EQUAL(b.xref_binary, a.xref_binary)
  TEST_EQUAL(b.units == a.units, true)
}
END_SECTION

START_SECTION(([ControlledVocabulary::CVTerm] CVTerm& operator=(const CVTerm &rhs)))
{
  ControlledVocabulary::CVTerm a,b;

  a.name = "test_cvterm";
  a.id = "test_id";
  a.parents.insert("test_parent");
  a.children.insert("test_children");
  a.obsolete = true;
  a.description = "test_description";
  a.synonyms = StringList::create("test,synonyms");
  a.unparsed = StringList::create("test,unparsed");
  a.xref_type = ControlledVocabulary::CVTerm::XSD_DECIMAL;
  a.xref_binary = StringList::create("test,xref_binary");
  a.units.insert("units");

  b = a;

  TEST_STRING_EQUAL(b.name,a.name)
  TEST_STRING_EQUAL(b.id,a.id)
  TEST_EQUAL(b.parents == a.parents, true)
  TEST_EQUAL(b.children == a.children, true)
  TEST_EQUAL(b.obsolete, a.obsolete)
  TEST_STRING_EQUAL(b.description,a.description)
  TEST_EQUAL(b.synonyms, a.synonyms)
  TEST_EQUAL(b.unparsed, a.unparsed)
  TEST_EQUAL(b.xref_type == a.xref_type, true)
  TEST_EQUAL(b.xref_binary, a.xref_binary)
  TEST_EQUAL(b.units == a.units, true)

  a = a;

  TEST_STRING_EQUAL(a.name,a.name)
  TEST_STRING_EQUAL(a.id,a.id)
  TEST_EQUAL(a.parents == a.parents, true)
  TEST_EQUAL(a.children == a.children, true)
  TEST_EQUAL(a.obsolete, a.obsolete)
  TEST_STRING_EQUAL(a.description,a.description)
  TEST_EQUAL(a.synonyms, a.synonyms)
  TEST_EQUAL(a.unparsed, a.unparsed)
  TEST_EQUAL(a.xref_type == a.xref_type, true)
  TEST_EQUAL(a.xref_binary, a.xref_binary)
  TEST_EQUAL(a.units == a.units, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
