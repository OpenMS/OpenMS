// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/CVTermList.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVTermList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVTermList* ptr = 0;
START_SECTION(CVTermList())
{
	ptr = new CVTermList();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~CVTermList())
{
	delete ptr;
}
END_SECTION

START_SECTION((bool operator==(const CVTermList &cv_term_list) const ))
{
  CVTermList cv_term_list, cv_term_list2;
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
	cv_term_list.setMetaValue("blubb", "blubber");
	TEST_EQUAL(cv_term_list == cv_term_list2, false)
	cv_term_list2.setMetaValue("blubb", "blubber");
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
	CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
	CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	cv_term_list.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list == cv_term_list2, false)
	cv_term_list2.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
}
END_SECTION

START_SECTION((bool operator!=(const CVTermList &cv_term_list) const ))
{
  CVTermList cv_term_list, cv_term_list2;
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  cv_term_list.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
}
END_SECTION

START_SECTION((bool hasCVTerm(const String &accession) const ))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
	CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	CVTermList cv_term_list;
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
	cv_term_list.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
}
END_SECTION

/*
START_SECTION((bool checkCVTerms(const CVMappingRule &rule, const ControlledVocabulary &cv) const ))
{
}
END_SECTION
*/

START_SECTION((void setCVTerms(const std::vector< CVTerm > &terms)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	CVTerm cv_term2("my_accession2", "my_name2", "my_cv_identifier_ref2", "4.0", unit);
	CVTermList cv_term_list;
	vector<CVTerm> cv_terms;
	cv_terms.push_back(cv_term);
	cv_terms.push_back(cv_term2);
	cv_term_list.setCVTerms(cv_terms);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession2"), true);
}
END_SECTION

START_SECTION((const Map<String, std::vector<CVTerm> >& getCVTerms() const ))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession2", "my_name2", "my_cv_identifier_ref2", "4.0", unit);
  CVTermList cv_term_list;
  vector<CVTerm> cv_terms;
  cv_terms.push_back(cv_term);
  cv_terms.push_back(cv_term2);
  cv_term_list.setCVTerms(cv_terms);
  TEST_EQUAL(cv_term_list.getCVTerms().size(), 2);
  TEST_EQUAL(cv_term_list.getCVTerms().has("my_accession"), true);
  TEST_EQUAL(cv_term_list.getCVTerms().has("my_accession2"), true);
}
END_SECTION

START_SECTION((void addCVTerm(const CVTerm &term)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
}
END_SECTION

/*
START_SECTION(([EXTRA] bool checkCVTerms(const ControlledVocabulary &cv) const ))
{
// [Term]
// id: MS:1000132
// name: percent of base peak
// def: "The magnitude of a peak or measurement element expressed in terms of the percentage of the magnitude of the base peak intensity." [PSI:MS]
// xref: value-type:xsd\:float "The allowed value-type for this CV term."
// is_a: MS:1000043 ! intensity unit
  CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)

	ControlledVocabulary cv;
	cv.loadFromOBO("MS", "CV/psi-ms.obo");

	TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)

  CVTerm cv_term2("MS:1000132", "percent of base peaks wrong", "MS", "3.0", unit);
	cv_term_list.addCVTerm(cv_term2);
	TEST_EQUAL(cv_term_list.checkCVTerms(cv), false)
}
END_SECTION

START_SECTION(([EXTRA] void correctCVTermNames()))
{
  CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)

  ControlledVocabulary cv;
  cv.loadFromOBO("MS", "CV/psi-ms.obo");

  TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)

  CVTerm cv_term2("MS:1000132", "percent of base peaks wrong", "MS", "3.0", unit);
  cv_term_list.addCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.checkCVTerms(cv), false)

	cv_term_list.correctCVTermNames();
	TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)	
}
END_SECTION
*/

START_SECTION((CVTermList(const CVTermList &rhs)))
{
  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
	CVTermList cv_term_list2(cv_term_list);
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermList cv_term_list3(cv_term_list);
  TEST_EQUAL(cv_term_list == cv_term_list3, true)
}
END_SECTION

START_SECTION((CVTermList& operator=(const CVTermList &rhs)))
{
  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
  CVTermList cv_term_list2;
	cv_term_list2 = cv_term_list;
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermList cv_term_list3;
	cv_term_list3 = cv_term_list;
  TEST_EQUAL(cv_term_list == cv_term_list3, true)
}
END_SECTION

START_SECTION((bool empty() const))
{
	CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermList cv_term_list;
	TEST_EQUAL(cv_term_list.empty(), true)
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("MS:1000132"), true)
	TEST_EQUAL(cv_term_list.empty(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



