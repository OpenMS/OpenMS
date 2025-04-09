// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <map>

///////////////////////////
#include <OpenMS/METADATA/CVTermList.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVTermList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVTermList* ptr = nullptr;
CVTermList* nullPointer = nullptr;
START_SECTION(CVTermList())
{
  ptr = new CVTermList();
  TEST_NOT_EQUAL(ptr, nullPointer)
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
  TEST_TRUE(cv_term_list == cv_term_list2)
  cv_term_list.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.setMetaValue("blubb", "blubber");
  TEST_TRUE(cv_term_list == cv_term_list2)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.addCVTerm(cv_term);
  TEST_TRUE(cv_term_list == cv_term_list2)
}
END_SECTION

START_SECTION((bool operator!=(const CVTermList &cv_term_list) const ))
{
  CVTermList cv_term_list, cv_term_list2;
  TEST_TRUE(cv_term_list == cv_term_list2)
  cv_term_list.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.setMetaValue("blubb", "blubber");
  TEST_TRUE(cv_term_list == cv_term_list2)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.addCVTerm(cv_term);
  TEST_TRUE(cv_term_list == cv_term_list2)
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
  const auto& t = cv_term_list.getCVTerms();
  TEST_EQUAL(t.size(), 2);
  TEST_EQUAL(t.find("my_accession") != t.end(), true);
  TEST_EQUAL(t.find("my_accession2") != t.end(), true);
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

START_SECTION((void replaceCVTerm(const CVTerm &cv_term)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession").size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession")[0].getValue(), "3.0")
  CVTerm cv_term2("my_accession", "my_name", "my_cv_identifier_ref", "2.0", unit);
  cv_term_list.replaceCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession").size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession")[0].getValue(), "2.0")
}
END_SECTION

START_SECTION((void replaceCVTerms(const std::vector<CVTerm> &cv_terms)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession", "my_name", "my_cv_identifier_ref", "2.0", unit);
  std::vector<CVTerm> tmp;
  tmp.push_back(cv_term);
  tmp.push_back(cv_term2);
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerms(tmp, "my_accession");
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession").size(), 2)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession")[0].getValue(), "3.0")
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession")[1].getValue(), "2.0")
  cv_term_list.replaceCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession").size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession")[0].getValue(), "2.0")
}
END_SECTION

START_SECTION((void replaceCVTerms(const Map<String, vector<CVTerm> >& cv_term_map)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession2", "my_name", "my_cv_identifier_ref", "2.0", unit);
  std::vector<CVTerm> tmp;
  tmp.push_back(cv_term);
  std::vector<CVTerm> tmp2;
  tmp2.push_back(cv_term2);
  std::map<String, std::vector<CVTerm> >new_terms;
  new_terms["my_accession2"] = tmp2;
  TEST_EQUAL(new_terms.find("my_accession2") != new_terms.end(), true);

  // create CVTermList with old "my_accession"
  CVTermList cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerms(tmp, "my_accession");
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession").size(), 1)

  // replace the terms, delete "my_accession" and introduce "my_accession2"
  cv_term_list.replaceCVTerms(new_terms);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession2"), true)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession2").size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms().at("my_accession2")[0].getValue(), "2.0")
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

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((CVTermList(const CVTermList &rhs)))
{
  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
  CVTermList cv_term_list2(cv_term_list);
  TEST_TRUE(cv_term_list == cv_term_list2)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermList cv_term_list3(cv_term_list);
  TEST_TRUE(cv_term_list == cv_term_list3)
}
END_SECTION

START_SECTION((CVTermList(CVTermList &&rhs) noexcept))
{
  // Ensure that CVTermList has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(CVTermList(std::declval<CVTermList&&>())), true)

  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb2", "blubbe");

  CVTermList orig = cv_term_list;
  CVTermList cv_term_list2(std::move(cv_term_list));

  TEST_TRUE(orig == cv_term_list2)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list2.addCVTerm(cv_term);

  orig = cv_term_list2;
  CVTermList cv_term_list3(std::move(cv_term_list2));
  TEST_TRUE(orig == cv_term_list3)
  TEST_EQUAL(cv_term_list3.getCVTerms().size(), 1)
}
END_SECTION

START_SECTION((CVTermList& operator=(const CVTermList &rhs)))
{
  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
  CVTermList cv_term_list2;
  cv_term_list2 = cv_term_list;
  TEST_TRUE(cv_term_list == cv_term_list2)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermList cv_term_list3;
  cv_term_list3 = cv_term_list;
  TEST_TRUE(cv_term_list == cv_term_list3)
}
END_SECTION

START_SECTION((CVTermList& operator=(CVTermList &&rhs)))
{
  CVTermList cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");

  CVTermList orig = cv_term_list;

  CVTermList cv_term_list2;
  cv_term_list2 = std::move(cv_term_list);
  TEST_TRUE(orig == cv_term_list2)
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



