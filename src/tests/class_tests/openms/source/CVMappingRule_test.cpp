// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>

using namespace OpenMS;
using namespace std;

START_TEST(CVMappingRule, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappingRule* ptr = nullptr;
CVMappingRule* nullPointer = nullptr;
START_SECTION(CVMappingRule())
{
	ptr = new CVMappingRule();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVMappingRule())
{
	delete ptr;
}
END_SECTION

ptr = new CVMappingRule();

START_SECTION((CVMappingRule(const CVMappingRule &rhs)))
{
	CVMappingRule cvmr;
	cvmr.setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(CVMappingRule(cvmr).getIdentifier(), "my_test_identifier")

	cvmr.setElementPath("my_test_elementpath");
	TEST_STRING_EQUAL(CVMappingRule(cvmr).getElementPath(), "my_test_elementpath")

	cvmr.setRequirementLevel(CVMappingRule::MUST);
	TEST_EQUAL(CVMappingRule(cvmr).getRequirementLevel(), CVMappingRule::MUST);
	cvmr.setRequirementLevel(CVMappingRule::SHOULD);
	TEST_EQUAL(CVMappingRule(cvmr).getRequirementLevel(), CVMappingRule::SHOULD);

	cvmr.setCombinationsLogic(CVMappingRule::AND);
	TEST_EQUAL(CVMappingRule(cvmr).getCombinationsLogic(), CVMappingRule::AND);
	cvmr.setCombinationsLogic(CVMappingRule::XOR);
	TEST_EQUAL(CVMappingRule(cvmr).getCombinationsLogic(), CVMappingRule::XOR);

	cvmr.setScopePath("my_test_scopepath");
	TEST_STRING_EQUAL(CVMappingRule(cvmr).getScopePath(), "my_test_scopepath")

	CVMappingTerm term1, term2;
	term1.setAccession("BLA:1");
	term2.setAccession("BLA:2");
	vector<CVMappingTerm> terms;
	terms.push_back(term1);
	terms.push_back(term2);
	cvmr.setCVTerms(terms);
	TEST_EQUAL(CVMappingRule(cvmr).getCVTerms() == terms, true)
}
END_SECTION

START_SECTION((CVMappingRule& operator=(const CVMappingRule &rhs)))
{
  CVMappingRule cvmr, cvmr_copy;
  cvmr.setIdentifier("my_test_identifier");
	cvmr_copy = cvmr;
  TEST_STRING_EQUAL(cvmr_copy.getIdentifier(), "my_test_identifier")

  cvmr.setElementPath("my_test_elementpath");
	cvmr_copy = cvmr;
  TEST_STRING_EQUAL(cvmr_copy.getElementPath(), "my_test_elementpath")

  cvmr.setRequirementLevel(CVMappingRule::MUST);
	cvmr_copy = cvmr;
  TEST_EQUAL(cvmr_copy.getRequirementLevel(), CVMappingRule::MUST);
  cvmr.setRequirementLevel(CVMappingRule::SHOULD);
	cvmr_copy = cvmr;
  TEST_EQUAL(cvmr_copy.getRequirementLevel(), CVMappingRule::SHOULD);

  cvmr.setCombinationsLogic(CVMappingRule::AND);
	cvmr_copy = cvmr;
  TEST_EQUAL(cvmr_copy.getCombinationsLogic(), CVMappingRule::AND);
  cvmr.setCombinationsLogic(CVMappingRule::XOR);
	cvmr_copy = cvmr;
  TEST_EQUAL(cvmr_copy.getCombinationsLogic(), CVMappingRule::XOR);

  cvmr.setScopePath("my_test_scopepath");
	cvmr_copy = cvmr;
  TEST_STRING_EQUAL(cvmr_copy.getScopePath(), "my_test_scopepath")

  CVMappingTerm term1, term2;
  term1.setAccession("BLA:1");
  term2.setAccession("BLA:2");
  vector<CVMappingTerm> terms;
  terms.push_back(term1);
  terms.push_back(term2);
  cvmr.setCVTerms(terms);
	cvmr_copy = cvmr;
  TEST_EQUAL(cvmr_copy.getCVTerms() == terms, true)
}
END_SECTION

START_SECTION((bool operator != (const CVMappingRule& rhs) const))
{
  CVMappingRule cvmr, cvmr_copy;
  cvmr.setIdentifier("my_test_identifier");
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)

  cvmr.setElementPath("my_test_elementpath");
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)

  cvmr.setRequirementLevel(CVMappingRule::MUST); // default
	TEST_EQUAL(cvmr != cvmr_copy, false)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)
  cvmr.setRequirementLevel(CVMappingRule::SHOULD);
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)

  cvmr.setCombinationsLogic(CVMappingRule::AND);
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)
  cvmr.setCombinationsLogic(CVMappingRule::XOR);
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)

  cvmr.setScopePath("my_test_scopepath");
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)

  CVMappingTerm term1, term2;
  term1.setAccession("BLA:1");
  term2.setAccession("BLA:2");
  vector<CVMappingTerm> terms;
  terms.push_back(term1);
  terms.push_back(term2);
  cvmr.setCVTerms(terms);
	TEST_FALSE(cvmr == cvmr_copy)
	cvmr_copy = cvmr;
	TEST_EQUAL(cvmr != cvmr_copy, false)
}
END_SECTION

START_SECTION((bool operator == (const CVMappingRule& rhs) const))
{
  CVMappingRule cvmr, cvmr_copy;
  cvmr.setIdentifier("my_test_identifier");
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)

  cvmr.setElementPath("my_test_elementpath");
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)

  cvmr.setRequirementLevel(CVMappingRule::MUST); // default
  TEST_TRUE(cvmr == cvmr_copy)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)
  cvmr.setRequirementLevel(CVMappingRule::SHOULD);
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)

  cvmr.setCombinationsLogic(CVMappingRule::AND);
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)
  cvmr.setCombinationsLogic(CVMappingRule::XOR);
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)

  cvmr.setScopePath("my_test_scopepath");
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)

  CVMappingTerm term1, term2;
  term1.setAccession("BLA:1");
  term2.setAccession("BLA:2");
  vector<CVMappingTerm> terms;
  terms.push_back(term1);
  terms.push_back(term2);
  cvmr.setCVTerms(terms);
  TEST_EQUAL(cvmr == cvmr_copy, false)
  cvmr_copy = cvmr;
  TEST_TRUE(cvmr == cvmr_copy)
}
END_SECTION

START_SECTION((void setIdentifier(const String &identifier)))
{
  ptr->setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(ptr->getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((const String& getIdentifier() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setElementPath(const String &element_path)))
{
  ptr->setElementPath("my_test_elementpath");
	TEST_STRING_EQUAL(ptr->getElementPath(), "my_test_elementpath")
}
END_SECTION

START_SECTION((const String& getElementPath() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setRequirementLevel(RequirementLevel level)))
{
  ptr->setRequirementLevel(CVMappingRule::MUST);
	TEST_EQUAL(ptr->getRequirementLevel(), CVMappingRule::MUST)
	ptr->setRequirementLevel(CVMappingRule::MAY);
	TEST_EQUAL(ptr->getRequirementLevel(), CVMappingRule::MAY)
}
END_SECTION

START_SECTION((RequirementLevel getRequirementLevel() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setCombinationsLogic(CombinationsLogic combinations_logic)))
{
  ptr->setCombinationsLogic(CVMappingRule::AND);
	TEST_EQUAL(ptr->getCombinationsLogic(), CVMappingRule::AND)
	ptr->setCombinationsLogic(CVMappingRule::XOR);
	TEST_EQUAL(ptr->getCombinationsLogic(), CVMappingRule::XOR)
}
END_SECTION

START_SECTION((CombinationsLogic getCombinationsLogic() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setScopePath(const String &path)))
{
  ptr->setScopePath("my_test_scopepath");
	TEST_STRING_EQUAL(ptr->getScopePath(), "my_test_scopepath")
}
END_SECTION

START_SECTION((const String& getScopePath() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setCVTerms(const std::vector< CVMappingTerm > &cv_terms)))
{
  CVMappingTerm cv_term1, cv_term2;
	cv_term1.setAccession("BLA:1");
	cv_term2.setAccession("BLA:2");
	vector<CVMappingTerm> terms;
	terms.push_back(cv_term1);
	terms.push_back(cv_term2);
	ptr->setCVTerms(terms);
	TEST_EQUAL(ptr->getCVTerms() == terms, true)
}
END_SECTION

START_SECTION((const std::vector<CVMappingTerm>& getCVTerms() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void addCVTerm(const CVMappingTerm &cv_terms)))
{
  TEST_EQUAL(ptr->getCVTerms().size(), 2)
	CVMappingTerm cv_term;
	cv_term.setAccession("BLA:3");
	ptr->addCVTerm(cv_term);
	TEST_EQUAL(ptr->getCVTerms().size(), 3)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



