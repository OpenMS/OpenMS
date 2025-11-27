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
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>

using namespace OpenMS;
using namespace std;

START_TEST(CVMappings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappings* ptr = nullptr;
CVMappings* nullPointer = nullptr;
START_SECTION(CVMappings())
{
	ptr = new CVMappings();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVMappings())
{
	delete ptr;
}
END_SECTION

ptr = new CVMappings();

START_SECTION((CVMappings(const CVMappings &rhs)))
{
	CVMappings cvm;
  CVMappingRule r1, r2;
  vector<CVMappingRule> rules;
  rules.push_back(r1);
  rules.push_back(r2);
	cvm.setMappingRules(rules);
	TEST_EQUAL(CVMappings(cvm).getMappingRules() == rules, true)

	CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
	cvm.setCVReferences(refs);
	TEST_EQUAL(CVMappings(cvm).getCVReferences() == refs, true)
}
END_SECTION

START_SECTION((CVMappings& operator=(const CVMappings &rhs)))
{
  CVMappings cvm, cvm_copy;
  CVMappingRule r1, r2;
  vector<CVMappingRule> rules;
  rules.push_back(r1);
  rules.push_back(r2);
  cvm.setMappingRules(rules);
	cvm_copy = cvm;
  TEST_EQUAL(cvm_copy.getMappingRules() == rules, true)

  CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
  cvm.setCVReferences(refs);
	cvm_copy = cvm;
  TEST_EQUAL(cvm_copy.getCVReferences() == refs, true)
}
END_SECTION

START_SECTION((bool operator == (const CVMappings& rhs) const))
{
  CVMappings cvm, cvm_copy;
  CVMappingRule r1, r2;
  vector<CVMappingRule> rules;
  rules.push_back(r1);
  rules.push_back(r2);
  cvm.setMappingRules(rules);
	TEST_EQUAL(cvm == cvm_copy, false)
  cvm_copy = cvm;
	TEST_TRUE(cvm == cvm_copy)

  CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
  cvm.setCVReferences(refs);
	TEST_EQUAL(cvm == cvm_copy, false)
  cvm_copy = cvm;
	TEST_TRUE(cvm == cvm_copy)
}
END_SECTION

START_SECTION((bool operator != (const CVMappings& rhs) const))
{
  CVMappings cvm, cvm_copy;
  CVMappingRule r1, r2;
  vector<CVMappingRule> rules;
  rules.push_back(r1);
  rules.push_back(r2);
  cvm.setMappingRules(rules);
  TEST_FALSE(cvm == cvm_copy)
  cvm_copy = cvm;
  TEST_EQUAL(cvm != cvm_copy, false)

  CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
  cvm.setCVReferences(refs);
  TEST_FALSE(cvm == cvm_copy)
  cvm_copy = cvm;
  TEST_EQUAL(cvm != cvm_copy, false)
}
END_SECTION

START_SECTION((void setMappingRules(const std::vector< CVMappingRule > &cv_mapping_rules)))
{
  CVMappingRule r1, r2;
	vector<CVMappingRule> rules;
	rules.push_back(r1);
	rules.push_back(r2);
	TEST_EQUAL(ptr->getMappingRules().size(), 0)
	ptr->setMappingRules(rules);
	TEST_EQUAL(ptr->getMappingRules() == rules, true)
}
END_SECTION

START_SECTION((const std::vector<CVMappingRule>& getMappingRules() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void addMappingRule(const CVMappingRule &cv_mapping_rule)))
{
  CVMappingRule r1;
	TEST_EQUAL(ptr->getMappingRules().size(), 2)
	ptr->addMappingRule(r1);
	TEST_EQUAL(ptr->getMappingRules().size(), 3)
}
END_SECTION

START_SECTION((void setCVReferences(const std::vector< CVReference > &cv_references)))
{
  CVReference r1, r2;
	r1.setIdentifier("Ref1");
	r2.setIdentifier("Ref2");
	vector<CVReference> refs;
	refs.push_back(r1);
	refs.push_back(r2);
	TEST_EQUAL(ptr->getCVReferences().size(), 0)
	ptr->setCVReferences(refs);
	TEST_EQUAL(ptr->getCVReferences() == refs, true)
}
END_SECTION

START_SECTION((const std::vector<CVReference>& getCVReferences() const))
	NOT_TESTABLE
END_SECTION

START_SECTION((void addCVReference(const CVReference &cv_reference)))
{
  CVReference r1;
	r1.setIdentifier("Ref2,5");
	TEST_EQUAL(ptr->getCVReferences().size(), 2)
	ptr->addCVReference(r1);
	TEST_EQUAL(ptr->getCVReferences().size(), 3)
}
END_SECTION

START_SECTION((bool hasCVReference(const String &identifier)))
{
 	TEST_EQUAL(ptr->hasCVReference("Ref1"), true)
	TEST_EQUAL(ptr->hasCVReference("Ref2"), true)
	TEST_EQUAL(ptr->hasCVReference("Ref2,5"), true)
	TEST_EQUAL(ptr->hasCVReference("Ref3"), false)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



