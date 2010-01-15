// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVMappings, "$Id: CVMappings_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappings* ptr = 0;
START_SECTION(CVMappings())
{
	ptr = new CVMappings();
	TEST_NOT_EQUAL(ptr, 0)
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
	TEST_EQUAL(cvm == cvm_copy, true)

  CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
  cvm.setCVReferences(refs);
	TEST_EQUAL(cvm == cvm_copy, false)
  cvm_copy = cvm;
	TEST_EQUAL(cvm == cvm_copy, true)
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
  TEST_EQUAL(cvm != cvm_copy, true)
  cvm_copy = cvm;
  TEST_EQUAL(cvm != cvm_copy, false)

  CVReference ref1, ref2;
  ref1.setIdentifier("Ref1");
  ref2.setIdentifier("Ref2");
  vector<CVReference> refs;
  refs.push_back(ref1);
  refs.push_back(ref2);
  cvm.setCVReferences(refs);
  TEST_EQUAL(cvm != cvm_copy, true)
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



