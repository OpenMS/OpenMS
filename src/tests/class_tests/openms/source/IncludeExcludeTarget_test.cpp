// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
///////////////////////////

#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;
using namespace std;

START_TEST(IncludeExcludeTarget, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IncludeExcludeTarget* ptr = nullptr;
IncludeExcludeTarget* null_ptr = nullptr;
START_SECTION(IncludeExcludeTarget())
{
	ptr = new IncludeExcludeTarget();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IncludeExcludeTarget())
{
	delete ptr;
}
END_SECTION

START_SECTION((IncludeExcludeTarget(const IncludeExcludeTarget &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~IncludeExcludeTarget()))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorMZ(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductMZ(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getProductMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getProductCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setInterpretations(const std::vector< CVTermList > &interpretations)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getInterpretations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInterpretation(const CVTermList &interpretation)))
{
  // TODO
}
END_SECTION

START_SECTION((void setConfigurations(const std::vector< Configuration > &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Configuration>& getConfigurations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addConfiguration(const Configuration &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrediction(const CVTermList &prediction)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPredictionTerm(const CVTerm &prediction)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrediction() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRetentionTime(RetentionTime rt)))
{
  // TODO
}
END_SECTION

START_SECTION((const RetentionTime& getRetentionTime() const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IncludeExcludeTarget &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IncludeExcludeTarget &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((IncludeExcludeTarget& operator=(const IncludeExcludeTarget &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([EXTRA] std::hash<IncludeExcludeTarget>))
{
  // Test that equal targets have equal hashes
  IncludeExcludeTarget t1, t2;
  t1.setName("target1");
  t1.setPrecursorMZ(500.0);
  t1.setProductMZ(200.0);
  t1.setPeptideRef("peptide_ref_1");
  t1.setCompoundRef("compound_ref_1");

  t2.setName("target1");
  t2.setPrecursorMZ(500.0);
  t2.setProductMZ(200.0);
  t2.setPeptideRef("peptide_ref_1");
  t2.setCompoundRef("compound_ref_1");

  std::hash<IncludeExcludeTarget> hasher;
  TEST_EQUAL(hasher(t1), hasher(t2))

  // Test that hash changes when values change
  IncludeExcludeTarget t3;
  t3.setName("target2");  // Different name
  t3.setPrecursorMZ(500.0);
  t3.setProductMZ(200.0);
  t3.setPeptideRef("peptide_ref_1");
  t3.setCompoundRef("compound_ref_1");
  TEST_NOT_EQUAL(hasher(t1), hasher(t3))

  // Test that different precursor_mz produces different hash
  IncludeExcludeTarget t4;
  t4.setName("target1");
  t4.setPrecursorMZ(600.0);  // Different precursor_mz
  t4.setProductMZ(200.0);
  t4.setPeptideRef("peptide_ref_1");
  t4.setCompoundRef("compound_ref_1");
  TEST_NOT_EQUAL(hasher(t1), hasher(t4))

  // Test use in unordered_set
  std::unordered_set<IncludeExcludeTarget> target_set;
  target_set.insert(t1);
  TEST_EQUAL(target_set.size(), 1)
  target_set.insert(t2); // same as t1
  TEST_EQUAL(target_set.size(), 1) // should not increase
  target_set.insert(t3);
  TEST_EQUAL(target_set.size(), 2)

  // Test use in unordered_map
  std::unordered_map<IncludeExcludeTarget, int> target_map;
  target_map[t1] = 42;
  TEST_EQUAL(target_map[t1], 42)
  TEST_EQUAL(target_map[t2], 42) // t2 == t1, should get same value
  target_map[t3] = 99;
  TEST_EQUAL(target_map[t3], 99)
  TEST_EQUAL(target_map.size(), 2)

  // Test with retention time set
  IncludeExcludeTarget t5, t6;
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(100.0);
  t5.setRetentionTime(rt);
  t6.setRetentionTime(rt);
  TEST_EQUAL(hasher(t5), hasher(t6))

  // Test with different retention times
  IncludeExcludeTarget t7;
  TargetedExperimentHelper::RetentionTime rt2;
  rt2.setRT(200.0);
  t7.setRetentionTime(rt2);
  TEST_NOT_EQUAL(hasher(t5), hasher(t7))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



