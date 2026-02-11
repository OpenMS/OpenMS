// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>

#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;
using namespace std;

START_TEST(ChargePair, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChargePair* ptr = nullptr;
ChargePair* nullPointer = nullptr;
START_SECTION(ChargePair())
{
	ptr = new ChargePair();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ChargePair())
{
	delete ptr;
}
END_SECTION

Compomer cmp;
cmp.setID(99);

START_SECTION((ChargePair(const Size &index0, const Size &index1, const Int &charge0, const Int &charge1, const Compomer &compomer, const double &mass_diff, const bool active)))
{
	ChargePair cp(34,45, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp (cp2);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 1);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair& operator=(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp = cp2;
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 1);
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION


START_SECTION((Int getCharge(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...	
}
END_SECTION

START_SECTION((void setCharge(UInt pairID, Int e)))
{
  ChargePair cp;
	cp.setCharge(0,123);
	cp.setCharge(1,321);
	TEST_EQUAL(cp.getCharge(0), 123)
	TEST_EQUAL(cp.getCharge(1), 321)
}
END_SECTION

START_SECTION((Size getElementIndex(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setElementIndex(UInt pairID, Size e)))
{
  ChargePair cp;
	cp.setElementIndex(0,123);
	cp.setElementIndex(1,321);
	TEST_EQUAL(cp.getElementIndex(0), 123)
	TEST_EQUAL(cp.getElementIndex(1), 321)
}
END_SECTION

START_SECTION((const Compomer& getCompomer() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setCompomer(const Compomer &compomer)))
{
  ChargePair cp;
	cp.setCompomer(cmp);
	TEST_EQUAL(cp.getCompomer(), cmp)
}
END_SECTION

START_SECTION((double getMassDiff() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setMassDiff(double mass_diff)))
{
  ChargePair cp;
	cp.setMassDiff(123.432);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 123.432)
}
END_SECTION

START_SECTION((double getEdgeScore() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setEdgeScore(double score)))
{
  ChargePair cp;
	cp.setEdgeScore(1123.432f);
	TEST_REAL_SIMILAR(cp.getEdgeScore(), 1123.432)
}
END_SECTION
		

START_SECTION((bool isActive() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setActive(const bool active)))
{
  ChargePair cp;
	cp.setActive(true);
	TEST_EQUAL(cp.isActive(), true)
	cp.setActive(false);
	TEST_EQUAL(cp.isActive(), false)
}
END_SECTION

START_SECTION((virtual bool operator==(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp2(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp1==cp2, false);
	ChargePair cp3(34,15, 4,5, cmp, 12.34, true);
	ChargePair cp4(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp3==cp4, false);
	ChargePair cp5(34,15, 4,5, cmp, 12.34, false);
	ChargePair cp6(34,15, 4,5, cmp, 12.34, false);
	TEST_TRUE(cp5 == cp6);
	
}
END_SECTION

START_SECTION((virtual bool operator!=(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp2(34,15, 4,5, cmp, 12.34, false);
	TEST_FALSE(cp1 == cp2);
	ChargePair cp3(34,15, 4,5, cmp, 12.34, true);
	ChargePair cp4(34,15, 4,5, cmp, 12.34, false);
	TEST_FALSE(cp3 == cp4);
	ChargePair cp5(34,15, 4,5, cmp, 12.34, false);
	ChargePair cp6(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp5!=cp6, false);
}
END_SECTION

START_SECTION(([EXTRA] std::hash<ChargePair>))
{
  std::hash<ChargePair> hasher;

  // Test that equal objects have equal hashes
  ChargePair cp1(34,45, 4,5, cmp, 12.34, false);
  ChargePair cp2(34,45, 4,5, cmp, 12.34, false);
  TEST_TRUE(cp1 == cp2);
  TEST_EQUAL(hasher(cp1), hasher(cp2));

  // Test that different objects (likely) have different hashes
  ChargePair cp3(34,15, 4,5, cmp, 12.34, false);
  TEST_FALSE(cp1 == cp3);
  TEST_NOT_EQUAL(hasher(cp1), hasher(cp3));

  // Test with different charges
  ChargePair cp4(34,45, 3,5, cmp, 12.34, false);
  TEST_FALSE(cp1 == cp4);
  TEST_NOT_EQUAL(hasher(cp1), hasher(cp4));

  // Test with different active state
  ChargePair cp5(34,45, 4,5, cmp, 12.34, true);
  TEST_FALSE(cp1 == cp5);
  TEST_NOT_EQUAL(hasher(cp1), hasher(cp5));

  // Test with different mass_diff
  ChargePair cp6(34,45, 4,5, cmp, 99.99, false);
  TEST_FALSE(cp1 == cp6);
  TEST_NOT_EQUAL(hasher(cp1), hasher(cp6));

  // Test that score_ (not in operator==) does not affect hash
  ChargePair cp7(34,45, 4,5, cmp, 12.34, false);
  cp7.setEdgeScore(999.0);
  ChargePair cp8(34,45, 4,5, cmp, 12.34, false);
  cp8.setEdgeScore(1.0);
  TEST_TRUE(cp7 == cp8);
  TEST_EQUAL(hasher(cp7), hasher(cp8));

  // Test use in unordered_set
  std::unordered_set<ChargePair> cp_set;
  cp_set.insert(cp1);
  cp_set.insert(cp2); // duplicate, should not increase size
  cp_set.insert(cp3);
  TEST_EQUAL(cp_set.size(), 2);
  TEST_EQUAL(cp_set.count(cp1), 1);
  TEST_EQUAL(cp_set.count(cp3), 1);

  // Test use in unordered_map
  std::unordered_map<ChargePair, int> cp_map;
  cp_map[cp1] = 100;
  cp_map[cp3] = 200;
  TEST_EQUAL(cp_map.size(), 2);
  TEST_EQUAL(cp_map[cp1], 100);
  TEST_EQUAL(cp_map[cp3], 200);
  cp_map[cp2] = 150; // cp2 == cp1, should overwrite
  TEST_EQUAL(cp_map.size(), 2);
  TEST_EQUAL(cp_map[cp1], 150);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



