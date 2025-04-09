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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



