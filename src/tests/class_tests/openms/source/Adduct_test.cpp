// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Adduct.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Adduct, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Adduct* ptr = nullptr;
Adduct* nullPointer = nullptr;
START_SECTION(Adduct())
{
	ptr = new Adduct();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~Adduct())
{
	delete ptr;
}
END_SECTION

START_SECTION((Adduct(Int charge)))
{
	Adduct a(123);
	TEST_EQUAL(a.getCharge(), 123);
}
END_SECTION

START_SECTION((Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, const String label="")))
{
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
  TEST_REAL_SIMILAR(a.getRTShift(), -10);
  TEST_EQUAL(a.getLabel(), "");

	Adduct a2(123, 43, 123.456f, "S", -0.3453, -10, "testlabel");
  TEST_EQUAL(a2.getLabel(), "testlabel");	
}
END_SECTION

START_SECTION([EXTRA] friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b))
{
	
	Adduct a(123,  3, 123.456f, "S", -0.3453f, 0);
	Adduct b(a);

	TEST_TRUE(a == b);
	a.setAmount(22);
	TEST_EQUAL(a==b, false);
	
}
END_SECTION

START_SECTION((const Int& getCharge() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setCharge(const Int &charge)))
{
	Adduct a;
	a.setCharge(123);
	TEST_EQUAL(a.getCharge(), 123);
}
END_SECTION

START_SECTION((const Int& getAmount() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setAmount(const Int &amount)))
{
	Adduct a;
	a.setAmount(43);
  TEST_EQUAL(a.getAmount(), 43);
}
END_SECTION

START_SECTION((const double& getSingleMass() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setSingleMass(const double &singleMass)))
{
	Adduct a;
	a.setSingleMass(43.21);
  TEST_REAL_SIMILAR(a.getSingleMass(), 43.21);
}
END_SECTION

START_SECTION((const double& getLogProb() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setLogProb(const double &log_prob)))
{
	Adduct a;
	a.setLogProb(43.21f);
  TEST_REAL_SIMILAR(a.getLogProb(), 43.21);
}
END_SECTION

START_SECTION((const String& getFormula() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setFormula(const String &formula)))
	Adduct a;
	a.setFormula("S");
  TEST_EQUAL(a.getFormula()=="S1", true);
END_SECTION

START_SECTION((const double& getRTShift() const))
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
  TEST_REAL_SIMILAR(a.getRTShift(), -10);
	Adduct a1(123, 43, 123.456f, "S", -0.3453, 11);
  TEST_REAL_SIMILAR(a1.getRTShift(), 11);
END_SECTION

START_SECTION((const String& getLabel() const ))
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
  TEST_EQUAL(a.getLabel(), "");
	Adduct a1(123, 43, 123.456f, "S", -0.3453, 11, "mylabel");
  TEST_EQUAL(a1.getLabel(), "mylabel");
END_SECTION



START_SECTION((Adduct operator *(const Int m) const))
{
	Adduct a_p(123, 43, 123.456, "S", -0.3453, 0);
	Adduct a = a_p*4;
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43*4);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456f);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
}
END_SECTION

START_SECTION((Adduct operator+(const Adduct &rhs)))
{
	Adduct a_p(123, 43, 123.456f, "S", -0.3453f, 0);
	Adduct a_p2(123, 40, 123.456f, "S", -0.3453f, 0);
	Adduct a = a_p + a_p2;
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43+40);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
}
END_SECTION

START_SECTION((void operator+=(const Adduct &rhs)))
{
	Adduct a_p(123, 43, 123.456f, "S", -0.3453f, 0);
	Adduct a(a_p);
	a.setAmount(10);
	a	+= a_p;
	TEST_EQUAL(a.getAmount(), 43+10);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



