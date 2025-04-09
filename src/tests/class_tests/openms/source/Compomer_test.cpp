// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $ 
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(Compomer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Compomer* ptr = nullptr;
Compomer* nullPointer = nullptr;
START_SECTION(Compomer())
{
	ptr = new Compomer();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~Compomer())
{
	delete ptr;
}
END_SECTION

START_SECTION((Compomer(Int net_charge, double mass, double log_p)))
{
	Compomer c(34, 45.32f, 12.34f);
	TEST_EQUAL(c.getNetCharge(), 34);
	TEST_REAL_SIMILAR(c.getMass(), 45.32);	
	TEST_REAL_SIMILAR(c.getLogP(), 12.34);
}
END_SECTION


START_SECTION((Compomer(const Compomer& p) ))
{
	Compomer c(34, 45.32f, 12.34f);
	Adduct a1(123,  3, 123.456, "S", -0.3453, 0);
	Adduct b1(3, -2, 1.456, "H", -0.13, 0);
	c.setID(434);
	c.add(a1, Compomer::RIGHT);
	c.add(b1, Compomer::LEFT);

	Compomer c2(c);
	TEST_EQUAL(c2.getNetCharge(), c.getNetCharge());
	TEST_REAL_SIMILAR(c2.getMass(), c.getMass());	
	TEST_EQUAL(c2.getPositiveCharges(), c.getPositiveCharges());	
	TEST_EQUAL(c2.getNegativeCharges(), c.getNegativeCharges());
	TEST_REAL_SIMILAR(c2.getLogP(), c.getLogP());
	TEST_EQUAL(c2.getID(), c.getID());
	
}
END_SECTION

START_SECTION((Compomer& operator=(const Compomer &source)))
{
	Compomer c(34, 45.32f, 12.34f);
	Adduct a1(123,  3, 123.456, "S", -0.3453, 0);
	Adduct b1(3, -2, 1.456, "H", -0.13, 0);
	c.setID(434);
	c.add(a1, Compomer::RIGHT);
	c.add(b1, Compomer::LEFT);

	Compomer c2 = c;
	TEST_EQUAL(c2.getNetCharge(), c.getNetCharge());
	TEST_REAL_SIMILAR(c2.getMass(), c.getMass());	
	TEST_EQUAL(c2.getPositiveCharges(), c.getPositiveCharges());	
	TEST_EQUAL(c2.getNegativeCharges(), c.getNegativeCharges());
	TEST_REAL_SIMILAR(c2.getLogP(), c.getLogP());
	TEST_EQUAL(c2.getID(), c.getID());
}
END_SECTION


START_SECTION([EXTRA] friend OPENMS_DLLAPI bool operator==(const Compomer& a, const  Compomer& b))
{
	Compomer c(34, 45.32f, 12.34f);
	Adduct a1(123,  3, 123.456, "S", -0.3453, 0);
	Adduct b1(3, -2, 1.456, "H", -0.13, 0);
	c.setID(434);
	c.add(a1, Compomer::RIGHT);
	
	Compomer c2(c);
	TEST_TRUE(c == c2);
	c.setID(2);
	TEST_EQUAL(c==c2, false);
	
}
END_SECTION

START_SECTION((void add(const Adduct &a, UInt side)))
{
	//Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob
	Adduct a1(123, 43, 123.456, "S", -0.3453, 0);
	Adduct a2(123,  3, 123.456, "S", -0.3453, 0);

	Adduct b1(3, -2, 1.456, "H", -0.13, 0);

	Compomer c;
	c.add(a1, Compomer::RIGHT);
	//PRECISION(0.0001);
	TEST_EQUAL(c.getNetCharge(), 123*43);
	TEST_REAL_SIMILAR(c.getMass(), 123.456*43);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*43);
	TEST_EQUAL(c.getPositiveCharges(), 123*43);
	TEST_EQUAL(c.getNegativeCharges(), 0);
	
	c.add(a2, Compomer::RIGHT);
	TEST_EQUAL(c.getNetCharge(), 123*46);
	TEST_REAL_SIMILAR(c.getMass(), 123.456*46);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*46);
	TEST_EQUAL(c.getPositiveCharges(), 123*46);
	TEST_EQUAL(c.getNegativeCharges(), 0);
	
	c.add(b1, Compomer::RIGHT);
	TEST_EQUAL(c.getNetCharge(), 123*46+ 3*(-2));
	TEST_REAL_SIMILAR(c.getMass(), 123.456*46 - 2*1.456);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*46 -0.13*2);
	TEST_EQUAL(c.getPositiveCharges(), 123*46);
	TEST_EQUAL(c.getNegativeCharges(), 6);	
	
		
}
END_SECTION

START_SECTION(bool isConflicting(const Compomer &cmp, UInt side_this, UInt side_other) const)
{
	EmpiricalFormula ef("H");
	Adduct default_adduct(1, 1, ef.getMonoWeight(), ef.toString(), log(0.7f), 0);

	{
	Adduct a1(1, 1, 1.007, "H1", -0.13, 0);
	Adduct a2(1, 2, 123.456, "NH4", -0.3453, 0);

	
	Compomer c,d;
	c.add(a1, Compomer::RIGHT);
	d.add(a1, Compomer::RIGHT);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), true);

	// this should not change the result	
	c.add(a1, Compomer::RIGHT);
	d.add(a1, Compomer::RIGHT);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), true);

	// this neither
	c.add(a2, Compomer::LEFT);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), true);
	}
	
	{
  Adduct a1(1, -2, 123.456f, "NH4", -0.3453f, 0);
	Adduct a2(1, 1, 1.007, "H1", -0.13f, 0);
	Adduct b1(1, 2, 1.007, "H1", -0.13, 0);

	Compomer c,d;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	d.add(b1, Compomer::RIGHT);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::LEFT), false);
	}

	{
  Adduct a1(1, 3, 123.456, "NH4", -0.3453, 0);
	Adduct a2(1, 3, 1.007, "H1", -0.13, 0);

	Compomer c,d;
	c.add(a1, Compomer::RIGHT);
	d.add(a1, Compomer::LEFT);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::LEFT), true);

	c.add(a1, Compomer::LEFT);
	c.add(a2, Compomer::RIGHT);
	d.add(a1, Compomer::LEFT);
	d.add(a2, Compomer::RIGHT);
	//		C										D
	//a1				a1a2	; 	a1a1	a2
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::LEFT), true);

	c.add(a1, Compomer::RIGHT);
	d.add(a2, Compomer::LEFT);
	
	d.add(a1, Compomer::RIGHT);
	d.add(a1, Compomer::RIGHT);
	//		C										D
	//a1				a1a2a1	; 	a1a1a2	a2a1a1
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::LEFT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::RIGHT,Compomer::RIGHT), false);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::RIGHT), true);
	TEST_EQUAL(c.isConflicting(d,Compomer::LEFT,Compomer::LEFT), true);

	}

	
}
END_SECTION

START_SECTION((void setID(const Size &id)))
{
  NOT_TESTABLE //well.. tested below...	
}
END_SECTION

START_SECTION((const Size& getID() const))
{
  Compomer c;
	c.setID(123);
	TEST_EQUAL(c.getID(), 123)
}
END_SECTION

START_SECTION((const Int& getNetCharge() const))
{
  Compomer c(-123,1.23,-0.12);
	TEST_EQUAL(c.getNetCharge(), -123)
}
END_SECTION

START_SECTION((const double& getMass() const))
{
  Compomer c(1,-123.12, 0.23);
	TEST_REAL_SIMILAR(c.getMass(), -123.12)
}
END_SECTION


START_SECTION((const Int& getPositiveCharges() const))
{
  Compomer c;
	Adduct a1(3, -2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(6, 1, 1.007, "H1", -0.13, 0);

	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	TEST_EQUAL(c.getPositiveCharges(), 6)
}
END_SECTION

START_SECTION((const Int& getNegativeCharges() const))
{
  Compomer c;
	Adduct a1(3, -2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(6, 1, 1.007, "H1", -0.13, 0);

	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	TEST_EQUAL(c.getNegativeCharges(), 6)
}
END_SECTION


START_SECTION((const double& getLogP() const))
{
  Compomer c(1,1,-123.12);
	TEST_REAL_SIMILAR(c.getLogP(), -123.12)
}
END_SECTION

START_SECTION((const double& getRTShift() const))
{
  Compomer c(1,1,-123.12);
  Adduct a(123, 43, 123.456, "S", -0.3453, -10.12);  
  c.add(a,0);
	TEST_REAL_SIMILAR(c.getRTShift(), 435.16)
}
END_SECTION

START_SECTION((StringList getLabels(const UInt side) const))
{
  Compomer c(1,1,-123.12);
	TEST_EQUAL(c.getLabels(0).size(), 0)
  Adduct a(123, 43, 123.456, "S", -0.3453, -10.12, "testlabel");  
  c.add(a,0);
	TEST_EQUAL(c.getLabels(0).size(), 1)
	TEST_EQUAL(c.getLabels(1).size(), 0)
}
END_SECTION

START_SECTION((String getAdductsAsString() const))
{
  Adduct a1(1, 2, 123.456f, "NH4", -0.3453f, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	TEST_EQUAL(c.getAdductsAsString(), "() --> (H-1H8N2)");
	c.add(a1, Compomer::LEFT);
	TEST_EQUAL(c.getAdductsAsString(), "(H8N2) --> (H-1H8N2)");
}
END_SECTION

START_SECTION((String getAdductsAsString(UInt side) const))
{
  Adduct a1(1, 2, 123.456f, "NH4", -0.3453f, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	TEST_EQUAL(c.getAdductsAsString(Compomer::LEFT), "");
	TEST_EQUAL(c.getAdductsAsString(Compomer::RIGHT), "H-1H8N2");
	c.add(a1, Compomer::LEFT);
	TEST_EQUAL(c.getAdductsAsString(Compomer::LEFT), "H8N2");
	TEST_EQUAL(c.getAdductsAsString(Compomer::RIGHT), "H-1H8N2");
}
END_SECTION

START_SECTION((const CompomerComponents& getComponent() const))
{
  Adduct a1(1, 2, 123.456f, "NH4", -0.3453f, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	Compomer::CompomerComponents comp(2);
	TEST_EQUAL(c.getComponent()==comp, true);

	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	c.add(a1, Compomer::LEFT);
	comp[Compomer::RIGHT][a1.getFormula()] = a1;
	comp[Compomer::RIGHT][a2.getFormula()] = a2;
	comp[Compomer::LEFT][a1.getFormula()] = a1;
	TEST_EQUAL(c.getComponent()==comp, true);
}
END_SECTION

START_SECTION((Compomer removeAdduct(const Adduct& a) const))
{
  Adduct a1(1, 2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	c.add(a1, Compomer::LEFT);
	Compomer tmp = c.removeAdduct(a1);
	TEST_EQUAL(tmp.getAdductsAsString(), "() --> (H-1)");
}
END_SECTION

START_SECTION((Compomer removeAdduct(const Adduct& a, const UInt side) const))
{
  Adduct a1(1, 2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	c.add(a1, Compomer::LEFT);
	Compomer tmp = c.removeAdduct(a1, Compomer::RIGHT);
	TEST_EQUAL(tmp.getAdductsAsString(), "(H8N2) --> (H-1)");
					 tmp = c.removeAdduct(a1, Compomer::LEFT);
	TEST_EQUAL(tmp.getAdductsAsString(), "() --> (H-1H8N2)");	
}
END_SECTION

START_SECTION((void add(const CompomerSide& add_side, UInt side)))
{
  Adduct a1(1, 2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	c.add(a1, Compomer::LEFT);
	TEST_EQUAL(c.getAdductsAsString(), "(H8N2) --> (H-1H8N2)");
	Compomer tmp = c;
	tmp.add(c.getComponent()[Compomer::RIGHT], Compomer::RIGHT);
	TEST_EQUAL(tmp.getAdductsAsString(), "(H8N2) --> (H-2H16N4)");
	tmp.add(c.getComponent()[Compomer::RIGHT], Compomer::LEFT);
	TEST_EQUAL(tmp.getAdductsAsString(), "(H-1H16N4) --> (H-2H16N4)");	
}
END_SECTION

START_SECTION((bool isSingleAdduct(Adduct &a, const UInt side) const))
  Adduct a1(1, 2, 123.456, "NH4", -0.3453, 0);
	Adduct a2(1, -1, 1.007, "H1", -0.13, 0);
	Compomer c;
	c.add(a1, Compomer::RIGHT);
	c.add(a2, Compomer::RIGHT);
	c.add(a1, Compomer::LEFT);
	TEST_EQUAL(c.isSingleAdduct(a1,Compomer::LEFT), true);
	TEST_EQUAL(c.isSingleAdduct(a2,Compomer::LEFT), false);
	TEST_EQUAL(c.isSingleAdduct(a1,Compomer::RIGHT), false);
	TEST_EQUAL(c.isSingleAdduct(a2,Compomer::RIGHT), false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



