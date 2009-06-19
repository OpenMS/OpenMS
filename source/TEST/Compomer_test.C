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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(Compomer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Compomer* ptr = 0;
START_SECTION(Compomer())
{
	ptr = new Compomer();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~Compomer())
{
	delete ptr;
}
END_SECTION

START_SECTION((Compomer(Int net_charge, DoubleReal mass, DoubleReal log_p)))
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
	Adduct a1(123,  3, 123.456f, "SECRET", -0.3453f);
	Adduct b1(3, -2, 1.456f, "H", -0.13f);
	c.setID(434);
	c.add(a1);
	c.add(b1);

	Compomer c2(c);
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
	Adduct a1(123,  3, 123.456f, "SECRET", -0.3453f);
	Adduct b1(3, -2, 1.456f, "H", -0.13f);
	c.setID(434);
	c.add(a1);
	
	Compomer c2(c);
	TEST_EQUAL(c==c2, true);
	c.setID(2);
	TEST_EQUAL(c==c2, false);
	
}
END_SECTION

START_SECTION((void add(const Adduct &a)))
{
	//Adduct(Int charge, Int amount, DoubleReal singleMass, String formula, DoubleReal log_prob
	Adduct a1(123, 43, 123.456f, "SECRET", -0.3453f);
	Adduct a2(123,  3, 123.456f, "SECRET", -0.3453f);

	Adduct b1(3, -2, 1.456f, "H", -0.13f);

	Compomer c;
	c.add(a1);
	//PRECISION(0.0001);
	TEST_EQUAL(c.getNetCharge(), 123*43);
	TEST_REAL_SIMILAR(c.getMass(), 123.456*43);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*43);
	TEST_EQUAL(c.getPositiveCharges(), 123*43);
	TEST_EQUAL(c.getNegativeCharges(), 0);
	
	c.add(a2);
	TEST_EQUAL(c.getNetCharge(), 123*46);
	TEST_REAL_SIMILAR(c.getMass(), 123.456*46);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*46);
	TEST_EQUAL(c.getPositiveCharges(), 123*46);
	TEST_EQUAL(c.getNegativeCharges(), 0);
	
	c.add(b1);
	TEST_EQUAL(c.getNetCharge(), 123*46+ 3*(-2));
	TEST_REAL_SIMILAR(c.getMass(), 123.456*46 - 2*1.456);	
	TEST_REAL_SIMILAR(c.getLogP(), -0.3453*46 -0.13*2);
	TEST_EQUAL(c.getPositiveCharges(), 123*46);
	TEST_EQUAL(c.getNegativeCharges(), 6);	
	
		
}
END_SECTION

START_SECTION((bool isConflicting(const Compomer &cmp, bool left_this, bool left_other) const))
{
	EmpiricalFormula ef("H");
	Adduct default_adduct(1, 1, ef.getMonoWeight(), ef.getString(), log(0.7f));

	{
	Adduct a1(1, 2, 123.456f, "NH4", -0.3453f);
	Adduct a2(1, -1, 1.007f, "H1", -0.13f);

	Adduct b1(1, -1, 1.007f, "H1", -0.13f);

	Compomer c,d;
	c.add(a1);
	c.add(a2);
	d.add(b1);
	TEST_EQUAL(c.isConflicting(d,true,true,default_adduct*6,default_adduct*6), false);
	TEST_EQUAL(c.isConflicting(d,true,true,default_adduct*2,default_adduct), true);

	TEST_EQUAL(c.isConflicting(d,true,false,default_adduct,default_adduct), true);
	TEST_EQUAL(c.isConflicting(d,false,true,default_adduct,default_adduct), true);
	TEST_EQUAL(c.isConflicting(d,false,false,default_adduct,default_adduct), true);
	}
	
	{
  Adduct a1(1, -2, 123.456f, "NH4", -0.3453f);
	Adduct a2(1, 1, 1.007f, "H1", -0.13f);

	Adduct b1(1, 2, 1.007f, "H1", -0.13f);

	Compomer c,d;
	c.add(a1);
	c.add(a2);
	d.add(b1);
	TEST_EQUAL(c.isConflicting(d,true,true,default_adduct*5,default_adduct*4), true);
	TEST_EQUAL(c.isConflicting(d,true,false,default_adduct,default_adduct), true);
	TEST_EQUAL(c.isConflicting(d,false,true,default_adduct,default_adduct), true);
	TEST_EQUAL(c.isConflicting(d,false,false,default_adduct,default_adduct), true);
	TEST_EQUAL(c.isConflicting(d,false,false,default_adduct*5,default_adduct*4), false);

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

START_SECTION((const DoubleReal& getMass() const))
{
  Compomer c(1,-123.12, 0.23);
	TEST_REAL_SIMILAR(c.getMass(), -123.12)
}
END_SECTION


START_SECTION((const Int& getPositiveCharges() const))
{
  Compomer c;
	Adduct a1(3, -2, 123.456f, "NH4", -0.3453f);
	Adduct a2(6, 1, 1.007f, "H1", -0.13f);

	c.add(a1);
	c.add(a2);
	TEST_EQUAL(c.getPositiveCharges(), 6)
}
END_SECTION

START_SECTION((const Int& getNegativeCharges() const))
{
  Compomer c;
	Adduct a1(3, -2, 123.456f, "NH4", -0.3453f);
	Adduct a2(6, 1, 1.007f, "H1", -0.13f);

	c.add(a1);
	c.add(a2);
	TEST_EQUAL(c.getNegativeCharges(), 6)
}
END_SECTION


START_SECTION((const DoubleReal& getLogP() const))
{
  Compomer c(1,1,-123.12);
	TEST_REAL_SIMILAR(c.getLogP(), -123.12)
}
END_SECTION

START_SECTION((String getAdductsAsString(Int side=0)))
{
  Adduct a1(1, 2, 123.456f, "NH4", -0.3453f);
	Adduct a2(1, -1, 1.007f, "H1", -0.13f);
	Compomer c;
	c.add(a1);
	c.add(a2);
	TEST_EQUAL(c.getAdductsAsString(), "-1(H1)2(NH4)");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



