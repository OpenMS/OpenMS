// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(FactoryProduct, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

class TestProduct1 
	: public FactoryProduct
{
  public:
	TestProduct1()
		: FactoryProduct("TestProduct1")
	{
		defaults_.setValue("check",0,"desc1");
		defaults_.setValue("value",1,"desc2");
		
		defaultsToParam_();
	}
	
	TestProduct1(const TestProduct1& rhs)
		: FactoryProduct(rhs)
	{
		updateMembers_();
	}

	TestProduct1& operator=(const TestProduct1& rhs)
	{
		if (&rhs==this) return *this;
		
		FactoryProduct::operator=(rhs);
		updateMembers_();
		
		return *this;
	}
	
	void updateMembers_()
	{
		check = param_.getValue("check");
	}
	
	int check;
};

FactoryProduct* ptr = 0;
CHECK((FactoryProduct(const String& name)))
	ptr = new FactoryProduct("TEST");
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK([EXTRA] ~FactoryProduct())
	delete ptr;
RESULT

CHECK([EXTRA] const String& getName() const)
	TestProduct1 s;
  TEST_EQUAL(s.getName(), "TestProduct1")
RESULT

CHECK([EXTRA] const Param& getParameters() const)
	TestProduct1 s;
  Param p;
  p.setValue("value",1);
	p.setValue("check",0);
  TEST_EQUAL(s.getParameters(), p)
RESULT


CHECK([EXTRA] void setParameters(const Param& p))
	TestProduct1 s;
  Param p;
	p.setValue("value",1);
	p.setValue("check",0);
  TEST_EQUAL(s.getParameters(),p)

  Param q;
	q.setValue("value",2);
	s.setParameters(q);
	q.setValue("check",0);
  TEST_EQUAL(s.getParameters(), q)
RESULT

CHECK((FactoryProduct& operator = (const FactoryProduct& source)))
	TestProduct1 fp1;
  Param p;
  p.setValue("check",1);
  fp1.setParameters(p);

  TestProduct1 fp2;
  fp2 = fp1;
	
	TEST_EQUAL(fp1,fp2)
RESULT

CHECK((FactoryProduct(const FactoryProduct& source)))
	TestProduct1 fp1;
  Param p;
  p.setValue("check",1);
  fp1.setParameters(p);

  TestProduct1 fp2(fp1);

	TEST_EQUAL(fp1, fp2)
RESULT
 
CHECK((bool operator == (const FactoryProduct& rhs) const))
	TestProduct1 s,t;
  Param p;
  p.setValue("check",1);

  TEST_EQUAL(s==t, true)
  
  s.setParameters(p);

  TEST_EQUAL(s==t, false)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
