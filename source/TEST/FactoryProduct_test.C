// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff$
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

class TestProduct1 : public FactoryProduct
{
  public:
	TestProduct1(): FactoryProduct()
	{
		name_ = "TestProduct1";
		defaults_.setValue("check",0);
		defaults_.setValue("value",1);
		param_.setDefaults(defaults_,"",false);
	}
};

class TestProduct2 : public FactoryProduct
{
  public:
	TestProduct2(): FactoryProduct()
	{
		name_ = "TestProduct2";
		defaults_.setValue("value",0);
		param_.setDefaults(defaults_,"",false);
	}
};

// default ctor
FactoryProduct* ptr = 0;
CHECK(FactoryProduct())
	ptr = new FactoryProduct();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~FactoryProduct())
	delete ptr;
RESULT

// assignment operator
CHECK(FactoryProduct& operator = (const FactoryProduct& source))
	TestProduct1 fp1;
  Param p;
  p.setValue("check",1);
  fp1.setParam(p);

  TestProduct1 fp2;
  fp2 = fp1;

  TestProduct1 fp3;
  Param q;
  q.setValue("check",1);
  fp3.setParam(q);

  fp1 = TestProduct1();
	TEST_EQUAL(fp2,fp3)
RESULT

// copy constructor
CHECK(FactoryProduct(const FactoryProduct& source))
	TestProduct1 fp1;
  Param p;
  p.setValue("check",1);
  fp1.setParam(p);

  TestProduct1 fp2(fp1);

  TestProduct1 fp3;
  Param q;
  q.setValue("check",1);
  fp3.setParam(q);

  fp1 = TestProduct1();
	TEST_EQUAL(fp2, fp3)
RESULT
 
CHECK(bool operator == (const FactoryProduct& rhs) const)
	TestProduct1 s;
  Param p;
  p.setValue("check",1);
  s.setParam(p);

  TestProduct1 t;
  t.setParam(p);

  bool res = s==t;
  TEST_EQUAL(res, true)
RESULT


CHECK(bool operator != (const FactoryProduct& rhs) const)
	TestProduct1 s;
  Param p;
  p.setValue("check",1);
  s.setParam(p);

  TestProduct1 t;
  bool res = s!=t;
  TEST_EQUAL(res, true)

	TestProduct1 tp1;
  TestProduct2 tp2;
  res = tp1!=tp2;
  TEST_EQUAL(res,true)
RESULT
 
CHECK( void setParam(const Param& p))
	TestProduct1 s;
  Param p;
	p.setValue("value",1);
	p.setValue("check",0);
  TEST_EQUAL(s.getParam(),p)

  Param q;
	q.setValue("value",2);
	s.setParam(q);
	q.setValue("check",0);
  TEST_EQUAL(s.getParam(), q)
RESULT

CHECK(const Param& getParam() const)
	TestProduct1 s;
  Param p;
  p.setValue("value",1);
	p.setValue("check",0);
  TEST_EQUAL(s.getParam(), p)
RESULT

CHECK(const String& getName() const)
	TestProduct1 s;
  TEST_EQUAL(s.getName(), "TestProduct1")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
