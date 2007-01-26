// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/FactoryProduct2.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(FactoryProduct2, "$Id: FactoryProduct2_test.C 1371 2007-01-19 17:18:08Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

class TestProduct1 : public FactoryProduct2
{
  public:
		TestProduct1()
			: FactoryProduct2("TestProduct1")
		{
			check_defaults_ = false;
			
			defaults_.setValue("int",1);
			
			defaultsToParam_();
		}
	
};

FactoryProduct2* ptr = 0;
CHECK(FactoryProduct2(const String& name))
	ptr = new TestProduct1();
	TEST_NOT_EQUAL(ptr, 0)
	TEST_EQUAL(ptr->getName(),"TestProduct1")
RESULT

CHECK(~FactoryProduct2())
	delete ptr;
RESULT

Param p;
p.setValue("int",1);
p.setValue("string","bla");

CHECK((FactoryProduct2& operator = (const FactoryProduct2& source)))
	TestProduct1 fp1;
  fp1.setParameters(p);

  TestProduct1 fp2;
  fp2 = fp1;
	TEST_EQUAL(fp1.getParameters(),fp2.getParameters())
	
	fp2 = TestProduct1();
	TEST_EQUAL(fp2.getParameters().size(),1)
RESULT

CHECK((FactoryProduct2(const FactoryProduct2& source)))
	TestProduct1 fp1,fp4;
  fp1.setParameters(p);

  TestProduct1 fp2(fp1);
	TEST_EQUAL(fp1.getParameters(),fp2.getParameters())
	
  TestProduct1 fp3(fp4);
	TEST_EQUAL(fp3.getParameters().size(),1)
RESULT
 
CHECK((bool operator == (const FactoryProduct2& rhs) const))
	TestProduct1 s;
  TestProduct1 t;
  
  TEST_EQUAL(s==t, true)
  
  s.setParameters(p);

  TEST_EQUAL(s==t, false)
  
  s = t;

  TEST_EQUAL(s==t, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
