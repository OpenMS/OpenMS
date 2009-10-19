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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

///////////////////////////

START_TEST(DefaultParamHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

class TestHandler 
	: public DefaultParamHandler
{
  public:
  	
		TestHandler(const String& name)
			: DefaultParamHandler(name)
		{
			defaults_.setValue("int",0,"intdesc");
			defaults_.setValue("string","default","stingdesc");
			subsections_.push_back("ignore");
			
			defaultsToParam_();
		}
		
		TestHandler(const TestHandler& rhs)
			: DefaultParamHandler(rhs),
				string_var(rhs.string_var)
		{
			updateMembers_();
		}
		
		TestHandler& operator= (const TestHandler& rhs)
		{
			if (&rhs == this) return *this;
			
			DefaultParamHandler::operator=(rhs);
			string_var = rhs.string_var;
			
			updateMembers_();
			
			return *this;
		}
		
		void updateMembers_()
		{
			string_var = (string)(param_.getValue("string"));
		}
		
		String string_var;
};

DefaultParamHandler* ptr = 0;
START_SECTION((DefaultParamHandler(const String& name)))
	ptr = new DefaultParamHandler("dummy");
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~DefaultParamHandler()))
	delete ptr;
END_SECTION

START_SECTION((const String& getName() const))
	DefaultParamHandler s("dummy2");
  TEST_EQUAL(s.getName(), "dummy2")
END_SECTION

START_SECTION((void setName(const String& name)))
	DefaultParamHandler s("dummy2");
	s.setName("SetName");
  TEST_EQUAL(s.getName(), "SetName")
END_SECTION

START_SECTION((const std::vector<String>& getSubsections() const))
	DefaultParamHandler s("dummy2");
  TEST_EQUAL(s.getSubsections().size(),0)
END_SECTION

START_SECTION((const Param& getDefaults() const))
	DefaultParamHandler s("dummy2");
  TEST_EQUAL(s.getDefaults().size(),0)
	TestHandler t("dummy2");
  TEST_EQUAL(t.getDefaults().size(),2)
END_SECTION

START_SECTION((const Param& getParameters() const))
	TestHandler s("dummy");
	Param empty;
  TEST_EQUAL(s.getParameters().size(),2)
  TEST_EQUAL((int)(s.getParameters().getValue("int")),0)
  TEST_EQUAL((string)(s.getParameters().getValue("string")),"default")
  TEST_EQUAL(s.string_var, "default")
END_SECTION

START_SECTION((void setParameters(const Param &param)))
	Param p;
	p.setValue("int",1);
	p.setValue("string","test");
	p.setValue("ignore:bli",4711);
	
	TestHandler s("dummy");
  s.setParameters(p);
	
  TEST_EQUAL((int)(s.getParameters().getValue("int")), 1)
	TEST_EQUAL((string)(s.getParameters().getValue("string")), "test")
	TEST_EQUAL(s.string_var, "test")
END_SECTION

START_SECTION((bool operator == (const DefaultParamHandler& rhs) const))
	TestHandler empty("dummy");
  TestHandler h("dummy");
  TEST_EQUAL(empty==h,true);
  
  Param p;
	p.setValue("int",1);
  h.setParameters(p);
	TEST_EQUAL(empty==h,false);
END_SECTION

START_SECTION((DefaultParamHandler & operator=(const DefaultParamHandler &rhs)))
	Param p;
	p.setValue("int",1);
	p.setValue("string","test");
	p.setValue("ignore:bli",4711);
	
	TestHandler s("dummy");
  s.setParameters(p);
	
	TestHandler s2("dummy2");
	s2 = s;
  TEST_EQUAL((int)(s2.getParameters().getValue("int")), 1)
	TEST_EQUAL((string)(s2.getParameters().getValue("string")), "test")
	TEST_EQUAL(s2.string_var, "test")
	
	s2 = TestHandler("dummy");
	TEST_EQUAL(s2==TestHandler("dummy"),true)
END_SECTION

START_SECTION((DefaultParamHandler(const DefaultParamHandler &rhs)))
	Param p;
	p.setValue("int",1);
	p.setValue("string","test");
	p.setValue("ignore:bli",4711);
	
	TestHandler s("dummy");
  s.setParameters(p);
	
	TestHandler s2(s);
	
  TEST_EQUAL((int)(s2.getParameters().getValue("int")), 1)
	TEST_EQUAL((string)(s2.getParameters().getValue("string")), "test")
	TEST_EQUAL(s2.string_var, "test")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
