// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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
		
		void updateMembers_() override
		{
			string_var = (string)(param_.getValue("string"));
		}
		
		String string_var;
};

DefaultParamHandler* ptr = nullptr;
DefaultParamHandler* nullPointer = nullptr;
START_SECTION((DefaultParamHandler(const String& name)))
	ptr = new DefaultParamHandler("dummy");
	TEST_NOT_EQUAL(ptr, nullPointer)
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
