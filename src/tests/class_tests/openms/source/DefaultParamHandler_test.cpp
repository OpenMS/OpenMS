// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

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
    
    explicit TestHandler(const String& name)
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

START_SECTION(static void writeParametersToMetaValues(const Param& write_this, MetaInfoInterface& write_here, const String& prefix = ""))
  MetaInfoInterface meta_values;
  Param p;
  p.setValue("int", 1);
  p.setValue("string", "test");
  p.setValue("ignore:bli", 4711);
  DefaultParamHandler::writeParametersToMetaValues(p, meta_values);
  DefaultParamHandler::writeParametersToMetaValues(p, meta_values, "prefix");
  ABORT_IF(!meta_values.metaValueExists("int"))
  ABORT_IF(!meta_values.metaValueExists("string"))
  ABORT_IF(!meta_values.metaValueExists("bli"))
  TEST_EQUAL(meta_values.getMetaValue("int"), 1)
  TEST_EQUAL(meta_values.getMetaValue("string"), "test")
  TEST_EQUAL(meta_values.getMetaValue("bli"), 4711)
  ABORT_IF(!meta_values.metaValueExists("prefix:int"))
  ABORT_IF(!meta_values.metaValueExists("prefix:string"))
  ABORT_IF(!meta_values.metaValueExists("prefix:bli"))
  TEST_EQUAL(meta_values.getMetaValue("prefix:int"), 1)
  TEST_EQUAL(meta_values.getMetaValue("prefix:string"), "test")
  TEST_EQUAL(meta_values.getMetaValue("prefix:bli"), 4711)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
