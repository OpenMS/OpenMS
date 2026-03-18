// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(ModelDescription<2>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

ModelDescription<2>* ptr = nullptr;
ModelDescription<2>* nullPointer = nullptr;
START_SECTION((ModelDescription()))
	ptr = new ModelDescription<2>();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ModelDescription()))
	delete ptr;
END_SECTION

START_SECTION( virtual bool operator==(const ModelDescription &rhs) const )
	ModelDescription<2> fp1,fp2;	
	TEST_EQUAL(fp1==fp2,true)
  
  fp1.setName("halligalli2000");
	TEST_EQUAL(fp1==fp2,false)
  
  fp2.setName("halligalli2000");
  TEST_EQUAL(fp1==fp2,true)
  
  Param param;
  param.setValue("bla","bluff");
	fp1.setParam(param);
	TEST_EQUAL(fp1==fp2,false)
  
	fp2.setParam(param);
  TEST_EQUAL(fp1==fp2,true)
END_SECTION

START_SECTION( virtual bool operator!=(const ModelDescription &rhs) const )
	ModelDescription<2> fp1, fp2;
	TEST_EQUAL(fp1!=fp2,false)
	
  fp1.setName("halligalli2000");
	TEST_EQUAL(fp1!=fp2,true)

  fp2.setName("halligalli2000");
	TEST_EQUAL(fp1!=fp2,false)
END_SECTION

START_SECTION((virtual ModelDescription& operator=(const ModelDescription &source)))
	ModelDescription<2> tm1;
  tm1.setName("halligalli");
	Param param;
	param.setValue("test","test");
	tm1.setParam(param);
	
  ModelDescription<2> tm2;
  tm2 = tm1;

	TEST_EQUAL(tm1==tm2,true)
END_SECTION

START_SECTION((ModelDescription(const ModelDescription &source)))
	ModelDescription<2> tm1;
  tm1.setName("halligalli");
	Param param;
	param.setValue("test","test");
	tm1.setParam(param);
	
  ModelDescription<2> tm2(tm1);

	TEST_EQUAL(tm1==tm2,true)
END_SECTION

START_SECTION((const String& getName() const ))
  const ModelDescription<2> m;
  TEST_EQUAL(m.getName(), "")
END_SECTION

START_SECTION((void setName(const String &name)))
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
END_SECTION

START_SECTION((const Param& getParam() const ))
  ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

START_SECTION( String& getName() )
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
END_SECTION

START_SECTION( Param& getParam() )
	ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

START_SECTION( void setParam(const Param &param) )
	ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
