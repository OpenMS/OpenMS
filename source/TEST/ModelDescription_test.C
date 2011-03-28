// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(ModelDescription<2>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

ModelDescription<2>* ptr = 0;
ModelDescription<2>* nullPointer = 0;
START_SECTION((ModelDescription()))
	ptr = new ModelDescription<2>();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ModelDescription()))
	delete ptr;
END_SECTION

START_SECTION(BaseModel<D>* createModel())
	BaseModel<2>* ptr = ModelDescription<2>().createModel();
	TEST_EQUAL(ptr, 0)	// no name is set, should be zero pointer
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

START_SECTION( ModelDescription(const BaseModel< D > *model) )
	const BaseModel<1> * bm = new IsotopeModel();

  ModelDescription<1> md(bm);
	
	BaseModel<1>* ptr = md.createModel();
	TEST_EQUAL( *ptr == *bm, true)	
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
