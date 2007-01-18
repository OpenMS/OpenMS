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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(ModelDescription<2>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;


// default ctor
ModelDescription<2>* ptr = 0;
CHECK(ModelDescription<2>())
	ptr = new ModelDescription<2>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(ModelDescription<2>())
	delete ptr;
RESULT

// assignment operator
CHECK(ModelDescription<2>& operator = (const ModelDescription<2>& source))
	ModelDescription<2> tm1;
  tm1.setName("halligalli");
  ModelDescription<2> tm2;
  tm2 = tm1;

  ModelDescription<2> tm3;
	tm3.setName("halligalli");
  
// 	tm1 = ModelDescription<2>();
	TEST_EQUAL(tm3==tm2,true)
RESULT

// copy constructor
CHECK(ModelDescription<2>(const ModelDescription<2>& source))
	ModelDescription<2> fp1;	
  fp1.setName("halligalli2000");

  ModelDescription<2> fp2(fp1);

  ModelDescription<2> fp3;
  fp3.setName("halligalli2000");

//   fp1 = ModelDescription<2>();
	TEST_EQUAL(fp2==fp3,true)
RESULT

CHECK(String getName() const)
  const ModelDescription<2> m;
  TEST_EQUAL(m.getName(), "")
RESULT

CHECK(String setName() )
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
RESULT

CHECK(String getParam() const)
  ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
RESULT

CHECK(String setName() )
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
RESULT 


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
