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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>

///////////////////////////

START_TEST(<Factory>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


// Factory is singleton, therefore we don't test the constructor
CHECK(static FactoryProduct* create(const String& name))
	DataReducer* p = Factory<DataReducer>::create("max_reducer");
	MaxReducer reducer;
	TEST_EQUAL(*p,reducer);
RESULT

CHECK( static Factory* instance() )
	Factory<DataReducer>* factory = Factory<DataReducer>::instance();
	TEST_NOT_EQUAL(factory, 0)
RESULT

CHECK( static void registerProduct(const String& name, const FunctionType creator) )
	Factory<DataReducer>::registerProduct(MaxReducer::getProductName(), &MaxReducer::create);
	DataReducer* ext = Factory<DataReducer>::create("max_reducer");
	TEST_NOT_EQUAL(ext, 0)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


