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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/SingletonRegistry.h>
#include <OpenMS/CONCEPT/Factory.h>


#include <OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>

///////////////////////////

START_TEST(<SingletonRegistry>, "$Id:$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


CHECK(static FactoryBase* getFactory(const String& name))
	Factory<DataReducer>::create("max_reducer");
	String myName = typeid(Factory<DataReducer>).name();

	TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), 0)
RESULT


CHECK(static void registerFactory(const String& name, FactoryBase* chosenOne))
	String myName = typeid(FactoryBase).name();
	FactoryBase* fb = new FactoryBase;
	SingletonRegistry::registerFactory(myName, fb);
	TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), 0)
RESULT

CHECK(static bool isRegistered(String name))
	TEST_EQUAL(SingletonRegistry::isRegistered(typeid(Factory<DataReducer>).name()), true)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


