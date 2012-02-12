// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Chris Bielow  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/SingletonRegistry.h>
#include <OpenMS/CONCEPT/Factory.h>


#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(<SingletonRegistry>, "$Id$")
/////////////////////////////////////////////////////////////
FactoryBase* nullPointer = 0;

START_SECTION(static FactoryBase* getFactory(const String& name))
	Factory<FilterFunctor>::create("TICFilter");
	String myName = typeid(Factory<FilterFunctor>).name();


  TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), nullPointer)
END_SECTION


START_SECTION(static void registerFactory(const String& name, FactoryBase* instance))
	String myName = typeid(FactoryBase).name();
	FactoryBase* fb = new FactoryBase;
	SingletonRegistry::registerFactory(myName, fb);
  TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), nullPointer)
END_SECTION

START_SECTION(static bool isRegistered(String name))
	TEST_EQUAL(SingletonRegistry::isRegistered(typeid(Factory<FilterFunctor>).name()), true)
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


