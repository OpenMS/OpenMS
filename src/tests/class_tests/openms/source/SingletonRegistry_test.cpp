// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/SingletonRegistry.h>
#include <OpenMS/CONCEPT/Factory.h>


#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(<SingletonRegistry>, "$Id$")
/////////////////////////////////////////////////////////////
FactoryBase* nullPointer = nullptr;

START_SECTION(static FactoryBase* getFactory(const String& name))
	auto ptr = Factory<FilterFunctor>::create("TICFilter");
	String myName = typeid(Factory<FilterFunctor>).name();

  TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), nullPointer)
  delete ptr;
END_SECTION


START_SECTION(static void registerFactory(const String& name, FactoryBase* instance))
	String myName = typeid(FactoryBase).name();
	FactoryBase* fb = new FactoryBase;
	SingletonRegistry::registerFactory(myName, fb);
  TEST_NOT_EQUAL(SingletonRegistry::getFactory(myName), nullPointer)
  delete fb;
END_SECTION

START_SECTION(static bool isRegistered(String name))
	TEST_EQUAL(SingletonRegistry::isRegistered(typeid(Factory<FilterFunctor>).name()), true)
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


