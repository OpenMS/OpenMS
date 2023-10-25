// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>

#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(<Factory>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Factory is singleton, therefore we don't test the constructor
START_SECTION(static FactoryProduct* create(const String& name))
	FilterFunctor* p = Factory<FilterFunctor>::create("TICFilter");
	TICFilter reducer;
	TEST_EQUAL(*p==reducer,true);
	delete p;
END_SECTION

START_SECTION( static void registerProduct(const String& name, const FunctionType creator) )
	Factory<FilterFunctor>::registerProduct(TICFilter::getProductName(), &TICFilter::create);
	FilterFunctor* ext = Factory<FilterFunctor>::create("TICFilter");
  FilterFunctor* nullPointer = nullptr;
  TEST_NOT_EQUAL(ext, nullPointer)
  delete ext;
END_SECTION

START_SECTION(static bool isRegistered(const String& name))
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter"), true)
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter_bla_bluff"), false)
END_SECTION

START_SECTION(static std::vector<String> registeredProducts())
	vector<String> list = Factory<FilterFunctor>::registeredProducts();
	TEST_EQUAL(list.size(),6)
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{

   int nr_iterations (1e2);
   int test = 0;
#pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    FilterFunctor* p = Factory<FilterFunctor>::create("TICFilter");
    TICFilter reducer;
    test += (*p == reducer);
    delete p;
  }
  TEST_EQUAL(test, nr_iterations)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


