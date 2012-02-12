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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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
END_SECTION

START_SECTION( static void registerProduct(const String& name, const FunctionType creator) )
	Factory<FilterFunctor>::registerProduct(TICFilter::getProductName(), &TICFilter::create);
	FilterFunctor* ext = Factory<FilterFunctor>::create("TICFilter");
  FilterFunctor* nullPointer = 0;
  TEST_NOT_EQUAL(ext, nullPointer)
END_SECTION

START_SECTION(static bool isRegistered(const String& name))
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter"), true)
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter_bla_bluff"), false)
END_SECTION

START_SECTION(static std::vector<String> registeredProducts())
	vector<String> list = Factory<FilterFunctor>::registeredProducts();
	TEST_EQUAL(list.size(),6)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


