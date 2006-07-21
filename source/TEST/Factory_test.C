// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>

///////////////////////////

START_TEST(<Factory>, "$Id: Factory_test.C,v 1.4 2006/03/28 12:53:13 ole_st Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


// Factory is singleton, therefore we don't test the constructor
CHECK(Factory<DummyProduct>::create(name))
	BaseSeeder* p = Factory<BaseSeeder>::create("SimpleSeeder");
	SimpleSeeder seeder;
	TEST_EQUAL(*p,seeder);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


