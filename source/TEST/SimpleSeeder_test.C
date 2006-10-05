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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/KERNEL/DPeakArray.h>

///////////////////////////

START_TEST(SimpleSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

enum DimensionId
{
	RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
	MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
};


// default ctor
SimpleSeeder* ptr = 0;
CHECK(Simple())
	ptr = new SimpleSeeder();
  TEST_EQUAL(ptr->getName(), "SimpleSeeder")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SimpleSeeder())
	delete ptr;
RESULT

CHECK(nextSeed())
  SimpleSeeder seeder;
  FeaFiTraits* traits = new FeaFiTraits();
  DPeakArray<2> peak_array;
  
  double mzs[] = {675, 675.5, 676, 676.5, 677};
	double rts[] = {1261, 1261, 1261, 1261, 1261};
	double its[] = {5, 10, 7, 3, 15};
	
	const Size num = 5;
	
	for (unsigned int i=0; i < num; i++)
	{
		DPeak<2> p;
		p.getPosition()[MZ] = mzs[i];
		p.getPosition()[RT] = rts[i];
		p.getIntensity()    = its[i];
		peak_array.push_back(p);
	}
	
	DPeakArray<2>::const_iterator citer1 = peak_array.begin();
	DPeakArray<2>::const_iterator citer2 = peak_array.end();
	
	traits->setData(citer1,citer2);
	
	seeder.setTraits(traits);
	
	Index peak = seeder.nextSeed();
	TEST_EQUAL(traits->getPeakIntensity(peak),15);	
	
	peak = seeder.nextSeed();
	TEST_EQUAL(traits->getPeakIntensity(peak),10);					
	
	peak = seeder.nextSeed();
	TEST_EQUAL(traits->getPeakIntensity(peak),7);		
	
	peak = seeder.nextSeed();
	TEST_EQUAL(traits->getPeakIntensity(peak),5);		
  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


