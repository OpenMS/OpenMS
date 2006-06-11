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
// $Id: SimpleExtender_test.C,v 1.6 2006/04/08 18:24:13 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/IndexSet.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/DPeak.h>

#include <OpenMS/FORMAT/Param.h>

///////////////////////////

START_TEST(SimpleExtender, "$Id: SimpleExtender_test.C,v 1.6 2006/04/08 18:24:13 ole_st Exp $")

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
SimpleExtender* ptr = 0;
CHECK(Simple())
	ptr = new SimpleExtender();
  TEST_EQUAL(ptr->getName(), "SimpleExtender")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SimpleExtender())
	delete ptr;
RESULT

CHECK(nextSeed())
  SimpleExtender extender;
  FeaFiTraits* traits = new FeaFiTraits();
  DPeakArrayNonPolymorphic<2> peak_array;
  
  double mzs[] = {675, 675.5, 676, 676.5, 677};
	double rts[] = {1260, 1260, 1260, 1260, 1260};
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
	
	DPeakArrayNonPolymorphic<2>::const_iterator citer1 = peak_array.begin();
	DPeakArrayNonPolymorphic<2>::const_iterator citer2 = peak_array.end();
	
	traits->setData(citer1,citer2);
	
	extender.setTraits(traits);

	Param param;
	param.setValue("dist_mz_up",10);
	param.setValue("dist_mz_down",10);
	param.setValue("dist_rt_down",10);
	param.setValue("dist_rt_up",10);
	param.setValue("priority_thr",-5);
	
  extender.setParam(param);
  IndexSet region = extender.extend(4);
  
  TEST_NOT_EQUAL(region.size(),0);
  
  IndexSet::const_iterator ind_iter = region.begin();
  
    
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


