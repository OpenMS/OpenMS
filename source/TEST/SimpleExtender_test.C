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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/DPeak.h>

#include <OpenMS/FORMAT/Param.h>

///////////////////////////

START_TEST(SimpleExtender, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

enum DimensionId
{
	RT = DimensionDescription < LCMS_Tag >::RT,
	MZ = DimensionDescription < LCMS_Tag >::MZ
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
  MSExperiment<DPeak<1> > exp;
	MSExperiment<DPeak<1> >::SpectrumType spec;
	spec.setRetentionTime(1260);
  
  double mzs[] = {675, 675.5, 676, 676.5, 677};
	double its[] = {5, 10, 7, 3, 15};
	
	const Size num = 5;
	
	for (unsigned int i=0; i < num; i++)
	{
		DPeak<1> p;
		p.getPosition()[0] = mzs[i];
		p.getIntensity() = its[i];
		
		spec.push_back(p);
	}
	
	exp.push_back(spec);
	
	traits->setData(exp.begin(), exp.end(),100);
	
	extender.setTraits(traits);
 	
	FeaFiModule::IndexSet  set;
	for (UnsignedInt i=0; i< 5; ++i)
	{
		set.insert(std::make_pair(0,i));
	}
  FeaFiModule::IndexSet  region = extender.extend(set);
  
  TEST_NOT_EQUAL(region.size(),0);
  
    
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


