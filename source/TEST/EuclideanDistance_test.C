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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EuclideanDistance.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

///////////////////////////

START_TEST(EuclideanDistance, "$Id EuclideanDistance_test.C 139 2006-07-14 10:08:39Z ole_st $")

using namespace OpenMS;
using namespace std;
using namespace OpenMS::Math;

typedef ProductModel<2> ProductModel;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EuclideanDistance* ptr = 0;
CHECK((EuclideanDistance()))
	ptr = new EuclideanDistance();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~EuclideanDistance()))
	delete ptr;
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)))
	
	EuclideanDistance dist;
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(2.5);
	stat.setVariance(3.0);
	gm1.setScalingFactor(5.0);
	gm1.setInterpolationStep(0.3);
	gm1.setParam(stat,1,5);
	
	FeaFiTraits traits;
	DPeakArray<2> peak_array;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1;
	p1.getPosition()[1] = 1;
	p1.getIntensity()    = 0;
	peak_array.push_back(p1);
	
	DPeak<2> p2;
	p2.getPosition()[0] = 2;
	p2.getPosition()[1] = 2;
	p2.getIntensity()    = 3;
	peak_array.push_back(p2);
	
	DPeak<2> p3;
	p3.getPosition()[0] = 3;
	p3.getPosition()[1] = 3;
	p3.getIntensity()    = 5;
	peak_array.push_back(p3);
	
	DPeak<2> p4;
	p4.getPosition()[0] = 4;
	p4.getPosition()[1] = 4;
	p4.getIntensity()    = 3;
	peak_array.push_back(p4);
	
	DPeak<2> p5;
	p5.getPosition()[0] = 5;
	p5.getPosition()[1] = 5;
	p5.getIntensity()    = 0;
	peak_array.push_back(p5);
	
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);
		
	dist.setTraits(&traits);
	
	FeaFiModule::IndexSet set;
	for (UnsignedInt i=0; i<=4; ++i)
	{
		set.insert(std::make_pair(i,0));
	}
	
	// evaluate rt dimension
	double result = dist.evaluate(set, gm1,0);
	TEST_REAL_EQUAL(result,-4.54373)
	// evaluate mz dimension
	result = dist.evaluate(set, gm1,1);
	TEST_REAL_EQUAL(result,-4.54373)
	
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<2>& model)))

	EuclideanDistance dist;
	GaussModel * gm1 = new GaussModel();
	GaussModel * gm2 = new GaussModel();
		
	BasicStatistics<>  stat;
	stat.setMean(2.5);
	stat.setVariance(3.0);
	
	gm1->setScalingFactor(5.0);
	gm1->setInterpolationStep(0.3);
	gm1->setParam(stat,1,5);
	
	gm2->setScalingFactor(5.0);
	gm2->setInterpolationStep(0.3);
	gm2->setParam(stat,1,5);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	FeaFiTraits traits;
	DPeakArray<2> peak_array;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1;
	p1.getPosition()[1] = 1;
	p1.getIntensity()    = 0;
	peak_array.push_back(p1);
	
	DPeak<2> p2;
	p2.getPosition()[0] = 2;
	p2.getPosition()[1] = 2;
	p2.getIntensity()    = 3;
	peak_array.push_back(p2);
	
	DPeak<2> p3;
	p3.getPosition()[0] = 3;
	p3.getPosition()[1] = 3;
	p3.getIntensity()    = 5;
	peak_array.push_back(p3);
	
	DPeak<2> p4;
	p4.getPosition()[0] = 4;
	p4.getPosition()[1] = 4;
	p4.getIntensity()    = 3;
	peak_array.push_back(p4);
	
	DPeak<2> p5;
	p5.getPosition()[0] = 5;
	p5.getPosition()[1] = 5;
	p5.getIntensity()    = 0;
	peak_array.push_back(p5);
	
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);
	
	dist.setTraits(&traits);
	
	FeaFiModule::IndexSet  set;
	for (UnsignedInt i=0; i<=4; ++i)
	{
		set.insert(std::make_pair(i,0));
	}
	
	double result = dist.evaluate(set, pm1);
	double pval   = dist.getPvalue();
	TEST_REAL_EQUAL(result,-3.88326);
	TEST_REAL_EQUAL(pval,-1);	// euclidean distance does not have a p-value
RESULT

CHECK(const String getName())
	TEST_EQUAL(EuclideanDistance::getProductName(),"EuclideanDistance")
	TEST_EQUAL(EuclideanDistance().getName(),"EuclideanDistance")
RESULT

CHECK(static BaseQuality* create())
	TEST_NOT_EQUAL(EuclideanDistance::create(),0)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



