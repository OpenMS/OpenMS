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

CHECK((virtual ~EuclideanDistance()))
	delete ptr;
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<1>& model, UInt dim)))
	
	EuclideanDistance dist;
	GaussModel gm1;
	gm1.setScalingFactor(5.0);
	gm1.setInterpolationStep(0.3);
	Param tmp;
	tmp.setValue("bounding_box:min",1);
	tmp.setValue("bounding_box:max",5 );
	tmp.setValue("statistics:variance",3.0 );
	tmp.setValue("statistics:mean",2.5 );			
	gm1.setParameters(tmp);
	
	FeaFiTraits traits;
	DPeakArray<Peak2D> peak_array;
	
	Peak2D p1;
	p1.getPosition()[0] = 1;
	p1.getPosition()[1] = 1;
	p1.setIntensity(0);
	peak_array.push_back(p1);
	
	Peak2D p2;
	p2.getPosition()[0] = 2;
	p2.getPosition()[1] = 2;
	p2.setIntensity(3);
	peak_array.push_back(p2);
	
	Peak2D p3;
	p3.getPosition()[0] = 3;
	p3.getPosition()[1] = 3;
	p3.setIntensity(5);
	peak_array.push_back(p3);
	
	Peak2D p4;
	p4.getPosition()[0] = 4;
	p4.getPosition()[1] = 4;
	p4.setIntensity(3);
	peak_array.push_back(p4);
	
	Peak2D p5;
	p5.getPosition()[0] = 5;
	p5.getPosition()[1] = 5;
	p5.setIntensity(0);
	peak_array.push_back(p5);
	
	MSExperimentExtern<Peak1D > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);
		
	dist.setTraits(&traits);
	
	FeaFiModule::IndexSet set;
	for (UInt i=0; i<=4; ++i)
	{
		set.insert(std::make_pair(i,0));
	}
	
	// evaluate rt dimension
	double result = dist.evaluate(set, gm1,0);
	TEST_REAL_EQUAL(result,-6.10346)
	// evaluate mz dimension
	result = dist.evaluate(set, gm1,1);
	TEST_REAL_EQUAL(result,-6.10346);
	
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<2>& model)))

	EuclideanDistance dist;
	GaussModel * gm1 = new GaussModel();
	GaussModel * gm2 = new GaussModel();
			
	gm1->setScalingFactor(5.0);
	gm1->setInterpolationStep(0.3);
	Param tmp;
	tmp.setValue("bounding_box:min",1);
	tmp.setValue("bounding_box:max",5 );
	tmp.setValue("statistics:variance",3.0);
	tmp.setValue("statistics:mean",2.5 );			
	gm1->setParameters(tmp);
	
	gm2->setScalingFactor(5.0);
	gm2->setInterpolationStep(0.3);
	gm2->setParameters(tmp);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	FeaFiTraits traits;
	DPeakArray<Peak2D> peak_array;
	
	Peak2D p1;
	p1.getPosition()[0] = 1;
	p1.getPosition()[1] = 1;
	p1.setIntensity(0);
	peak_array.push_back(p1);
	
	Peak2D p2;
	p2.getPosition()[0] = 2;
	p2.getPosition()[1] = 2;
	p2.setIntensity(3);
	peak_array.push_back(p2);
	
	Peak2D p3;
	p3.getPosition()[0] = 3;
	p3.getPosition()[1] = 3;
	p3.setIntensity(5);
	peak_array.push_back(p3);
	
	Peak2D p4;
	p4.getPosition()[0] = 4;
	p4.getPosition()[1] = 4;
	p4.setIntensity(3);
	peak_array.push_back(p4);
	
	Peak2D p5;
	p5.getPosition()[0] = 5;
	p5.getPosition()[1] = 5;
	p5.setIntensity(0);
	peak_array.push_back(p5);
	
	MSExperimentExtern<Peak1D > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);
	
	dist.setTraits(&traits);
	
	FeaFiModule::IndexSet  set;
	for (UInt i=0; i<=4; ++i)
	{
		set.insert(std::make_pair(i,0));
	}
	
	double result = dist.evaluate(set, pm1);
	double pval   = dist.getPvalue();
	TEST_REAL_EQUAL(result,-6.42946);
	TEST_REAL_EQUAL(pval,-1);	// euclidean distance does not have a p-value
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(EuclideanDistance::getProductName(),"EuclideanDistance")
	TEST_EQUAL(EuclideanDistance().getName(),"EuclideanDistance")
RESULT

CHECK((static BaseQuality* create()))
	TEST_NOT_EQUAL(EuclideanDistance::create(),0)
RESULT

CHECK( double getPvalue() )
	TEST_EQUAL( EuclideanDistance().getPvalue() == (-1),true)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



