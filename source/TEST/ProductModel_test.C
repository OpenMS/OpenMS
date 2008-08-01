// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>


///////////////////////////

START_TEST(ProductModel<3>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

typedef ProductModel<2> ProductModel;

Param p1;
p1.setValue("bounding_box:min",1.0f);
p1.setValue("bounding_box:max",4.0f);
p1.setValue("statistics:mean",3.0f);
p1.setValue("statistics:variance",0.1f);

Param p2;
p2.setValue("bounding_box:min",5.0f);
p2.setValue("bounding_box:max",6.0f);
p2.setValue("statistics:mean",7.0f);
p2.setValue("statistics:variance",0.3f);

PRECISION(0.0001)

// default ctor
ProductModel* ptr = 0;
CHECK((ProductModel()))
	ptr = new ProductModel();
	TEST_EQUAL(ptr->getName(), "ProductModel2D")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~ProductModel()))
delete ptr;
RESULT

CHECK( static const String getProductName() )
	ptr = new ProductModel();
	TEST_EQUAL(ptr->getName(), "ProductModel2D")
	TEST_NOT_EQUAL(ptr, 0)
RESULT


// assignment operator
CHECK((virtual ProductModel& operator=(const ProductModel &source)))
GaussModel* gm1 = new GaussModel();
gm1->setParameters(p1);
GaussModel* gm2 = new GaussModel();
gm2->setParameters(p2);
GaussModel* gm3 = new GaussModel();
gm3->setParameters(p1);
GaussModel* gm4 = new GaussModel();
gm4->setParameters(p2);

ProductModel pm1;
pm1.setModel(0,gm1);
pm1.setModel(1,gm2);

ProductModel pm2;
pm2 = pm1;

ProductModel pm3;
pm3.setModel(0,gm3);
pm3.setModel(1,gm4);

pm1 = ProductModel();

TEST_EQUAL(pm2.getParameters(), pm3.getParameters())
RESULT


// copy ctor
CHECK((ProductModel(const ProductModel& source)))
GaussModel* gm1 = new GaussModel();
gm1->setParameters(p1);
GaussModel* gm2 = new GaussModel();
gm2->setParameters(p2);
GaussModel* gm3 = new GaussModel();
gm3->setParameters(p1);
GaussModel* gm4 = new GaussModel();
gm4->setParameters(p2);

ProductModel pm1;
pm1.setModel(0,gm1);
pm1.setModel(1,gm2);
ProductModel pm2(pm1);

ProductModel pm3;
pm3.setModel(0,gm3);
pm3.setModel(1,gm4);

pm1 = ProductModel();
TEST_EQUAL(pm3.getParameters(), pm2.getParameters())
RESULT

// ModelDescription
CHECK((static BaseModel<D>* create()))
GaussModel* gm1 = new GaussModel();
GaussModel* gm2 = new GaussModel();
GaussModel* gm3 = new GaussModel();
gm3->setParameters(p1);
GaussModel* gm4 = new GaussModel();
gm4->setParameters(p2);

ProductModel pm1;
pm1.setModel(0,gm1);
pm1.setModel(1,gm2);
pm1.setScale(4.0);
pm1.setCutOff(0.5);
gm1->setParameters(p1);
gm2->setParameters(p2);

ModelDescription<2> md(&pm1);
ProductModel* pm2 = static_cast< ProductModel* >(md.createModel());

ProductModel pm3;
pm3.setModel(0,gm3);
pm3.setModel(1,gm4);
pm3.setScale(4.0);
pm3.setCutOff(0.5);

pm1 = ProductModel();

//remove fitting data and compare
Param tmp1 = pm3.getParameters();
tmp1.remove("RT:bounding_box:");
tmp1.remove("RT:statistics:");
tmp1.remove("MZ:bounding_box:");
tmp1.remove("MZ:statistics:");
Param tmp2 = pm2->getParameters();
tmp2.remove("RT:bounding_box:");
tmp2.remove("RT:statistics:");
tmp2.remove("MZ:bounding_box:");
tmp2.remove("MZ:statistics:");
TEST_EQUAL(tmp1, tmp2)

DPosition<2> pos;
pos[0] = 3.5;
pos[1] = 7.5;
TEST_REAL_EQUAL(pm3.getIntensity(pos), pm2->getIntensity(pos))
RESULT

CHECK( IntensityType getIntensity(const PositionType &pos) const )

	PRECISION(0.1)	
	GaussModel* gm1 = new GaussModel();
	GaussModel* gm2 = new GaussModel();
	gm1->setParameters(p1);
	gm2->setParameters(p2);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	pm1.setScale(10.0);
	pm1.setCutOff(0.01);
	
	DPosition<2> pos;
	pos[0] = 2.5;
	pos[1] = 5.9;
	TEST_REAL_EQUAL(pm1.getIntensity(pos), 8.52587)
	pos[0] = 2.0;
	pos[1] = 5.9;
	TEST_REAL_EQUAL(pm1.getIntensity(pos), 0.200509)
	pos[0] = 1.8;
	pos[1] = 5.9;
	TEST_REAL_EQUAL(pm1.getIntensity(pos), 0.0222171)
RESULT

CHECK( void getSamples(SamplesType &cont) const )
		
		GaussModel* gm1 = new GaussModel();
		gm1->setParameters(p1);
		GaussModel* gm2 = new GaussModel();
		gm2->setParameters(p2);

		ProductModel pm1;
		pm1.setModel(0,gm1);
		pm1.setModel(1,gm2);		
		
		ProductModel pm2(pm1);

		TEST_EQUAL(pm1.getParameters(),pm2.getParameters())
		TEST_EQUAL(pm1.getModel(0)->getParameters(),pm2.getModel(0)->getParameters())
		TEST_EQUAL(pm1.getModel(1)->getParameters(),pm2.getModel(1)->getParameters())		
		TEST_EQUAL(pm1.getModel(0)->getName(),pm2.getModel(0)->getName())
		TEST_EQUAL(pm1.getModel(1)->getName(),pm2.getModel(1)->getName())	
		
		DPeakArray<DPeak<2> > dpa1;
		DPeakArray<DPeak<2> > dpa2;
		pm1.getSamples(dpa1);
		pm2.getSamples(dpa2);

		TEST_EQUAL(dpa1.size(),dpa2.size())
		ABORT_IF(dpa1.size()!=dpa2.size());
		for (UInt i=0; i<dpa1.size(); ++i)
		{
			TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
			TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
		}
RESULT

CHECK( void setScale(IntensityType scale) )
		ProductModel pm1;
		pm1.setScale(3.0);
		TEST_REAL_EQUAL(pm1.getScale(),3.0)	
RESULT

CHECK( IntensityType getScale() const )
		ProductModel pm1;
		pm1.setScale(66.6);
		TEST_REAL_EQUAL(pm1.getScale(),66.6)	
RESULT

CHECK( ProductModel& setModel(UInt dim, BaseModel< 1 > *dist) )
	GaussModel* gm1 = new GaussModel();
	gm1->setParameters(p1);
	GaussModel* gm2 = new GaussModel();
	gm2->setParameters(p2);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	TEST_EQUAL( pm1.getModel(0) == gm1, true)
	TEST_EQUAL( pm1.getModel(1) == gm2, true)
		
RESULT

CHECK( BaseModel<1>* getModel(UInt dim) const )
	GaussModel* gm1 = new GaussModel();
	gm1->setParameters(p1);
	GaussModel* gm2 = new GaussModel();
	gm2->setParameters(p2);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	TEST_EQUAL( pm1.getModel(0) == gm1, true)
	TEST_EQUAL( pm1.getModel(1) == gm2, true)
		
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
