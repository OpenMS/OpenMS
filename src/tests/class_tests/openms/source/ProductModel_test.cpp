// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>


///////////////////////////

START_TEST(ProductModel<2>, "$Id$")

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

TOLERANCE_ABSOLUTE(0.0001)

// default ctor
ProductModel* ptr = nullptr;
ProductModel* nullPointer = nullptr;
START_SECTION((ProductModel()))
	ptr = new ProductModel();
	TEST_EQUAL(ptr->getName(), "ProductModel2D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~ProductModel()))
delete ptr;
END_SECTION

START_SECTION( static const String getProductName() )
	ptr = new ProductModel();
	TEST_EQUAL(ptr->getName(), "ProductModel2D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


// assignment operator
START_SECTION((virtual ProductModel& operator=(const ProductModel &source)))
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
END_SECTION


// copy ctor
START_SECTION((ProductModel(const ProductModel& source)))
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
END_SECTION

// ModelDescription
START_SECTION((static BaseModel<D>* create()))
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
tmp1.removeAll("RT:bounding_box:");
tmp1.removeAll("RT:statistics:");
tmp1.removeAll("MZ:bounding_box:");
tmp1.removeAll("MZ:statistics:");
Param tmp2 = pm2->getParameters();
tmp2.removeAll("RT:bounding_box:");
tmp2.removeAll("RT:statistics:");
tmp2.removeAll("MZ:bounding_box:");
tmp2.removeAll("MZ:statistics:");
TEST_EQUAL(tmp1, tmp2)

DPosition<2> pos;
pos[0] = 3.5;
pos[1] = 7.5;
TEST_REAL_SIMILAR(pm3.getIntensity(pos), pm2->getIntensity(pos))
END_SECTION

START_SECTION( IntensityType getIntensity(const PositionType &pos) const )

	TOLERANCE_ABSOLUTE(0.1)	
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
	TEST_REAL_SIMILAR(pm1.getIntensity(pos), 8.52587)
	pos[0] = 2.0;
	pos[1] = 5.9;
	TEST_REAL_SIMILAR(pm1.getIntensity(pos), 0.200509)
	pos[0] = 1.8;
	pos[1] = 5.9;
	TEST_REAL_SIMILAR(pm1.getIntensity(pos), 0.0222171)
END_SECTION

START_SECTION( void getSamples(SamplesType &cont) const )
{
	GaussModel* gm1 = new GaussModel();
	gm1->setParameters(p1);
	GaussModel* gm2 = new GaussModel();
	gm2->setParameters(p2);

	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);		
		
	ProductModel pm2(pm1);

	TEST_EQUAL(pm1.getParameters(),pm2.getParameters());
	TEST_EQUAL(pm1.getModel(0)->getParameters(),pm2.getModel(0)->getParameters());
	TEST_EQUAL(pm1.getModel(1)->getParameters(),pm2.getModel(1)->getParameters());
	TEST_EQUAL(pm1.getModel(0)->getName(),pm2.getModel(0)->getName());
	TEST_EQUAL(pm1.getModel(1)->getName(),pm2.getModel(1)->getName());
		
	std::vector<Peak2D> dpa1;
	std::vector<Peak2D> dpa2;
	pm1.getSamples(dpa1);
	pm2.getSamples(dpa2);

	TEST_EQUAL(dpa1.size(),dpa2.size());
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0]);
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity());
	}
}
END_SECTION

START_SECTION( void setScale(IntensityType scale) )
		ProductModel pm1;
		pm1.setScale(3.0);
		TEST_REAL_SIMILAR(pm1.getScale(),3.0)	
END_SECTION

START_SECTION( IntensityType getScale() const )
		ProductModel pm1;
		pm1.setScale(66.6);
		TEST_REAL_SIMILAR(pm1.getScale(),66.6)	
END_SECTION

START_SECTION( ProductModel& setModel(UInt dim, BaseModel< 1 > *dist) )
	GaussModel* gm1 = new GaussModel();
	gm1->setParameters(p1);
	GaussModel* gm2 = new GaussModel();
	gm2->setParameters(p2);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	TEST_EQUAL( pm1.getModel(0) == gm1, true)
	TEST_EQUAL( pm1.getModel(1) == gm2, true)
		
END_SECTION

START_SECTION( BaseModel<1>* getModel(UInt dim) const )
	GaussModel* gm1 = new GaussModel();
	gm1->setParameters(p1);
	GaussModel* gm2 = new GaussModel();
	gm2->setParameters(p2);
	
	ProductModel pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	TEST_EQUAL( pm1.getModel(0) == gm1, true)
	TEST_EQUAL( pm1.getModel(1) == gm2, true)
		
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
