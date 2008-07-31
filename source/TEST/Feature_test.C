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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/Feature.h>

///////////////////////////

START_TEST(Feature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Feature* d_ptr = 0;
CHECK((Feature()))
{
	d_ptr = new Feature;
	TEST_NOT_EQUAL(d_ptr, 0);
}
RESULT

CHECK((~Feature()))
{
	delete d_ptr;
}
RESULT

CHECK((QualityType getOverallQuality() const))
	Feature p;
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
	p.setOverallQuality(123.456);
	TEST_REAL_EQUAL(p.getOverallQuality(), 123.456)
	p.setOverallQuality(-0.12345);
	TEST_REAL_EQUAL(p.getOverallQuality(), -0.12345)
	p.setOverallQuality(0.0);
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT

CHECK((void setOverallQuality(QualityType q)))
	Feature p;
	p.setOverallQuality(123.456);
	TEST_REAL_EQUAL(p.getOverallQuality(), 123.456)
	p.setOverallQuality(-0.12345);
	TEST_REAL_EQUAL(p.getOverallQuality(), -0.12345)
	p.setOverallQuality(0.0);
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT

CHECK((QualityType getQuality(UInt index) const))
	Feature p;
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	p.setQuality(0, 123.456);
	TEST_REAL_EQUAL(p.getQuality(0), 123.456)
	p.setQuality(0, -0.12345);
	TEST_REAL_EQUAL(p.getQuality(0), -0.12345)
	p.setQuality(0, 0.0);
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	TEST_REAL_EQUAL(p.getQuality(1), 0.0)
  TEST_EXCEPTION(Exception::Precondition, p.getQuality(10))
RESULT

CHECK((void setQuality(UInt index, QualityType q)))
	Feature p;
	p.setQuality(1, 123.456);
	TEST_REAL_EQUAL(p.getQuality(1), 123.456)
	p.setQuality(1, -0.12345);
	TEST_REAL_EQUAL(p.getQuality(1), -0.12345)
	p.setQuality(1, 0.0);
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	TEST_REAL_EQUAL(p.getQuality(1), 0.0)
  TEST_EXCEPTION(Exception::Precondition, p.setQuality(10,1.0))
RESULT


CHECK((const ModelDescription<2>& getModelDescription() const))
	const Feature p;
	TEST_EQUAL(p.getModelDescription().getName(), "")
	TEST_EQUAL(p.getModelDescription().getParam().empty(), true)
RESULT

CHECK((ModelDescription<2>& getModelDescription()))
	Feature p;
	TEST_EQUAL(p.getModelDescription().getName(), "")
	p.getModelDescription().setName("gauss");
	TEST_EQUAL(p.getModelDescription().getName(), "gauss")
	p.getModelDescription().setName("");
  TEST_EQUAL(p.getModelDescription().getName(), "")
RESULT

CHECK((void setModelDescription(const ModelDescription< 2 > &q)))
	Feature p;
  ModelDescription<2> desc;
  desc.setName("gauss");
	p.setModelDescription(desc);
	TEST_EQUAL(p.getModelDescription().getName(), "gauss")
	p.setModelDescription(ModelDescription<2>());
	TEST_EQUAL(p.getModelDescription().getName(), "")
RESULT

CHECK([EXTRA](IntensityType getIntensity() const))
	const Feature p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK([EXTRA](const PositionType& getPosition() const))
	const Feature	p;
	TEST_REAL_EQUAL(p.getPosition()[0], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[1], 0.0)
RESULT

CHECK([EXTRA](IntensityType& getIntensity()))
	Feature p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.setIntensity(123.456);
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.setIntensity(-0.12345);
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.setIntensity(0.0);
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK([EXTRA](PositionType& getPosition()))
	Feature::PositionType pos;
	Feature p;
	pos = p.getPosition();
	TEST_REAL_EQUAL(pos[0], 0.0)
	TEST_REAL_EQUAL(pos[1], 0.0)
	pos[0] = 1.0;
	pos[1] = 2.0;
	p.setPosition(pos);
	Feature::PositionType pos2(p.getPosition());
	TEST_REAL_EQUAL(pos2[0], 1.0)
	TEST_REAL_EQUAL(pos2[1], 2.0)
RESULT

CHECK((Feature(const Feature &feature)))
	Feature::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	Feature p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	p.setMetaValue("cluster_id",4711);
  p.setOverallQuality(0.9);
  p.setQuality(0, 0.1);
  p.setQuality(1, 0.2);
  ModelDescription<2> desc;
  desc.setName("gauss");
  p.setModelDescription(desc);

	Feature::PositionType pos2;
	Feature::IntensityType i2;

	Feature copy_of_p(p);
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();

	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)

	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4711));

  Feature::QualityType q2;
	q2 = copy_of_p.getOverallQuality();
	TEST_REAL_EQUAL(q2, 0.9)
	q2 = copy_of_p.getQuality(0);
	TEST_REAL_EQUAL(q2, 0.1)
	q2 = copy_of_p.getQuality(1);
	TEST_REAL_EQUAL(q2, 0.2)
	TEST_EQUAL(copy_of_p.getModelDescription().getName(), "gauss")
RESULT

CHECK((Feature& operator = (const Feature& rhs)))
	Feature::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	Feature p;
	p.setIntensity(123.456);
	p.setPosition(pos);
  p.setOverallQuality(0.9);
  p.setQuality(0, 0.1);
  p.setQuality(1, 0.2);
  ModelDescription<2> desc;
  desc.setName("gauss");
  p.setModelDescription(desc);
	p.setMetaValue("cluster_id",4712);

	Feature::PositionType pos2;
	Feature::IntensityType i2;

	Feature copy_of_p;
	copy_of_p = p;

	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();

  Feature::QualityType q2;

	TEST_REAL_EQUAL(i2, 123.456)
	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	q2 = copy_of_p.getOverallQuality();
	TEST_REAL_EQUAL(q2, 0.9)
	q2 = copy_of_p.getQuality(0);
	TEST_REAL_EQUAL(q2, 0.1)
	q2 = copy_of_p.getQuality(1);
	TEST_REAL_EQUAL(q2, 0.2)
	TEST_EQUAL(copy_of_p.getModelDescription().getName(), "gauss")
RESULT

CHECK((bool operator==(const Feature &rhs) const))
  ModelDescription<2> desc;
  desc.setName("gauss");

	Feature p1;
	Feature p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)

	p1.setIntensity(5);
  p1.setOverallQuality(0.9);
  p1.setQuality(0, 0.1);
  p1.setModelDescription(desc);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setIntensity(5);
  p2.setOverallQuality(0.9);
  p2.setQuality(0, 0.1);
  p2.setModelDescription(desc);
	TEST_REAL_EQUAL(p1==p2, true)

	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
RESULT

CHECK([EXTRA](Feature& operator != (const Feature& rhs)))
	Feature p1;
	Feature p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
RESULT

CHECK(([EXTRA]meta info with copy constructor))
	Feature p;
	p.setMetaValue(2,String("bla"));
 	Feature p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(([EXTRA]meta info with assignment))
	Feature p;
	p.setMetaValue(2,String("bla"));
 	Feature p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK((const std::vector<PeptideIdentification>& getPeptideIdentifications() const))
	Feature tmp;
	vector<PeptideIdentification> vec(tmp.getPeptideIdentifications());
	TEST_EQUAL(vec.size(),0);
RESULT

CHECK((void setPeptideIdentifications(const std::vector<PeptideIdentification>& identifications)))
	Feature tmp;
	vector<PeptideIdentification> vec;

	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),0);

	PeptideIdentification dbs;
	vec.push_back(dbs);
	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
RESULT

CHECK((std::vector<PeptideIdentification>& getPeptideIdentifications()))
	Feature tmp;
	vector<PeptideIdentification> vec;

	tmp.getPeptideIdentifications().resize(1);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
RESULT

//do not change these datastructures, they are used in the following tests...
std::vector< ConvexHull2D > hulls(2);
hulls[0].addPoint(DPosition<2>(1.0,2.0));
hulls[0].addPoint(DPosition<2>(3.0,4.0));
hulls[1].addPoint(DPosition<2>(0.5,0.0));
hulls[1].addPoint(DPosition<2>(1.0,1.0));

CHECK((const vector<ConvexHull2D>& getConvexHulls() const))
	Feature tmp;
	TEST_EQUAL(tmp.getConvexHulls().size(),0)
RESULT

CHECK((vector<ConvexHull2D>& getConvexHulls()))
	Feature tmp;
	tmp.setConvexHulls(hulls);
	TEST_EQUAL(tmp.getConvexHulls().size(),2)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[0][0],1.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[0][1],2.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[1][0],3.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[1][1],4.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[0][0],0.5)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[0][1],0.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[1][0],1.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[1][1],1.0)
RESULT

CHECK((void setConvexHulls(const vector<ConvexHull2D>& hulls)))
	Feature tmp;
	tmp.setConvexHulls(hulls);
	TEST_EQUAL(tmp.getConvexHulls().size(),2)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[0][0],1.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[0][1],2.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[1][0],3.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[0].getPoints()[1][1],4.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[0][0],0.5)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[0][1],0.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[1][0],1.0)
	TEST_REAL_EQUAL(tmp.getConvexHulls()[1].getPoints()[1][1],1.0)
RESULT



CHECK(ConvexHull2D& getConvexHull() const)
	Feature tmp;
	tmp.setConvexHulls(hulls);
	
	//check if the bounding box is ok
	DBoundingBox<2> bb = tmp.getConvexHull().getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],0.5)
	TEST_REAL_EQUAL(bb.min()[1],0.0)
	TEST_REAL_EQUAL(bb.max()[0],3.0)
	TEST_REAL_EQUAL(bb.max()[1],4.0)

	//check the convex hull points
	TEST_EQUAL(tmp.getConvexHull().getPoints().size(),3)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[0][0],0.5)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[0][1],0.0)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[1][0],3.0)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[1][1],4.0)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[2][0],1.0)
	TEST_REAL_EQUAL(tmp.getConvexHull().getPoints()[2][1],2.0)

RESULT

hulls[0].addPoint(DPosition<2>(3.0,2.0));
hulls[1].addPoint(DPosition<2>(2.0,1.0));
CHECK( bool encloses(DoubleReal rt, DoubleReal mz) const)
	Feature tmp;
	tmp.setConvexHulls(hulls);

	TEST_EQUAL(tmp.encloses(0.0,0.0), false);
	TEST_EQUAL(tmp.encloses(1.0,1.0), true);
	TEST_EQUAL(tmp.encloses(2.0,0.5), false);
	TEST_EQUAL(tmp.encloses(2.0,2.5), true);
	TEST_EQUAL(tmp.encloses(2.0,3.5), false);
	TEST_EQUAL(tmp.encloses(4.0,3.0), false);
	TEST_EQUAL(tmp.encloses(1.5,1.5), false);
RESULT


CHECK((const ChargeType& getCharge() const))
{
	Feature const tmp;
	TEST_EQUAL(tmp.getCharge(),0);
	// continued in setCharge()
}
RESULT

CHECK((void setCharge(const ChargeType &ch)))
{
	Feature tmp;
	TEST_EQUAL(tmp.getCharge(),0);
	tmp.setCharge(17);
	TEST_EQUAL(tmp.getCharge(),17);
}
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
