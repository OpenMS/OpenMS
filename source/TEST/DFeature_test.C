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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DFeature.h>

///////////////////////////

START_TEST(DFeature<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DFeature<10>* d10_ptr = 0;
CHECK(DFeature())
	d10_ptr = new DFeature<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~DFeature())
	delete d10_ptr;
RESULT

CHECK(const QualityType& getOverallQuality() const)
	const DFeature<10> p;
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT

CHECK(QualityType& getOverallQuality())
	DFeature<3> p;
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
	p.getOverallQuality() = 123.456;
	TEST_REAL_EQUAL(p.getOverallQuality(), 123.456)
	p.getOverallQuality() = -0.12345;
	TEST_REAL_EQUAL(p.getOverallQuality(), -0.12345)
	p.getOverallQuality() = 0.0;
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT

CHECK(setOverallQuality(QualityType))
	DFeature<3> p;
	p.setOverallQuality(123.456);
	TEST_REAL_EQUAL(p.getOverallQuality(), 123.456)
	p.setOverallQuality(-0.12345);
	TEST_REAL_EQUAL(p.getOverallQuality(), -0.12345)
	p.setOverallQuality(0.0);
	TEST_REAL_EQUAL(p.getOverallQuality(), 0.0)
RESULT


CHECK(const QualityType& getQuality(Position) const)
	const DFeature<10> p;
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	TEST_REAL_EQUAL(p.getQuality(1), 0.0)
	TEST_REAL_EQUAL(p.getQuality(2), 0.0)
	TEST_REAL_EQUAL(p.getQuality(3), 0.0)
	TEST_REAL_EQUAL(p.getQuality(4), 0.0)
	TEST_REAL_EQUAL(p.getQuality(5), 0.0)
	TEST_REAL_EQUAL(p.getQuality(6), 0.0)
	TEST_REAL_EQUAL(p.getQuality(7), 0.0)
	TEST_REAL_EQUAL(p.getQuality(8), 0.0)
	TEST_REAL_EQUAL(p.getQuality(9), 0.0)
#ifdef OPENMS_DEBUG
  TEST_EXCEPTION(Exception::Precondition, p.getQuality(10))
#endif
RESULT

CHECK(QualityType& getQuality(Position))
	DFeature<3> p;
	TEST_REAL_EQUAL(p.getQuality(2), 0.0)
	p.getQuality(2) = 123.456;
	TEST_REAL_EQUAL(p.getQuality(2), 123.456)
	p.getQuality(2) = -0.12345;
	TEST_REAL_EQUAL(p.getQuality(2), -0.12345)
	p.getQuality(2) = 0.0;
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	TEST_REAL_EQUAL(p.getQuality(1), 0.0)
	TEST_REAL_EQUAL(p.getQuality(2), 0.0)
#ifdef OPENMS_DEBUG
  TEST_EXCEPTION(Exception::Precondition, p.getQuality(10))
#endif
RESULT

CHECK(setQuality(QualityType,Position))
	DFeature<3> p;
	p.setQuality(2, 123.456);
	TEST_REAL_EQUAL(p.getQuality(2), 123.456)
	p.setQuality(2, -0.12345);
	TEST_REAL_EQUAL(p.getQuality(2), -0.12345)
	p.setQuality(2, 0.0);
	TEST_REAL_EQUAL(p.getQuality(0), 0.0)
	TEST_REAL_EQUAL(p.getQuality(1), 0.0)
	TEST_REAL_EQUAL(p.getQuality(2), 0.0)
#ifdef OPENMS_DEBUG
  TEST_EXCEPTION(Exception::Precondition, p.setQuality(10,1.0))
#endif
RESULT


CHECK(const QualityType& getModelDescription() const)
	const DFeature<10> p;
	TEST_EQUAL(p.getModelDescription().getName(), "")
	TEST_EQUAL(p.getModelDescription().getParam().empty(), true)
RESULT

CHECK(QualityType& getModelDescription())
	DFeature<3> p;
	TEST_EQUAL(p.getModelDescription().getName(), "")
	p.getModelDescription().getName() = "gauss";
	TEST_EQUAL(p.getModelDescription().getName(), "gauss")
	p.getModelDescription().getName() = "";
  TEST_EQUAL(p.getModelDescription().getName(), "")
RESULT

CHECK(setModelDescription(const ModelDescription&))
	DFeature<3> p;
  ModelDescription<3> desc;
  desc.setName("gauss");
	p.setModelDescription(desc);
	TEST_EQUAL(p.getModelDescription().getName(), "gauss")
	p.setModelDescription(ModelDescription<3>());
	TEST_EQUAL(p.getModelDescription().getName(), "")
RESULT

CHECK(const IntensityType& getIntensity() const)
	const DFeature<10> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(const PositionType& getPosition() const)
	const DFeature<10>	p;
	TEST_REAL_EQUAL(p.getPosition()[0], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[1], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[2], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[3], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[4], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[5], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[6], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[7], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[8], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[9], 0.0)
RESULT

CHECK(IntensityType& getIntensity())
	DFeature<3> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.getIntensity() = 123.456;
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.getIntensity() = -0.12345;
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.getIntensity() = 0.0;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(PositionType& getPosition())
	DFeature<3>::PositionType pos;
	DFeature<3> p;
	pos = p.getPosition();
	TEST_REAL_EQUAL(pos[0], 0.0)
	TEST_REAL_EQUAL(pos[1], 0.0)
	TEST_REAL_EQUAL(pos[2], 0.0)
	pos[0] = 1.0;
	pos[1] = 2.0;
	pos[2] = 3.0;
	p.getPosition() = pos;
	DFeature<3>::PositionType pos2(p.getPosition());
	TEST_REAL_EQUAL(pos2[0], 1.0)
	TEST_REAL_EQUAL(pos2[1], 2.0)
	TEST_REAL_EQUAL(pos2[2], 3.0)	
RESULT

CHECK(DFeature(const DFeature<D>& p))
	DFeature<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DFeature<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
	p.setMetaValue("cluster_id",4711);
  p.getOverallQuality() = 0.9;
  p.setQuality(0, 0.1);
  p.setQuality(1, 0.2);
  p.setQuality(2, 0.3);
  ModelDescription<3> desc;
  desc.setName("gauss");
  p.setModelDescription(desc);
	
	DFeature<3>::PositionType pos2;
	DFeature<3>::IntensityType i2;

	DFeature<3> copy_of_p(p);
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();

	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
	
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4711));

  DFeature<3>::QualityType q2;
	q2 = copy_of_p.getOverallQuality();
	TEST_REAL_EQUAL(q2, 0.9)
	q2 = copy_of_p.getQuality(0);
	TEST_REAL_EQUAL(q2, 0.1)
	q2 = copy_of_p.getQuality(1);
	TEST_REAL_EQUAL(q2, 0.2)
	q2 = copy_of_p.getQuality(2);
	TEST_REAL_EQUAL(q2, 0.3)
	TEST_EQUAL(copy_of_p.getModelDescription().getName(), "gauss")
RESULT

CHECK(DFeature& operator = (const DFeature& rhs))
	DFeature<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DFeature<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
  p.getOverallQuality() = 0.9;
  p.setQuality(0, 0.1);
  p.setQuality(1, 0.2);
  p.setQuality(2, 0.3);
  ModelDescription<3> desc;
  desc.setName("gauss");
  p.setModelDescription(desc);
	p.setMetaValue("cluster_id",4712);
	
	DFeature<3>::PositionType pos2;
	DFeature<3>::IntensityType i2;

	DFeature<3> copy_of_p;
	copy_of_p = p;
		
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();

  DFeature<3>::QualityType q2;

	TEST_REAL_EQUAL(i2, 123.456)
	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
	q2 = copy_of_p.getOverallQuality();
	TEST_REAL_EQUAL(q2, 0.9)
	q2 = copy_of_p.getQuality(0);
	TEST_REAL_EQUAL(q2, 0.1)
	q2 = copy_of_p.getQuality(1);
	TEST_REAL_EQUAL(q2, 0.2)
	q2 = copy_of_p.getQuality(2);
	TEST_REAL_EQUAL(q2, 0.3)
	TEST_EQUAL(copy_of_p.getModelDescription().getName(), "gauss")
RESULT

CHECK(DFeature& operator == (const DFeature& rhs))
  ModelDescription<1> desc;
  desc.setName("gauss");

	DFeature<1> p1;
	DFeature<1> p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getIntensity()=5;
  p1.getOverallQuality() = 0.9;
  p1.setQuality(0, 0.1);
  p1.setModelDescription(desc);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getIntensity()=5;
  p2.getOverallQuality() = 0.9;
  p2.setQuality(0, 0.1);
  p2.setModelDescription(desc);
	TEST_REAL_EQUAL(p1==p2, true)
		
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)	
RESULT

CHECK(DFeature& operator != (const DFeature& rhs))
	DFeature<1> p1;
	DFeature<1> p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
		
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)	
RESULT

CHECK(meta info with copy constructor)
	DFeature<1> p;
	p.setMetaValue(2,std::string("bla"));
 	DFeature<1> p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(meta info with assignment)
	DFeature<1> p;
	p.setMetaValue(2,std::string("bla"));
 	DFeature<1> p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK((const std::vector<Identification>& getIdentifications() const))
	DFeature<1> tmp;
	vector<Identification> vec(tmp.getIdentifications());
	TEST_EQUAL(vec.size(),0);
RESULT

CHECK((void setIdentifications(const std::vector<Identification>& identifications)))
	DFeature<1> tmp;
	vector<Identification> vec;
	
	tmp.setIdentifications(vec);
	TEST_EQUAL(tmp.getIdentifications().size(),0);
	
	Identification dbs;
	vec.push_back(dbs);
	tmp.setIdentifications(vec);
	TEST_EQUAL(tmp.getIdentifications().size(),1);
RESULT

CHECK((std::vector<Identification>& getIdentifications()))
	DFeature<1> tmp;
	vector<Identification> vec;
	
	tmp.getIdentifications().resize(1);
	TEST_EQUAL(tmp.getIdentifications().size(),1);
RESULT

//do not change these datastructures, they are used in the following tests...
std::vector< DConvexHull<2> > hulls(2);
hulls[0].addPoint(DPosition<2>(1.0,2.0));
hulls[0].addPoint(DPosition<2>(3.0,4.0));
hulls[1].addPoint(DPosition<2>(0.5,0.0));
hulls[1].addPoint(DPosition<2>(1.0,1.0));

CHECK(const ConvexHullVector& getConvexHulls() const)
	DFeature<2> tmp;
	TEST_EQUAL(tmp.getConvexHulls().size(),0)
RESULT

CHECK(ConvexHullVector& getConvexHulls())
	DFeature<2> tmp;
	tmp.getConvexHulls() = hulls;
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

CHECK(void setConvexHulls(const ConvexHullVector& hulls))
	DFeature<2> tmp;
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

CHECK((DBoundingBox<D> DFeature<D, Traits>::getBoundingBox() const))
	DFeature<2> tmp;
	tmp.setConvexHulls(hulls);
	DBoundingBox<2> bb = tmp.getBoundingBox();
	TEST_REAL_EQUAL(bb.min()[0],0.5)
	TEST_REAL_EQUAL(bb.min()[1],0.0)
	TEST_REAL_EQUAL(bb.max()[0],3.0)
	TEST_REAL_EQUAL(bb.max()[1],4.0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
