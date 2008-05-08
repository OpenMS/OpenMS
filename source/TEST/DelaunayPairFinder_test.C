// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition<2> PositionType;

START_TEST(DelaunayPairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DelaunayPairFinder<FeatureMap<> >* ptr = 0;
CHECK((DelaunayPairFinder()))
	ptr = new DelaunayPairFinder<FeatureMap<> >();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~DelaunayPairFinder()))
	delete ptr;
RESULT

DelaunayPairFinder<FeatureMap<> >::Point* ptr2 = 0;
CHECK(([EXTRA]Point()))
	ptr2 = new DelaunayPairFinder<FeatureMap<> >::Point();
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK(([EXTRA]~Point()))
	delete ptr2;
RESULT

CHECK(([EXTRA]Point& operator = (const Point& source)))
  DelaunayPairFinder<FeatureMap<> >::Point p(1,2);
  DelaunayPairFinder<FeatureMap<> >::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(Base::RT hx, Base::RT hy)))
  DelaunayPairFinder<FeatureMap<> >::Point p(1,2);
  DelaunayPairFinder<FeatureMap<> >::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(Base::RT hx, Base::RT hy, const PointType& f, UInt k=0)))
  DelaunayPairFinder<FeatureMap<> >::Point p(1,2);
  DelaunayPairFinder<FeatureMap<> >::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(const Base& cgal_point)))
  CGAL::Point_2< CGAL::Cartesian<double> > cp(1,2);
  DelaunayPairFinder<FeatureMap<> >::Point p(cp);
    	
  TEST_REAL_EQUAL(p.hx(),1)
  TEST_REAL_EQUAL(p.hy(),2)
RESULT

CHECK(([EXTRA]Point(const Point& source)))
  DelaunayPairFinder<FeatureMap<> >::Point p(1,2);
  DelaunayPairFinder<FeatureMap<> >::Point p_copy(p);
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point_2 operator()(const Circle_2& c) const))
  //
RESULT

CHECK((double getDiffIntercept(UInt dim)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  
  TEST_REAL_EQUAL(dpf.getDiffIntercept(0),1)
  TEST_REAL_EQUAL(dpf.getDiffIntercept(1),0.01)
RESULT

CHECK((float getMaxPairDistance(UInt dim)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  
  TEST_REAL_EQUAL(dpf.getMaxPairDistance(0),20)
  TEST_REAL_EQUAL(dpf.getMaxPairDistance(1),1)
RESULT

CHECK((float getPrecision(UInt dim)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  
  TEST_REAL_EQUAL(dpf.getPrecision(0),60)
  TEST_REAL_EQUAL(dpf.getPrecision(1),.5)
RESULT

CHECK((static BasePairFinder<PointMapType>* create()))
  // 
RESULT

CHECK((static const String getProductName()))
	DelaunayPairFinder<FeatureMap<> > dpf;
	
  TEST_EQUAL(dpf.getName() == "DelaunayPairFinder",true)
RESULT

CHECK((void findElementPairs()))
  FeatureMap<>scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  feat3.setPosition(pos3);
  feat3.setIntensity(400);
  scene.push_back(feat1);
  scene.push_back(feat2);
  scene.push_back(feat3);
  
  FeatureMap<>model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  feat5.setPosition(pos5);
  feat5.setIntensity(300);
  feat6.setPosition(pos6);
  feat6.setIntensity(400);
  model.push_back(feat4);
  model.push_back(feat5);
  model.push_back(feat6);
  
  DelaunayPairFinder< FeatureMap<> > dpf;
  dpf.setDiffIntercept(0,1.0);
  dpf.setDiffIntercept(1,1.0);
  dpf.setPrecision(0,5.0);
  dpf.setPrecision(1,5.0);
  dpf.setElementMap(0,model);
  dpf.setElementMap(1,scene);
  vector< ElementPair < Feature > > pairs;
  dpf.setElementPairs(pairs);
  dpf.findElementPairs();
  
  TEST_EQUAL((pairs.begin())->first == feat1, true)
  TEST_EQUAL((pairs.begin())->second == feat4, true)
  TEST_EQUAL((pairs.begin()+1)->first == feat2, true)
  TEST_EQUAL((pairs.begin()+1)->second == feat5, true)
  TEST_EQUAL((pairs.begin()+2)->first == feat3,true)
  TEST_EQUAL((pairs.begin()+2)->second == feat6,true)
RESULT

CHECK((template< typename ResultMapType > void computeConsensusMap(const PointMapType& first_map, ResultMapType& second_map)))
  ConsensusMap<ConsensusFeature<FeatureMap<> > > scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  feat3.setPosition(pos3);
  feat3.setIntensity(400); 
  ConsensusFeature<FeatureMap<> > cons1(0,0,feat1);
  ConsensusFeature<FeatureMap<> > cons2(0,1,feat2);
  ConsensusFeature<FeatureMap<> > cons3(0,2,feat3);
  scene.push_back(cons1);
  scene.push_back(cons2);
  scene.push_back(cons3);
  
  ConsensusMap<ConsensusFeature<FeatureMap<> > > model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  feat5.setPosition(pos5);
  feat5.setIntensity(300);
  feat6.setPosition(pos6);
  feat6.setIntensity(400);
  ConsensusFeature<FeatureMap<> > cons4(1,0,feat4);
  ConsensusFeature<FeatureMap<> > cons5(1,1,feat5);
  ConsensusFeature<FeatureMap<> > cons6(1,2,feat6);
  model.push_back(cons4);
  model.push_back(cons5);
  model.push_back(cons6);
  
  DelaunayPairFinder<ConsensusMap<ConsensusFeature<FeatureMap<> > > > dpf;
  dpf.setDiffIntercept(0,1.0);
  dpf.setDiffIntercept(1,1.0);
  dpf.setPrecision(0,5.0);
  dpf.setPrecision(1,5.0);
  dpf.computeConsensusMap(scene,model);
  Group group1 = model.begin()->getFeatures();
  Group group2 = (model.begin()+1)->getFeatures();
  Group group3 = (model.begin()+2)->getFeatures();
  
  IndexTuple ind1(0,0,feat1.getIntensity(),feat1.getPosition());
  IndexTuple ind2(0,1,feat2.getIntensity(),feat2.getPosition());
  IndexTuple ind3(0,2,feat3.getIntensity(),feat3.getPosition());
  IndexTuple ind4(1,0,feat4.getIntensity(),feat4.getPosition());
  IndexTuple ind5(1,1,feat5.getIntensity(),feat5.getPosition());
  IndexTuple ind6(1,2,feat6.getIntensity(),feat6.getPosition());

  Group::const_iterator it = group1.begin();
  TEST_EQUAL(*(it) == ind1, true)
	++it;
  TEST_EQUAL(*(it) == ind4, true)
	it = group2.begin();
  TEST_EQUAL(*(it) == ind2, true)
	++it;
  TEST_EQUAL(*(it) == ind5, true)
  it = group3.begin();
  TEST_EQUAL(*(it) == ind3, true)
	++it;
  TEST_EQUAL(*(it) == ind6, true)
RESULT

CHECK((void setDiffIntercept(UInt dim, DoubleReal intercept)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  dpf.setDiffIntercept(0,2);
  dpf.setDiffIntercept(1,2);
  
  TEST_REAL_EQUAL(dpf.getDiffIntercept(0),2)
  TEST_REAL_EQUAL(dpf.getDiffIntercept(1),2)
RESULT

CHECK((void setMaxPairDistance(UInt dim, Real max_pair_distance)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  dpf.setMaxPairDistance(0,2);
  dpf.setMaxPairDistance(1,2);
  
  TEST_REAL_EQUAL(dpf.getMaxPairDistance(0),2)
  TEST_REAL_EQUAL(dpf.getMaxPairDistance(1),2)
RESULT

CHECK((void setPrecision(UInt dim, Real precision)))
  DelaunayPairFinder<FeatureMap<> > dpf;
  dpf.setPrecision(0,2);
  dpf.setPrecision(1,2);
  
  TEST_REAL_EQUAL(dpf.getPrecision(0),2)
  TEST_REAL_EQUAL(dpf.getPrecision(1),2)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



